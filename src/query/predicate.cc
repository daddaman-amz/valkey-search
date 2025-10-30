/*
 * Copyright (c) 2025, valkey-search contributors
 * All rights reserved.
 * SPDX-License-Identifier: BSD 3-Clause
 *
 */

#include "src/query/predicate.h"

#include <memory>
#include <string>
#include <utility>

#include "absl/container/flat_hash_set.h"
#include "absl/strings/match.h"
#include "absl/strings/string_view.h"
#include "src/indexes/numeric.h"
#include "src/indexes/tag.h"
#include "vmsdk/src/log.h"
#include "vmsdk/src/managed_pointers.h"

#include "src/indexes/text/text_index.h"
#include "src/indexes/text/text_iterator.h"
#include "src/indexes/text/proximity.h"
#include "src/indexes/text.h"

namespace valkey_search::query {


namespace {
  constexpr bool kEnablePredicateEvaluatorProximity = true;  // Set to false to disable

  // Returns which term violated constraints, or nullopt if valid
  std::optional<size_t> FindViolatingTerm(
      const std::vector<valkey_search::indexes::text::PositionRange>& positions,
      uint32_t slop, bool inorder) {
    
    const size_t n = positions.size();
    
    if (inorder) {
      for (size_t i = 0; i < n - 1; ++i) {
        // Check overlap/ordering violations
        if (positions[i].end >= positions[i + 1].start) {
          return i + 1;  // Advance next term
        }
        // Check slop violations
        if (positions[i + 1].start - positions[i].end - 1 > slop) {
          return i;  // Advance current term
        }
      }
      return std::nullopt;  // Valid!
    }
    
    // Unordered: sort and check
    std::vector<std::pair<uint32_t, size_t>> sorted_positions;
    for (size_t i = 0; i < n; ++i) {
      sorted_positions.push_back({positions[i].start, i});
    }
    std::sort(sorted_positions.begin(), sorted_positions.end());
    
    for (size_t i = 0; i < n - 1; ++i) {
      size_t curr_idx = sorted_positions[i].second;
      size_t next_idx = sorted_positions[i + 1].second;
      
      if (positions[curr_idx].end >= positions[next_idx].start) {
        return next_idx;
      }
      if (positions[next_idx].start - positions[curr_idx].end - 1 > slop) {
        return curr_idx;
      }
    }
    return std::nullopt;
  }


    // Returns all position ranges that satisfy proximity constraints
  // CHANGE: New function to return matching ranges instead of just bool
  std::vector<valkey_search::indexes::text::PositionRange> FindMatchingProximityRanges(
      const std::vector<std::vector<valkey_search::indexes::text::PositionRange>>& term_positions,
      uint32_t slop, bool inorder) {
    
    std::vector<valkey_search::indexes::text::PositionRange> matching_ranges;
    
    if (term_positions.empty()) return matching_ranges;
    
    for (const auto& positions : term_positions) {
      if (positions.empty()) return matching_ranges;
    }
    
    const size_t num_terms = term_positions.size();
    std::vector<size_t> indices(num_terms, 0);
    std::vector<valkey_search::indexes::text::PositionRange> current_combo(num_terms);
    bool should_advance = false;
    
    // DEBUG: Count iterations
    size_t iteration_count = 0;
    size_t total_combinations = 1;
    for (const auto& positions : term_positions) {
      total_combinations *= positions.size();
    }
    // DEBUG: End counting setup

    while (true) {
      iteration_count++; // DEBUG
      
      // Build current combination
      for (size_t i = 0; i < num_terms; ++i) {
        current_combo[i] = term_positions[i][indices[i]];
      }
      
      if (should_advance) {
        should_advance = false;
        
        auto violating_term = FindViolatingTerm(current_combo, slop, inorder);
        
        if (violating_term.has_value()) {
          size_t advance_idx = violating_term.value();
          indices[advance_idx]++;
          if (indices[advance_idx] >= term_positions[advance_idx].size()) {
            // DEBUG: Log on completion
            VMSDK_LOG(WARNING, nullptr) << "FindMatchingProximityRanges: found " 
                                        << matching_ranges.size() << " matches after testing " 
                                        << iteration_count << " combinations out of " 
                                        << total_combinations << " possible";
            // DEBUG: End log
            return matching_ranges; // CHANGE: Return collected ranges
          }
        } else {
          size_t advance_idx = 0;
          while (advance_idx < num_terms) {
            indices[advance_idx]++;
            if (indices[advance_idx] < term_positions[advance_idx].size()) {
              break;
            }
            indices[advance_idx] = 0;
            advance_idx++;
          }
          if (advance_idx >= num_terms) {
            // DEBUG: Log on completion
            VMSDK_LOG(WARNING, nullptr) << "FindMatchingProximityRanges: found " 
                                        << matching_ranges.size() << " matches after testing " 
                                        << iteration_count << " combinations out of " 
                                        << total_combinations << " possible";
            // DEBUG: End log
            return matching_ranges; // CHANGE: Return collected ranges
          }
        }
        continue;
      }
      
      if (!FindViolatingTerm(current_combo, slop, inorder).has_value()) {
        // CHANGE: Valid combination found - compute the span and add to results
        uint32_t min_start = current_combo[0].start;
        uint32_t max_end = current_combo[0].end;
        for (const auto& pos : current_combo) {
          min_start = std::min(min_start, pos.start);
          max_end = std::max(max_end, pos.end);
        }
        matching_ranges.push_back({min_start, max_end});
        // CHANGE: End of span computation
        
        // Continue searching for more matches
        should_advance = true;
      } else {
        should_advance = true;
      }
    }
  }

  // Formatting for positions
  std::string FormatPositionRanges(const std::vector<valkey_search::indexes::text::PositionRange>& ranges) {
    if (ranges.empty()) return "[]";
    
    std::string result = "[";
    for (size_t i = 0; i < ranges.size(); ++i) {
      if (i > 0) result += ", ";
      result += "(start=" + std::to_string(ranges[i].start) + 
                " end=" + std::to_string(ranges[i].end) + ")";
    }
    result += "]";
    return result;
  }
}  

EvaluationResult NegatePredicate::Evaluate(Evaluator& evaluator) const {
  EvaluationResult result = predicate_->Evaluate(evaluator);
  return EvaluationResult(!result.matches);
}

TermPredicate::TermPredicate(const indexes::Text* index,
                             absl::string_view identifier,
                             absl::string_view alias, std::string term)
    : TextPredicate(),
      index_(index),
      identifier_(vmsdk::MakeUniqueValkeyString(identifier)),
      alias_(alias),
      term_(term) {}

EvaluationResult TermPredicate::Evaluate(Evaluator& evaluator) const {
  // call dynamic dispatch on the evaluator
  return evaluator.EvaluateText(*this);
}

EvaluationResult TermPredicate::Evaluate(const valkey_search::indexes::text::TextIndex& text_index) const {
  // Compute field_mask from own index (ignore parameter for now)
  size_t field_number = index_->GetTextFieldNumber();
  uint64_t field_mask = 1ULL << field_number;

  VMSDK_LOG(WARNING, nullptr) << "TermPredicate::Evaluate with index - field: '" 
                              << vmsdk::ToStringView(identifier_.get()) 
                              << "' search term: '" << term_ << "'"
                              << "' computed field mask: " << field_mask;
  
  auto word_iter = text_index.prefix_.GetWordIterator(term_);
  if (word_iter.Done()) {
    VMSDK_LOG(WARNING, nullptr) << "Word '" << term_ << "' not found in prefix tree";
    return EvaluationResult(false);
  }
  
  auto postings = word_iter.GetTarget();
  if (!postings) {
    VMSDK_LOG(WARNING, nullptr) << "No postings for word '" << term_ << "'";
    return EvaluationResult(false);
  }
  
  auto key_iter = postings->GetKeyIterator();
  bool contains_fields = key_iter.IsValid() && key_iter.ContainsFields(field_mask);
  if (!contains_fields) {
    VMSDK_LOG(WARNING, nullptr) << "Word '" << term_ << "' - no valid fields for the key";
  }
  
  VMSDK_LOG(WARNING, nullptr) << "Word '" << term_ << "' postings - key_count: " 
                              << postings->GetKeyCount() << " postings_count: " 
                              << postings->GetPostingCount() << " contains_fields: " 
                              << (contains_fields ? "TRUE" : "FALSE");
  // Extract positions
  std::vector<valkey_search::indexes::text::PositionRange> ranges;
  auto pos_iter = key_iter.GetPositionIterator();
  while (pos_iter.IsValid()) {
      // Check if this position belongs to the requested field
      if ((pos_iter.GetFieldMask() & field_mask) != 0) {
          auto pos = pos_iter.GetPosition();
          ranges.emplace_back(pos, pos);
      }
      pos_iter.NextPosition();
  }

  VMSDK_LOG(WARNING, nullptr) << "Word '" << term_ << "' found " << ranges.size() 
                            << " positions: " << FormatPositionRanges(ranges);
  if (ranges.empty()) {
    VMSDK_LOG(WARNING, nullptr) << "Word '" << term_ << "' found but no positions in requested field";
    return EvaluationResult(false);
  }
  return EvaluationResult(true, std::move(ranges));
}

PrefixPredicate::PrefixPredicate(const indexes::Text* index,
                                 absl::string_view identifier,
                                 absl::string_view alias, std::string term)
    : TextPredicate(),
      index_(index),
      identifier_(vmsdk::MakeUniqueValkeyString(identifier)),
      alias_(alias),
      term_(term) {}

EvaluationResult PrefixPredicate::Evaluate(Evaluator& evaluator) const {
  return evaluator.EvaluateText(*this);
}

EvaluationResult PrefixPredicate::Evaluate(const valkey_search::indexes::text::TextIndex& text_index) const {
  size_t field_number = index_->GetTextFieldNumber();
  uint64_t field_mask = 1ULL << field_number;
  
  VMSDK_LOG(WARNING, nullptr) << "PrefixPredicate::Evaluate with index - field: '" 
                              << vmsdk::ToStringView(identifier_.get()) 
                              << "' search term: '" << term_ << "'"
                              << " field_mask: " << field_mask;
  
  std::vector<valkey_search::indexes::text::PositionRange> all_ranges;
  auto word_iter = text_index.prefix_.GetWordIterator(term_);
  int matches_found = 0;
  
  while (!word_iter.Done()) {
    std::string_view word = word_iter.GetWord();
    if (!word.starts_with(term_)) break;
    
    auto postings = word_iter.GetTarget();
    if (postings) {
      auto key_iter = postings->GetKeyIterator();
      if (key_iter.IsValid() && key_iter.ContainsFields(field_mask)) {
        auto pos_iter = key_iter.GetPositionIterator();
        while (pos_iter.IsValid()) {
          // Check if this position belongs to the requested field
          if ((pos_iter.GetFieldMask() & field_mask) != 0) {
              auto pos = pos_iter.GetPosition();
              all_ranges.emplace_back(pos, pos);
          }
          pos_iter.NextPosition();
        }
        matches_found++;
        VMSDK_LOG(WARNING, nullptr) << "Prefix match: word '" << word << "' matches prefix '" << term_ << "'";
      }
    }
    word_iter.Next();
  }

  VMSDK_LOG(WARNING, nullptr) << "Prefix '" << term_ << "' found " << all_ranges.size() 
                            << " positions: " << FormatPositionRanges(all_ranges)
                            << " from " << matches_found << " words";
  return all_ranges.empty() ? EvaluationResult(false) : EvaluationResult(true, std::move(all_ranges));
}

SuffixPredicate::SuffixPredicate(const indexes::Text* index,
                                 absl::string_view identifier,
                                 absl::string_view alias, std::string term)
    : TextPredicate(),
      index_(index),
      identifier_(vmsdk::MakeUniqueValkeyString(identifier)),
      alias_(alias),
      term_(term) {}

EvaluationResult SuffixPredicate::Evaluate(Evaluator& evaluator) const {
  return evaluator.EvaluateText(*this);
}

EvaluationResult SuffixPredicate::Evaluate(const valkey_search::indexes::text::TextIndex& text_index) const {
  // TODO: Implement suffix logic with index
  VMSDK_LOG(WARNING, nullptr) << "SuffixPredicate::Evaluate >> Need impl return true for now '";
  return EvaluationResult(true); 
}

InfixPredicate::InfixPredicate(const indexes::Text* index,
                               absl::string_view identifier,
                               absl::string_view alias, std::string term)
    : TextPredicate(),
      index_(index),
      identifier_(vmsdk::MakeUniqueValkeyString(identifier)),
      alias_(alias),
      term_(term) {}

EvaluationResult InfixPredicate::Evaluate(Evaluator& evaluator) const {
  return evaluator.EvaluateText(*this);
}

EvaluationResult InfixPredicate::Evaluate(const valkey_search::indexes::text::TextIndex& text_index) const {
  // TODO: Implement infix logic with index
  VMSDK_LOG(WARNING, nullptr) << "InfixPredicate::Evaluate >> Need impl return true for now '";
  return EvaluationResult(true); 
}

FuzzyPredicate::FuzzyPredicate(const indexes::Text* index,
                               absl::string_view identifier,
                               absl::string_view alias, std::string term,
                               uint32_t distance)
    : TextPredicate(),
      index_(index),
      identifier_(vmsdk::MakeUniqueValkeyString(identifier)),
      alias_(alias),
      term_(term),
      distance_(distance) {}

EvaluationResult FuzzyPredicate::Evaluate(Evaluator& evaluator) const {
  return evaluator.EvaluateText(*this);
}

EvaluationResult FuzzyPredicate::Evaluate(const valkey_search::indexes::text::TextIndex& text_index) const {
  // TODO: Implement fuzzy logic with index
  VMSDK_LOG(WARNING, nullptr) << "FuzzyPredicate::Evaluate >> Need impl return true for now '";
  return EvaluationResult(true); 
}

ProximityPredicate::ProximityPredicate(
    std::vector<std::unique_ptr<TextPredicate>> terms, uint32_t slop,
    bool inorder)
    : TextPredicate(),
      terms_(std::move(terms)),
      inorder_(inorder),
      slop_(slop) {}

EvaluationResult ProximityPredicate::Evaluate(Evaluator& evaluator) const {
  return evaluator.EvaluateText(*this);
}
      
EvaluationResult ProximityPredicate::Evaluate(
    const valkey_search::indexes::text::TextIndex& text_index) const {
  
  if (!kEnablePredicateEvaluatorProximity) {
    VMSDK_LOG(WARNING, nullptr) << "ProximityPredicate::Evaluate - DISABLED";
    return EvaluationResult(true);
  }
  
  VMSDK_LOG(WARNING, nullptr) << "ProximityPredicate::Evaluate - checking " 
                              << terms_.size() << " terms with slop=" << slop_ 
                              << " inorder=" << inorder_;
  
  // Step 1: Evaluate each term and collect their positions
  std::vector<std::vector<valkey_search::indexes::text::PositionRange>> term_positions;
  term_positions.reserve(terms_.size());
  
  for (const auto& term : terms_) {
    auto result = term->Evaluate(text_index);
    if (!result.matches) {
      VMSDK_LOG(WARNING, nullptr) << "ProximityPredicate: term not found";
      return EvaluationResult(false);
    }
    term_positions.push_back(std::move(result.position_ranges));
  }

  // CHANGE START: Use new function and return position ranges
  // Step 2: Find valid position combinations and return matching ranges
  auto matching_ranges = FindMatchingProximityRanges(term_positions, slop_, inorder_);
  
  if (matching_ranges.empty()) {
    VMSDK_LOG(WARNING, nullptr) << "ProximityPredicate: proximity NOT satisfied";
    return EvaluationResult(false);
  }
  
  // LOG: Add this line to log the matching ranges
  VMSDK_LOG(WARNING, nullptr) << "ProximityPredicate: proximity SATISFIED with " 
                              << matching_ranges.size() << " matching ranges: "
                              << FormatPositionRanges(matching_ranges);
                              
  return EvaluationResult(true, std::move(matching_ranges));  // FIX: Return with ranges!

}

NumericPredicate::NumericPredicate(const indexes::Numeric* index,
                                   absl::string_view alias,
                                   absl::string_view identifier, double start,
                                   bool is_inclusive_start, double end,
                                   bool is_inclusive_end)
    : Predicate(PredicateType::kNumeric),
      index_(index),
      alias_(alias),
      identifier_(vmsdk::MakeUniqueValkeyString(identifier)),
      start_(start),
      is_inclusive_start_(is_inclusive_start),
      end_(end),
      is_inclusive_end_(is_inclusive_end) {}

EvaluationResult NumericPredicate::Evaluate(Evaluator& evaluator) const {
  return evaluator.EvaluateNumeric(*this);
}

EvaluationResult NumericPredicate::Evaluate(const double* value) const {
  if (!value) {
    return EvaluationResult(false);
  }
  bool matches = (((*value > start_ || (is_inclusive_start_ && *value == start_)) &&
           (*value < end_)) ||
          (is_inclusive_end_ && *value == end_));
  return EvaluationResult(matches);        
}

TagPredicate::TagPredicate(const indexes::Tag* index, absl::string_view alias,
                           absl::string_view identifier,
                           absl::string_view raw_tag_string,
                           const absl::flat_hash_set<absl::string_view>& tags)
    : Predicate(PredicateType::kTag),
      index_(index),
      alias_(alias),
      identifier_(vmsdk::MakeUniqueValkeyString(identifier)),
      raw_tag_string_(raw_tag_string),
      tags_(tags.begin(), tags.end()) {}

EvaluationResult TagPredicate::Evaluate(Evaluator& evaluator) const {
  return evaluator.EvaluateTags(*this);
}

EvaluationResult TagPredicate::Evaluate(
    const absl::flat_hash_set<absl::string_view>* in_tags,
    bool case_sensitive) const {
           VMSDK_LOG(WARNING, nullptr) << "In predicate.cc >> TagPredicate::Evaluate";  //msremove     
  if (!in_tags) {
    return EvaluationResult(false);
  }

  for (const auto& in_tag : *in_tags) {
    for (const auto& tag : tags_) {
      absl::string_view left_hand_side = in_tag;
      absl::string_view right_hand_side = tag;
      if (right_hand_side.back() == '*') {
        if (left_hand_side.length() < right_hand_side.length() - 1) {
          continue;
        }
        left_hand_side = left_hand_side.substr(0, right_hand_side.length() - 1);
        right_hand_side =
            right_hand_side.substr(0, right_hand_side.length() - 1);
      }
      if (case_sensitive) {
        if (left_hand_side == right_hand_side) {
          return EvaluationResult(true);
        }
      } else {
        if (absl::EqualsIgnoreCase(left_hand_side, right_hand_side)) {
          return EvaluationResult(true);
        }
      }
    }
  }
  return EvaluationResult(false);
}

ComposedPredicate::ComposedPredicate(std::unique_ptr<Predicate> lhs_predicate,
                                     std::unique_ptr<Predicate> rhs_predicate,
                                     LogicalOperator logical_op)
    : Predicate(logical_op == LogicalOperator::kAnd
                    ? PredicateType::kComposedAnd
                    : PredicateType::kComposedOr),
      lhs_predicate_(std::move(lhs_predicate)),
      rhs_predicate_(std::move(rhs_predicate)) {}

EvaluationResult ComposedPredicate::Evaluate(Evaluator& evaluator) const {
  if (GetType() == PredicateType::kComposedAnd) {
    EvaluationResult lhs = lhs_predicate_->Evaluate(evaluator);
    VMSDK_LOG(DEBUG, nullptr) << "Inline evaluate AND predicate lhs: " << lhs.matches;
    
    // Short-circuit for AND
    if (!lhs.matches) {
      return EvaluationResult(false);
    }
    
    EvaluationResult rhs = rhs_predicate_->Evaluate(evaluator);
    VMSDK_LOG(DEBUG, nullptr) << "Inline evaluate AND predicate rhs: " << rhs.matches;
    
    // TODO: Implement position-aware AND logic here for complex proximity
    // For now, just combine boolean results
    return EvaluationResult(lhs.matches && rhs.matches);
  }

  // OR logic
  EvaluationResult lhs = lhs_predicate_->Evaluate(evaluator);
  VMSDK_LOG(DEBUG, nullptr) << "Inline evaluate OR predicate lhs: " << lhs.matches;
  EvaluationResult rhs = rhs_predicate_->Evaluate(evaluator);
  VMSDK_LOG(DEBUG, nullptr) << "Inline evaluate OR predicate rhs: " << rhs.matches;
  
  // TODO: Implement position-aware OR logic here for complex proximity
  return EvaluationResult(lhs.matches || rhs.matches);
}

}  // namespace valkey_search::query
