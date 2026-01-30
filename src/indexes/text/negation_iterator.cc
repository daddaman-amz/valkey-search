/*
 * Copyright (c) 2025, valkey-search contributors
 * All rights reserved.
 * SPDX-License-Identifier: BSD 3-Clause
 */

#include "src/indexes/text/negation_iterator.h"
#include "absl/log/check.h"
#include "vmsdk/src/log.h"

namespace valkey_search::indexes::text {

NegationTextIterator::NegationTextIterator(
    InternedStringSet tracked_keys,
    InternedStringSet matched_keys,
    InternedStringSet untracked_keys,
    FieldMaskPredicate field_mask)
    : tracked_keys_(std::move(tracked_keys)),
      matched_keys_(std::move(matched_keys)),
      untracked_keys_(std::move(untracked_keys)),
      field_mask_(field_mask),
      in_untracked_phase_(false) {
  current_iter_ = tracked_keys_.begin();
  tracked_end_ = tracked_keys_.end();
  untracked_begin_ = untracked_keys_.begin();
  untracked_end_ = untracked_keys_.end();
  AdvanceToNextValid();
}

FieldMaskPredicate NegationTextIterator::QueryFieldMask() const {
  VMSDK_LOG(NOTICE, nullptr) << " Inside NegationTextIterator::QueryFieldMask() >> field_mask_ is " << field_mask_;
  return field_mask_;
}

bool NegationTextIterator::DoneKeys() const {
  return !current_key_;
}


const Key& NegationTextIterator::CurrentKey() const {
  // CHECK(current_key_) << "CurrentKey called on invalid iterator";
  return current_key_;
}

void NegationTextIterator::AdvanceToNextValid() {
  if (!in_untracked_phase_) {
    while (current_iter_ != tracked_end_) {
      if (!matched_keys_.contains(*current_iter_)) {
        current_key_ = *current_iter_;
        return;
      }
      ++current_iter_;
    }
    in_untracked_phase_ = true;
    current_iter_ = untracked_begin_;
  }
  
  if (current_iter_ != untracked_end_) {
    current_key_ = *current_iter_;
  } else {
    current_key_ = Key();  // Reset to empty/null key
  }
}

bool NegationTextIterator::NextKey() {
  if (current_key_) {
    ++current_iter_;
    AdvanceToNextValid();
  }
  return static_cast<bool>(current_key_);
}

bool NegationTextIterator::SeekForwardKey(const Key& target_key) {
  if (current_key_ && current_key_ >= target_key) {
    return true;
  }
  while (current_key_ && current_key_ < target_key) {
    if (!NextKey()) {
      return false;
    }
  }
  return static_cast<bool>(current_key_);
}

bool NegationTextIterator::DonePositions() const {
  return true;
}

const PositionRange& NegationTextIterator::CurrentPosition() const {
  CHECK(false) << "Negation iterator has no positions";
  static PositionRange dummy{0, 0};
  return dummy;
}

bool NegationTextIterator::NextPosition() {
  return false;
}

bool NegationTextIterator::SeekForwardPosition(Position target_position) {
  return false;
}

FieldMaskPredicate NegationTextIterator::CurrentFieldMask() const {
  return field_mask_;
}

bool NegationTextIterator::IsIteratorValid() const {
  return static_cast<bool>(current_key_);
}

}  // namespace valkey_search::indexes::text
