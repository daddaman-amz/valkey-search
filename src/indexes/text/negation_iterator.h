/*
 * Copyright (c) 2025, valkey-search contributors
 * All rights reserved.
 * SPDX-License-Identifier: BSD 3-Clause
 */

#ifndef VALKEY_SEARCH_INDEXES_TEXT_NEGATION_ITERATOR_H_
#define VALKEY_SEARCH_INDEXES_TEXT_NEGATION_ITERATOR_H_

#include "src/indexes/text/posting.h"
#include "src/indexes/text/text_iterator.h"
#include "src/utils/string_interning.h"

namespace valkey_search::indexes::text {

class NegationTextIterator : public TextIterator {
 public:
  NegationTextIterator(InternedStringSet tracked_keys,
                       InternedStringSet matched_keys,
                       InternedStringSet untracked_keys,
                       FieldMaskPredicate field_mask);

  FieldMaskPredicate QueryFieldMask() const override;
  bool DoneKeys() const override;
  const Key& CurrentKey() const override;
  bool NextKey() override;
  bool SeekForwardKey(const Key& target_key) override;
  bool DonePositions() const override;
  const PositionRange& CurrentPosition() const override;
  bool NextPosition() override;
  bool SeekForwardPosition(Position target_position) override;
  FieldMaskPredicate CurrentFieldMask() const override;
  bool IsIteratorValid() const override;

 private:
  InternedStringSet tracked_keys_;
  InternedStringSet matched_keys_;
  InternedStringSet untracked_keys_;
  FieldMaskPredicate field_mask_;
  
  InternedStringSet::const_iterator current_iter_;
  InternedStringSet::const_iterator tracked_end_;
  InternedStringSet::const_iterator untracked_begin_;
  InternedStringSet::const_iterator untracked_end_;
  
  Key current_key_;
  bool in_untracked_phase_;
  
  void AdvanceToNextValid();
};

}  // namespace valkey_search::indexes::text

#endif
