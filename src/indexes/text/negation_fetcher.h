#ifndef VALKEYSEARCH_SRC_INDEXES_TEXT_NEGATION_FETCHER_H_
#define VALKEYSEARCH_SRC_INDEXES_TEXT_NEGATION_FETCHER_H_

#include <memory>
#include "src/indexes/index_base.h"
#include "src/indexes/text/text_iterator.h"
#include "src/indexes/text/text_fetcher.h"

namespace valkey_search::query {

class NegationFetcher : public indexes::EntriesFetcherBase {
 public:
  NegationFetcher(std::unique_ptr<indexes::text::TextIterator> iter, size_t size)
      : iter_(std::move(iter)), size_(size) {}
  
  size_t Size() const override { return size_; }
  
  std::unique_ptr<indexes::EntriesFetcherIteratorBase> Begin() override {
    return std::make_unique<indexes::text::TextFetcher>(std::move(iter_));
  }

 private:
  std::unique_ptr<indexes::text::TextIterator> iter_;
  size_t size_;
};

}  // namespace valkey_search::query

#endif  // VALKEYSEARCH_SRC_INDEXES_TEXT_NEGATION_FETCHER_H_
