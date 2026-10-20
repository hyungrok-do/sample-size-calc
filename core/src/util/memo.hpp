// Values computed once per key and shared between threads and runs: the first thread to ask for a key
// computes it, and others asking meanwhile wait for that result instead of computing it again.  For
// quantities that depend on a few inputs but are needed for every scenario and every point of the power
// curve (a group-sequential design, the win probabilities of a hierarchical endpoint).
#pragma once

#include <cstddef>
#include <future>
#include <map>
#include <memory>
#include <mutex>

#include "util/parallel.hpp"

namespace ssc {

template <class Key, class Value>
class OnceCache {
 public:
  using Ptr = std::shared_ptr<const Value>;

  /// At most `capacity` values are kept (all are dropped when it is reached).
  explicit OnceCache(std::size_t capacity = 64) : capacity_(capacity) {}

  /// The value of `key`, from compute() the first time.  An exception of compute() reaches every thread
  /// waiting for it, and the key is computed again next time; a computation stopped because its own run
  /// was cancelled is redone for a waiting run that is not (`still_wanted` throws Cancelled if it is too).
  template <class Compute, class StillWanted>
  Ptr get(const Key& key, Compute&& compute, StillWanted&& still_wanted) {
    while (true) {
      std::promise<Ptr> promise;
      std::shared_future<Ptr> fut;
      bool owner = false;
      {
        std::lock_guard<std::mutex> lock(mu_);
        auto it = cache_.find(key);
        if (it != cache_.end()) {
          fut = it->second;
        } else {
          if (cache_.size() >= capacity_) cache_.clear();
          fut = promise.get_future().share();
          cache_.emplace(key, fut);
          owner = true;
        }
      }
      if (owner) {
        try {
          promise.set_value(std::make_shared<const Value>(compute()));
        } catch (...) {
          {
            std::lock_guard<std::mutex> lock(mu_);
            cache_.erase(key);
          }
          promise.set_exception(std::current_exception());
        }
      }
      try {
        return fut.get();
      } catch (const Cancelled&) {
        still_wanted();
      }
    }
  }

  template <class Compute>
  Ptr get(const Key& key, Compute&& compute) {
    return get(key, std::forward<Compute>(compute), [] { throw Cancelled(); });
  }

 private:
  std::mutex mu_;
  std::map<Key, std::shared_future<Ptr>> cache_;
  std::size_t capacity_;
};

}  // namespace ssc
