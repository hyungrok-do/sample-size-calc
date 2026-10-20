// Parallel execution utilities.
//
// parallel_for() distributes independent work items (scenarios, curve points,
// simulation replicates, ...) over threads with dynamic scheduling: the
// calling thread and the workers of a pool started once and kept.  A nested
// call (a scenario's simulation inside the loop over scenarios) is shared with
// the workers that are idle, so a few scenarios still use every core, and the
// machine is never oversubscribed.  A Token allows cooperative cancellation
// and progress reporting across threads.
#pragma once

#include <atomic>
#include <cstddef>
#include <functional>
#include <stdexcept>

namespace ssc {

/// Thrown by Token::check() once cancellation has been requested.
struct Cancelled : std::runtime_error {
  Cancelled() : std::runtime_error("calculation cancelled") {}
};

class Token {
 public:
  void cancel() { cancelled_.store(true, std::memory_order_relaxed); }
  bool cancelled() const { return cancelled_.load(std::memory_order_relaxed); }
  void check() const {
    if (cancelled()) throw Cancelled();
  }
  void set_total(long n) {
    total_.store(n > 0 ? n : 1, std::memory_order_relaxed);
    done_.store(0, std::memory_order_relaxed);
  }
  void add_done(long n = 1) { done_.fetch_add(n, std::memory_order_relaxed); }
  double progress() const {
    const double t = static_cast<double>(total_.load(std::memory_order_relaxed));
    const double d = static_cast<double>(done_.load(std::memory_order_relaxed));
    return d >= t ? 1.0 : d / t;
  }

 private:
  std::atomic<bool> cancelled_{false};
  std::atomic<long> total_{1};
  std::atomic<long> done_{0};
};

/// Number of hardware threads (at least 1).
int hardware_threads();

/// Resolve a requested thread count: <= 0 means "all hardware threads".
int resolve_threads(int requested);

/// Run fn(i) for i in [0, n) on up to `threads` threads (the calling thread
/// participates; nested calls get the workers that are free).  The first
/// exception thrown by any task is rethrown after all workers have stopped.
/// If `token` is cancelled, remaining items are skipped and Cancelled is thrown.
void parallel_for(std::size_t n, int threads, const std::function<void(std::size_t)>& fn,
                  const Token* token = nullptr);

/// True when called from inside a parallel_for task.
bool in_parallel_region();

}  // namespace ssc
