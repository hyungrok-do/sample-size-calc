#include "parallel.hpp"

#include <algorithm>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>
#if !defined(_WIN32)
#include <pthread.h>
#endif

namespace ssc {

namespace {
thread_local bool tl_in_parallel = false;

class Pool;

// One parallel_for in progress.  The thread that called it works on its items, and so do the pool's workers
// that join it (up to the threads the call asked for, less the caller).
struct Job {
  std::size_t n = 0, chunk = 1;  // items, and how many a thread takes at a time
  const std::function<void(std::size_t)>* fn = nullptr;
  const Token* token = nullptr;
  int max_helpers = 0;
  int helpers = 0;               // workers on it now (under the pool's mutex)
  std::exception_ptr error;      // the first exception of an item (under the pool's mutex)
  std::atomic<std::size_t> next{0};
  std::atomic<bool> stop{false};
  std::condition_variable left;  // the last worker on it has left
};

// Worker threads started when first needed and then kept, waiting for work, for the life of the process.  A
// worker joins the newest job that takes another helper: the work nested in a job (a scenario's simulation)
// comes before the job's next item, and nested work is shared with idle workers instead of running alone.
// The calling threads never wait for anything but their own job's items, so a nested job cannot deadlock.
class Pool {
 public:
  std::mutex mutex;
  std::condition_variable work;  // a job was published
  std::vector<Job*> jobs;        // published jobs, oldest first

  // Start workers until there are `count` (under the mutex).
  void ensure(int count) {
    while (workers_ < count) {
      try {
        std::thread([this] { loop(); }).detach();
      } catch (...) {
        return;  // no more threads: the caller and the workers there are do the work
      }
      ++workers_;
    }
  }

  void run(Job& job) {
    const bool was = tl_in_parallel;
    tl_in_parallel = true;
    while (!job.stop.load(std::memory_order_relaxed)) {
      const std::size_t first = job.next.fetch_add(job.chunk, std::memory_order_relaxed);
      if (first >= job.n) break;
      const std::size_t last = std::min(job.n, first + job.chunk);
      try {
        for (std::size_t i = first; i < last; ++i) {
          if (job.token && job.token->cancelled()) {
            job.stop.store(true);
            break;
          }
          (*job.fn)(i);
        }
      } catch (...) {
        std::lock_guard<std::mutex> lock(mutex);
        if (!job.error) job.error = std::current_exception();
        job.stop.store(true);
      }
    }
    tl_in_parallel = was;
  }

 private:
  int workers_ = 0;

  Job* pick() {
    for (auto it = jobs.rbegin(); it != jobs.rend(); ++it) {
      Job* j = *it;
      if (j->helpers < j->max_helpers && !j->stop.load(std::memory_order_relaxed) &&
          j->next.load(std::memory_order_relaxed) < j->n)
        return j;
    }
    return nullptr;
  }

  void loop() {
    std::unique_lock<std::mutex> lock(mutex);
    while (true) {
      Job* job = pick();
      if (!job) {
        work.wait(lock);
        continue;
      }
      ++job->helpers;
      lock.unlock();
      run(*job);
      lock.lock();
      if (--job->helpers == 0) job->left.notify_all();
    }
  }
};

// The pool is never destroyed: its workers may still be waiting on it while the process exits.  A child
// process made by fork() has none of the parent's workers, and starts a pool of its own.
Pool* g_pool = nullptr;
std::once_flag g_pool_once;

Pool& pool() {
  std::call_once(g_pool_once, [] {
    g_pool = new Pool();
#if !defined(_WIN32)
    pthread_atfork(nullptr, nullptr, [] { g_pool = new Pool(); });
#endif
  });
  return *g_pool;
}
}  // namespace

int hardware_threads() {
  const unsigned n = std::thread::hardware_concurrency();
  return n == 0 ? 1 : static_cast<int>(n);
}

int resolve_threads(int requested) {
  const int hw = hardware_threads();
  if (requested <= 0) return hw;
  return std::min(requested, 4 * hw);
}

bool in_parallel_region() { return tl_in_parallel; }

void parallel_for(std::size_t n, int threads, const std::function<void(std::size_t)>& fn, const Token* token) {
  if (n == 0) return;
  const int k = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(resolve_threads(threads)), n));
  if (k <= 1) {
    for (std::size_t i = 0; i < n; ++i) {
      if (token) token->check();
      fn(i);
    }
    return;
  }

  Pool& p = pool();
  Job job;
  job.n = n;
  job.fn = &fn;
  job.token = token;
  job.max_helpers = k - 1;
  // a few dozen chunks per thread: dynamic enough to balance uneven items, few enough to keep the shared
  // counter quiet when items are tiny
  job.chunk = std::max<std::size_t>(1, n / (32 * static_cast<std::size_t>(k)));
  {
    std::lock_guard<std::mutex> lock(p.mutex);
    p.ensure(k - 1);
    p.jobs.push_back(&job);
  }
  p.work.notify_all();
  p.run(job);
  {
    std::unique_lock<std::mutex> lock(p.mutex);
    p.jobs.erase(std::find(p.jobs.begin(), p.jobs.end(), &job));
    job.left.wait(lock, [&] { return job.helpers == 0; });
  }
  if (job.error) std::rethrow_exception(job.error);
  if (token) token->check();
}

}  // namespace ssc
