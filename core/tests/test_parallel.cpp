// The thread pool behind parallel_for: every item once, nested loops shared with idle workers, exceptions and
// cancellation reaching the caller from any depth, and many callers at once.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <mutex>
#include <numeric>
#include <set>
#include <stdexcept>
#include <thread>
#include <vector>

#include "testing.hpp"
#include "util/parallel.hpp"

using namespace ssc;

TEST_CASE("parallel: every item exactly once, on any number of threads") {
  for (int threads : {0, 1, 2, 3, 8, 64})
    for (std::size_t n : {0, 1, 2, 7, 64, 1000, 20000}) {
      std::vector<std::atomic<int>> seen(n);
      parallel_for(n, threads, [&](std::size_t i) { seen[i].fetch_add(1); });
      CHECK(std::all_of(seen.begin(), seen.end(), [](const std::atomic<int>& s) { return s.load() == 1; }));
    }
}

TEST_CASE("parallel: a nested loop is shared with the idle workers, with the same results") {
  // two outer items on four threads: before the pool, each inner loop ran alone on its outer item's thread
  std::mutex mu;
  std::set<std::thread::id> inner_threads;
  std::vector<std::vector<double>> got(2, std::vector<double>(200));
  parallel_for(2, 4, [&](std::size_t o) {
    parallel_for(got[o].size(), 4, [&](std::size_t i) {
      std::this_thread::sleep_for(std::chrono::microseconds(200));  // long enough for idle workers to join
      got[o][i] = static_cast<double>(o * 1000 + i) * 0.5;
      std::lock_guard<std::mutex> lock(mu);
      inner_threads.insert(std::this_thread::get_id());
    });
  });
  for (std::size_t o = 0; o < 2; ++o)
    for (std::size_t i = 0; i < got[o].size(); ++i) CHECK_EQ(got[o][i], static_cast<double>(o * 1000 + i) * 0.5);
  CHECK(inner_threads.size() > 2);
  CHECK(!in_parallel_region());
}

TEST_CASE("parallel: an exception in a nested loop reaches the caller, and the pool goes on") {
  bool caught = false;
  try {
    parallel_for(4, 4, [&](std::size_t o) {
      parallel_for(100, 4, [&](std::size_t i) {
        if (o == 2 && i == 37) throw std::runtime_error("item 37");
      });
    });
  } catch (const std::runtime_error& e) {
    caught = std::string(e.what()) == "item 37";
  }
  CHECK(caught);
  std::atomic<long> sum{0};
  parallel_for(1000, 4, [&](std::size_t i) { sum += static_cast<long>(i); });
  CHECK_EQ(sum.load(), 999L * 1000 / 2);
}

TEST_CASE("parallel: cancelling stops nested loops") {
  Token token;
  std::atomic<long> ran{0};
  bool cancelled = false;
  try {
    parallel_for(
        8, 4,
        [&](std::size_t) {
          parallel_for(
              1000, 4,
              [&](std::size_t) {
                if (ran.fetch_add(1) == 50) token.cancel();
                std::this_thread::sleep_for(std::chrono::microseconds(50));
              },
              &token);
        },
        &token);
  } catch (const Cancelled&) {
    cancelled = true;
  }
  CHECK(cancelled);
  CHECK(ran.load() < 8000);
}

TEST_CASE("parallel: many callers at once, each with nested loops") {
  constexpr int CALLERS = 8;
  std::vector<long> total(CALLERS, 0);
  std::vector<std::thread> callers;
  for (int c = 0; c < CALLERS; ++c)
    callers.emplace_back([&, c] {
      std::vector<long> part(16, 0);
      parallel_for(part.size(), 4, [&](std::size_t o) {
        std::vector<long> inner(500, 0);
        parallel_for(inner.size(), 4, [&](std::size_t i) { inner[i] = static_cast<long>(c + o + i); });
        part[o] = std::accumulate(inner.begin(), inner.end(), 0L);
      });
      total[static_cast<std::size_t>(c)] = std::accumulate(part.begin(), part.end(), 0L);
    });
  for (auto& t : callers) t.join();
  for (int c = 0; c < CALLERS; ++c) {
    long want = 0;
    for (long o = 0; o < 16; ++o)
      for (long i = 0; i < 500; ++i) want += c + o + i;
    CHECK_EQ(total[static_cast<std::size_t>(c)], want);
  }
}
