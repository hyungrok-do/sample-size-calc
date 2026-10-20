#include <chrono>
#include <cstring>

#include "testing.hpp"

namespace ssc::test {

namespace {
int g_checks = 0;
int g_failures = 0;
const char* g_current = "";
}  // namespace

std::vector<Case>& registry() {
  static std::vector<Case> cases;
  return cases;
}

void count_check() { ++g_checks; }

void report_failure(const std::string& msg, const char* file, int line) {
  ++g_failures;
  if (g_failures <= 200) std::cerr << "  FAIL [" << g_current << "] " << file << ":" << line << ": " << msg << "\n";
}

}  // namespace ssc::test

int main(int argc, char** argv) {
  using namespace ssc::test;
  const char* filter = argc > 1 ? argv[1] : nullptr;
  int ran = 0;
  const auto t0 = std::chrono::steady_clock::now();
  for (const auto& c : registry()) {
    if (filter && std::strstr(c.name.c_str(), filter) == nullptr) continue;
    g_current = c.name.c_str();
    const int before = g_failures;
    const auto s = std::chrono::steady_clock::now();
    try {
      c.fn();
    } catch (const std::exception& e) {
      report_failure(std::string("uncaught exception: ") + e.what(), "?", 0);
    }
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - s).count();
    std::cout << (g_failures == before ? "[ ok ] " : "[FAIL] ") << c.name << "  (" << ms << " ms)\n";
    ++ran;
  }
  const double total = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  std::cout << "\n" << ran << " test cases, " << g_checks << " checks, " << g_failures << " failures (" << total
            << " ms)\n";
  return g_failures == 0 ? 0 : 1;
}
