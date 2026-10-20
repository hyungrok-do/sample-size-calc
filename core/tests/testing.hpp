// Minimal self-contained test framework (no external dependencies).
#pragma once

#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace ssc::test {

struct Case {
  std::string name;
  std::function<void()> fn;
};

std::vector<Case>& registry();
void report_failure(const std::string& msg, const char* file, int line);
void count_check();

struct Registrar {
  Registrar(const char* name, std::function<void()> fn) { registry().push_back({name, std::move(fn)}); }
};

inline bool close(double got, double expected, double rtol, double atol) {
  if (std::isnan(expected)) return std::isnan(got);
  if (std::isinf(expected)) return got == expected;
  return std::fabs(got - expected) <= atol + rtol * std::fabs(expected);
}

}  // namespace ssc::test

#define SSC_CAT2(a, b) a##b
#define SSC_CAT(a, b) SSC_CAT2(a, b)

#define TEST_CASE(name)                                                            \
  static void SSC_CAT(test_fn_, __LINE__)();                                     \
  static ::ssc::test::Registrar SSC_CAT(test_reg_, __LINE__)(name, SSC_CAT(test_fn_, __LINE__)); \
  static void SSC_CAT(test_fn_, __LINE__)()

#define CHECK(cond)                                                                \
  do {                                                                             \
    ::ssc::test::count_check();                                                    \
    if (!(cond)) ::ssc::test::report_failure("CHECK(" #cond ")", __FILE__, __LINE__); \
  } while (0)

#define CHECK_CLOSE(got, expected, rtol, atol)                                     \
  do {                                                                             \
    ::ssc::test::count_check();                                                    \
    const double ssc_g_ = (got), ssc_e_ = (expected);                              \
    if (!::ssc::test::close(ssc_g_, ssc_e_, (rtol), (atol))) {                     \
      std::ostringstream ssc_os_;                                                  \
      ssc_os_.precision(17);                                                       \
      ssc_os_ << #got << " = " << ssc_g_ << ", expected " << ssc_e_                \
              << " (rel err " << std::fabs(ssc_g_ - ssc_e_) / std::fabs(ssc_e_) << ")"; \
      ::ssc::test::report_failure(ssc_os_.str(), __FILE__, __LINE__);              \
    }                                                                              \
  } while (0)

#define CHECK_EQ(got, expected)                                                    \
  do {                                                                             \
    ::ssc::test::count_check();                                                    \
    const auto ssc_g_ = (got);                                                     \
    const auto ssc_e_ = (expected);                                                \
    if (!(ssc_g_ == ssc_e_)) {                                                     \
      std::ostringstream ssc_os_;                                                  \
      ssc_os_.precision(17);                                                       \
      ssc_os_ << #got << " = " << ssc_g_ << ", expected " << ssc_e_;               \
      ::ssc::test::report_failure(ssc_os_.str(), __FILE__, __LINE__);              \
    }                                                                              \
  } while (0)
