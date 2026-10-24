// P1 from the effect that the tests of two proportions are given (the test of two proportions and its
// non-inferiority and equivalence tests): a difference, a ratio or an odds ratio with P2, or P1 itself.
#pragma once

#include <cmath>
#include <string>

#include "procedures/common.hpp"

namespace ssc::procs {

/// P1 from the actual effect ("es": "diff" with "d1", "ratio" with "r1", "or" with "or1", or "p1").  An odds ratio
/// of 1 gives P2 itself (o/(1 + o) gave 0.5999999999999999 for 0.6, which the answer showed beside P2 = 0.6).  One
/// that rounds to 1 is refused (the power would need 1 − P1 in the statistic at the expected proportions as well),
/// named with what 1 − P1 is: P2 = 0.5 with D1 = 0.49999999999999994, or an odds ratio of 1e16, puts it a rounding
/// error below 1, where it was said to be 1.
inline double p1_from(const Scenario& s, double p2) {
  double p1, q1 = NAN;  // q1: 1 − P1, as the difference, the ratio or the odds give it
  if (s.is("es", "diff")) {
    const double d1 = s("d1");
    p1 = p2 + d1;
    const double z = p1 - p2;
    q1 = (1 - p1) - ((p2 - (p1 - z)) + (d1 - z));
  } else if (s.is("es", "ratio")) {
    p1 = p2 * s("r1");
    q1 = (1 - p1) - std::fma(s("r1"), p2, -p1);
  } else if (s.is("es", "or")) {
    const double o = s("or1") * p2 / (1 - p2);
    p1 = s("or1") == 1 ? p2 : (std::isinf(o) ? 1 : o / (1 + o));  // odds past the largest double: P1 is 1
    q1 = (1 - p2) / (1 - p2 + s("or1") * p2);
  } else
    p1 = s("p1");
  if (p1 == 1 && q1 > 0)
    throw CalcError("the implied P1 = 1 − " + blocks::fnum(q1) + " rounds to 1; it must lie strictly between 0 and "
                    "1 by more than that");
  if (!(p1 > 0 && p1 < 1))
    throw CalcError("the implied P1 = " +
                    (p1 >= 1 && q1 < 0 && q1 > -1e-4 ? "1 + " + blocks::fnum(-q1) : blocks::fnum(p1)) +
                    " is not strictly between 0 and 1");
  return p1;
}

}  // namespace ssc::procs
