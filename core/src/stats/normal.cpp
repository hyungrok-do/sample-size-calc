// Standard normal distribution.
//   pnorm: W. J. Cody, "Rational Chebyshev approximations for the error
//          function", Math. Comp. 23 (1969) — near full precision in both tails.
//   qnorm: M. J. Wichura, Algorithm AS 241 (PPND16), Appl. Statist. 37 (1988).
#include "stats.hpp"

namespace ssc::stats {

namespace {
constexpr double A[5] = {2.2352520354606839287, 161.02823106855587881, 1067.6894854603709582,
                         18154.981253343561249, 0.065682337918207449113};
constexpr double B[4] = {47.20258190468824187, 976.09855173777669322, 10260.932208618978205,
                         45507.789335026729956};
constexpr double C[9] = {0.39894151208813466764, 8.8831497943883759412, 93.506656132177855979,
                         597.27027639480026226,  2494.5375852903726711, 6848.1904505362823326,
                         11602.651437647350124,  9842.7148383839780218, 1.0765576773720192317e-8};
constexpr double D[8] = {22.266688044328115691, 235.38790178262499861, 1519.377599407554805,
                         6485.558298266760755,  18615.571640885098091, 34900.952721145977266,
                         38912.003286093271411, 19685.429676859990727};
constexpr double P[6] = {0.21589853405795699,     0.1274011611602473639,  0.022235277870649807,
                         0.001421619193227893466, 2.9112874951168792e-5, 0.02307344176494017303};
constexpr double Q[5] = {1.28426009614491121, 0.468238212480865118, 0.0659881378689285515,
                         0.00378239633202758244, 7.29751555083966205e-5};
constexpr double SQRT32 = 5.656854249492380195206754896838;

/// exp(−y²/2) to full precision with one exponential (Cody, and R, use two, with y² split at a sixteenth):
/// y² = h + l exactly (Dekker's product of y's halves of 26 bits), and exp(−(h + l)/2) = exp(−h/2)(1 − l/2),
/// since |l| ≤ ulp(h)/2.  Against 40-digit values over |x| ≤ 38.4, Φ's smaller tail is within 5 ulp (mean
/// 0.67 ulp; with two exponentials 6 and 0.70), and Φ takes a quarter less time.
double exp_half_square(double y) {
  const double c = 134217729.0 * y;  // 2^27 + 1
  const double hi = c - (c - y), lo = y - hi;
  const double h = y * y, l = ((hi * hi - h) + 2 * hi * lo) + lo * lo;
  return std::exp(-0.5 * h) * (1 - 0.5 * l);
}
}  // namespace

std::pair<double, double> pnorm_both(double x) {
  if (std::isnan(x)) return {NaN, NaN};
  const double y = std::fabs(x);
  if (y <= 0.67448975) {
    double xnum = 0, xden = 0;
    if (y > DBL_EPSILON_ * 0.5) {
      const double xsq = x * x;
      xnum = A[4] * xsq;
      xden = xsq;
      for (int i = 0; i < 3; ++i) {
        xnum = (xnum + A[i]) * xsq;
        xden = (xden + B[i]) * xsq;
      }
    }
    const double temp = x * (xnum + A[3]) / (xden + B[3]);
    return {0.5 + temp, 0.5 - temp};
  }
  double r;
  if (y <= SQRT32) {
    double xnum = C[8] * y, xden = y;
    for (int i = 0; i < 7; ++i) {
      xnum = (xnum + C[i]) * y;
      xden = (xden + D[i]) * y;
    }
    const double temp = (xnum + C[7]) / (xden + D[7]);
    r = exp_half_square(y) * temp;
  } else if (y < 50) {
    const double xsq0 = 1.0 / (x * x);
    double xnum = P[5] * xsq0, xden = xsq0;
    for (int i = 0; i < 4; ++i) {
      xnum = (xnum + P[i]) * xsq0;
      xden = (xden + Q[i]) * xsq0;
    }
    double temp = xsq0 * (xnum + P[4]) / (xden + Q[4]);
    temp = (INV_SQRT_2PI - temp) / y;
    r = exp_half_square(y) * temp;
  } else {
    r = 0;
  }
  return x > 0 ? std::make_pair(1.0 - r, r) : std::make_pair(r, 1.0 - r);
}

double pnorm(double x, bool lower) {
  const auto pr = pnorm_both(x);
  return lower ? pr.first : pr.second;
}

double dnorm(double x) { return INV_SQRT_2PI * std::exp(-0.5 * x * x); }

double qnorm(double p, bool lower) {
  if (std::isnan(p) || p < 0 || p > 1) return NaN;
  if (p == 0) return lower ? -INF : INF;
  if (p == 1) return lower ? INF : -INF;
  const double q = lower ? p - 0.5 : 0.5 - p;
  if (std::fabs(q) <= 0.425) {
    const double r = 0.180625 - q * q;
    return q *
           (((((((2509.0809287301226727 * r + 33430.575583588128105) * r + 67265.770927008700853) * r +
                45921.953931549871457) * r + 13731.693765509461125) * r + 1971.5909503065514427) * r +
             133.14166789178437745) * r + 3.387132872796366608) /
           (((((((5226.495278852545925 * r + 28729.085735721942674) * r + 39307.89580009271061) * r +
                21213.794301586595867) * r + 5394.1960214247511077) * r + 687.1870074920579083) * r +
             42.313330701600911252) * r + 1.0);
  }
  // probability in the smaller tail
  const double s = q < 0 ? (lower ? p : 1 - p) : (lower ? 1 - p : p);
  double r = std::sqrt(-std::log(s));
  double val;
  if (r <= 5.0) {
    r -= 1.6;
    val = (((((((7.7454501427834140764e-4 * r + 0.0227238449892691845833) * r + 0.24178072517745061177) * r +
               1.27045825245236838258) * r + 3.64784832476320460504) * r + 5.7694972214606914055) * r +
            4.6303378461565452959) * r + 1.42343711074968357734) /
          (((((((1.05075007164441684324e-9 * r + 5.475938084995344946e-4) * r + 0.0151986665636164571966) * r +
               0.14810397642748007459) * r + 0.68976733498510000455) * r + 1.6763848301838038494) * r +
            2.05319162663775882187) * r + 1.0);
  } else {
    r -= 5.0;
    val = (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r + 0.0012426609473880784386) * r +
               0.026532189526576123093) * r + 0.29656057182850489123) * r + 1.7848265399172913358) * r +
            5.4637849111641143699) * r + 6.6579046435011037772) /
          (((((((2.04426310338993978564e-15 * r + 1.4215117583164458887e-7) * r + 1.8463183175100546818e-5) * r +
               7.868691311456132591e-4) * r + 0.0148753612908506148525) * r + 0.13692988092273580531) * r +
            0.59983220655588793769) * r + 1.0);
  }
  return q < 0 ? -val : val;
}

}  // namespace ssc::stats
