// Numerical core: special functions and probability distributions.
//
// Every function here is pure and thread-safe (no hidden global state), so the
// engine can evaluate many scenarios concurrently.  Accuracy targets are close
// to double precision for the central distributions and ~1e-13 absolute for the
// noncentral ones; see tests/ for the reference comparisons.
#pragma once

#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace ssc::stats {

constexpr double LN_SQRT_2PI = 0.918938533204672741780329736406;  // ln √(2π)
constexpr double LN_2PI = 1.837877066409345483560659472811;       // ln 2π
constexpr double SQRT_2PI = 2.506628274631000502415765284811;     // √(2π)
constexpr double INV_SQRT_2PI = 0.398942280401432677939946059934; // 1/√(2π)
constexpr double PI = 3.141592653589793238462643383280;
constexpr double DBL_EPSILON_ = 2.220446049250313e-16;
constexpr double DBL_MIN_ = 2.2250738585072014e-308;
constexpr double INF = std::numeric_limits<double>::infinity();
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

// ---------------------------------------------------------------- special
/// Stirling remainder δ(x) = lnΓ(x) − [(x − ½)ln x − x + ln√(2π)], x > 0.
double stirlerr(double x);
/// lnΓ(x) (thread-safe replacement for std::lgamma).
double lgammafn(double x);
/// lnΓ(b + a) − lnΓ(b) without cancellation for large b.
double lgamma_diff(double b, double a);
/// ln B(a, b).
double lbeta(double a, double b);
/// Deviance term x ln(x/m) + m − x, accurate when x ≈ m.
double bd0(double x, double m);
/// λ^x e^{−λ}/Γ(x + 1) for real x ≥ 0 (saddle-point form).
double dpois_raw(double x, double lambda);
/// C(n, x) p^x q^{n−x} for real 0 ≤ x ≤ n (saddle-point form, q = 1 − p).
double dbinom_raw(double x, double n, double p, double q);
/// x^a y^b / B(a, b) with y = 1 − x, to full relative precision.
double beta_prefactor(double a, double b, double x, double y);
/// ln C(n, k) for real arguments.
double lchoose(double n, double k);

inline double clamp01(double p) { return p < 0 ? 0 : (p > 1 ? 1 : p); }

// ---------------------------------------------------------------- normal
/// {Φ(x), 1 − Φ(x)}, each to near full relative precision (Cody 1969).
std::pair<double, double> pnorm_both(double x);
double pnorm(double x, bool lower = true);
double dnorm(double x);
/// Normal quantile (Wichura AS 241). lower=false: z with 1 − Φ(z) = p.
double qnorm(double p, bool lower = true);

// ---------------------------------------------------------------- gamma family
/// Regularized incomplete gamma P(a, x) (lower) or Q(a, x).
double pgamma_raw(double x, double a, bool lower = true);
/// P(a, a + d) (or Q): the offset form keeps x − a exact when a is huge.
double pgamma_offset(double a, double d, bool lower = true);
/// Gamma(shape a, scale 1) density.
double dgamma_std(double x, double a);
/// Gamma(shape a, scale 1) quantile.
double qgamma(double p, double a, bool lower = true);
double pchisq(double x, double df, bool lower = true);
double qchisq(double p, double df, bool lower = true);
double dchisq(double x, double df);
double dpois(double k, double lambda);
double ppois(double k, double lambda, bool lower = true);
double qpois(double p, double lambda, bool lower = true);

// ---------------------------------------------------------------- beta family
/// Regularized incomplete beta I_x(a, b) (lower) or 1 − I_x(a, b); y = 1 − x.
double pbeta_raw(double x, double a, double b, bool lower, double y);
double pbeta(double x, double a, double b, bool lower = true);
double dbeta(double x, double a, double b, double y);
/// Beta quantile as the pair {x, 1 − x}, both to full relative precision.
std::pair<double, double> qbeta_pair(double p, double a, double b, bool lower = true);
double qbeta(double p, double a, double b, bool lower = true);
/// Student t (df may be non-integer or +∞).
double pt(double t, double df, bool lower = true);
double dt(double t, double df);
double qt(double p, double df, bool lower = true);
/// F distribution.
double pf(double f, double df1, double df2, bool lower = true);
double qf(double p, double df1, double df2, bool lower = true);
/// Binomial.
double dbinom(double k, double n, double p);
double pbinom(double k, double n, double p, bool lower = true);
/// Smallest k with P(X ≤ k) ≥ p.
double qbinom(double p, double n, double prob);

// ---------------------------------------------------------------- noncentral
/// Noncentral beta CDF Σ Pois(j; λ/2) I_x(a + j, b).
double pnbeta(double x, double a, double b, double ncp, bool lower, double y);
double pnf(double f, double df1, double df2, double ncp, bool lower = true);
double pnchisq(double x, double df, double ncp, bool lower = true);
/// Noncentral t CDF.
double pnt(double t, double df, double ncp, bool lower = true);

// ---------------------------------------------------------------- hypergeometric
/// P(X = x) for X ~ Hypergeometric(m white, n black, k drawn).
double dhyper(double x, double m, double n, double k);
/// P(X ≤ x) (lower) or P(X > x).
double phyper(double x, double m, double n, double k, bool lower = true);

// ---------------------------------------------------------------- bivariate normal
/// P(X > h, Y > k) for standard bivariate normal with correlation r (Genz 2004).
double bvn_upper(double h, double k, double r);
/// P(X ≤ h, Y ≤ k).
/// Exact distribution of the Pearson correlation of n bivariate normal pairs with correlation rho.
double dcorr(double r, double n, double rho);
double pcorr(double x, double n, double rho, bool lower = true);

double pbvnorm(double h, double k, double r);

// ---------------------------------------------------------------- quadrature
struct GLRule {
  std::vector<double> x;  // nodes on [−1, 1]
  std::vector<double> w;  // weights
};
/// Gauss–Legendre rule with n nodes (cached; thread-safe).
const GLRule& gauss_legendre(int n);

}  // namespace ssc::stats
