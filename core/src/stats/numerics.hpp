// Header-only numerical utilities: root finding and quadrature.
// Templates so that lambdas inline into the hot loops.
#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "stats.hpp"

namespace ssc::numerics {

/// Brent's method (zeroin): root of f in [a, b] given f(a), f(b) of opposite sign.
template <class F>
double brent(F&& f, double a, double b, double fa, double fb, double tol = 1e-12, int max_iter = 300) {
  if (fa == 0) return a;
  if (fb == 0) return b;
  if ((fa > 0) == (fb > 0)) throw std::domain_error("root is not bracketed");
  double c = a, fc = fa, d = b - a, e = d;
  for (int iter = 0; iter < max_iter; ++iter) {
    if ((fb > 0) == (fc > 0)) {
      c = a;
      fc = fa;
      d = e = b - a;
    }
    if (std::fabs(fc) < std::fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    const double tol1 = 2 * stats::DBL_EPSILON_ * std::fabs(b) + 0.5 * tol;
    const double xm = 0.5 * (c - b);
    if (std::fabs(xm) <= tol1 || fb == 0) return b;
    if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
      const double s = fb / fa;
      double p, q;
      if (a == c) {
        p = 2 * xm * s;
        q = 1 - s;
      } else {
        const double qq = fa / fc, r = fb / fc;
        p = s * (2 * xm * qq * (qq - r) - (b - a) * (r - 1));
        q = (qq - 1) * (r - 1) * (s - 1);
      }
      if (p > 0) q = -q;
      p = std::fabs(p);
      if (2 * p < std::min(3 * xm * q - std::fabs(tol1 * q), std::fabs(e * q))) {
        e = d;
        d = p / q;
      } else {
        d = xm;
        e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    b += (std::fabs(d) > tol1) ? d : (xm > 0 ? tol1 : -tol1);
    fb = f(b);
  }
  return b;
}

template <class F>
double brent(F&& f, double a, double b, double tol = 1e-12, int max_iter = 300) {
  return brent(f, a, b, f(a), f(b), tol, max_iter);
}

/// 15-point Gauss–Kronrod rule on [a, b]; err receives |K15 − G7|.
template <class F>
double gk15(F& f, double a, double b, double& err) {
  static constexpr double xgk[8] = {0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
                                    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
                                    0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
                                    0.207784955007898467600689403773245, 0.0};
  static constexpr double wgk[8] = {0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
                                    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
                                    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
                                    0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
  static constexpr double wg[4] = {0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
                                   0.381830050505118944950369775488975, 0.417959183673469387755102040816327};
  const double c = 0.5 * (a + b), h = 0.5 * (b - a);
  const double fc = f(c);
  double resk = fc * wgk[7], resg = fc * wg[3];
  for (int j = 0; j < 7; ++j) {
    const double dx = h * xgk[j];
    const double s = f(c - dx) + f(c + dx);
    resk += wgk[j] * s;
    if (j % 2 == 1) resg += wg[j / 2] * s;
  }
  err = std::fabs((resk - resg) * h);
  return resk * h;
}

namespace detail {
template <class F>
double adapt(F& f, double a, double b, double eps_abs, double eps_rel, int depth, double whole) {
  double err;
  const double r = gk15(f, a, b, err);
  // A NaN at a node of a removable point clears in a split or two; a NaN integrand is not split to max_depth.
  if (std::isnan(r)) depth = std::min(depth, 4);
  if (depth <= 0 || err <= std::max(eps_abs, eps_rel * std::fabs(whole)) || !(b - a > 1e-300)) return r;
  const double m = 0.5 * (a + b);
  return adapt(f, a, m, eps_abs * 0.5, eps_rel, depth - 1, whole) +
         adapt(f, m, b, eps_abs * 0.5, eps_rel, depth - 1, whole);
}
}  // namespace detail

/// Adaptive Gauss–Kronrod integration of f over [a, b].
template <class F>
double integrate(F&& f, double a, double b, double eps_abs = 1e-13, double eps_rel = 1e-11, int max_depth = 30) {
  if (a == b) return 0;
  double err;
  const double whole = gk15(f, a, b, err);
  if (err <= std::max(eps_abs, eps_rel * std::fabs(whole))) return whole;
  const double m = 0.5 * (a + b);
  return detail::adapt(f, a, m, eps_abs * 0.5, eps_rel, max_depth, whole) +
         detail::adapt(f, m, b, eps_abs * 0.5, eps_rel, max_depth, whole);
}

/// Fixed-order Gauss–Legendre integration of f over [a, b].
template <class F>
double gauss_legendre_integrate(F&& f, double a, double b, int n) {
  const stats::GLRule& gl = stats::gauss_legendre(n);
  const double c = 0.5 * (a + b), h = 0.5 * (b - a);
  double s = 0;
  for (int i = 0; i < n; ++i) s += gl.w[i] * f(c + h * gl.x[i]);
  return s * h;
}

}  // namespace ssc::numerics
