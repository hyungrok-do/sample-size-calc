#include "procedures/groupseq.hpp"

#include <algorithm>
#include <cmath>
#include <optional>

#include "engine/procedure.hpp"
#include "stats/numerics.hpp"
#include "stats/stats.hpp"

namespace ssc::gs {

using stats::pnorm;
using stats::qnorm;

SpendFn spending(const std::string& name, double level, double param) {
  if (name == "obf") {
    // level(t) = 2 − 2Φ(z_{level/2} / √t)
    const double z = qnorm(level / 2, false);
    return [=](double t) { return t <= 0 ? 0.0 : (t >= 1 ? level : 2 * pnorm(z / std::sqrt(t), false)); };
  }
  if (name == "pocock") {
    return [=](double t) { return t <= 0 ? 0.0 : (t >= 1 ? level : level * std::log1p((std::exp(1.0) - 1) * t)); };
  }
  if (name == "power") {
    const double rho = param;
    return [=](double t) { return t <= 0 ? 0.0 : (t >= 1 ? level : level * std::pow(t, rho)); };
  }
  if (name == "hsd") {
    const double g = param;
    return [=](double t) {
      if (t <= 0) return 0.0;
      if (t >= 1) return level;
      if (std::fabs(g) < 1e-10) return level * t;
      return level * std::expm1(-g * t) / std::expm1(-g);
    };
  }
  throw CalcError("unknown spending function: " + name);
}

SpendFn alpha_spending(const std::string& name, double alpha, bool two_sided, double param) {
  if (!two_sided) return spending(name, alpha, param);
  SpendFn side = spending(name, alpha / 2, param);
  return [side](double t) { return 2 * side(t); };
}

namespace {

constexpr double INV_SQRT_2PI = 0.398942280401432677939946059934;
constexpr double TINY = 1e-300;

// Jennison–Turnbull grid on (lo, hi) for a normal variable with mean mu, with
// Simpson weights.  Points are dense near the mean and sparse in the tails.
void make_grid(double mu, double lo, double hi, int r, std::vector<double>& z, std::vector<double>& w) {
  std::vector<double> x;
  x.reserve(static_cast<std::size_t>(6 * r));
  for (int i = 1; i <= 6 * r - 1; ++i) {
    double v;
    if (i < r)
      v = mu - 3 - 4 * std::log(static_cast<double>(r) / i);
    else if (i <= 5 * r)
      v = mu - 3 + 3.0 * (i - r) / (2.0 * r);
    else
      v = mu + 3 + 4 * std::log(static_cast<double>(r) / (6 * r - i));
    x.push_back(v);
  }
  const double L = std::max(lo, x.front()), H = std::min(hi, x.back());
  if (!(H > L)) {
    z.assign(1, std::isfinite(L + H) ? 0.5 * (L + H) : mu);
    w.assign(1, 0.0);
    return;
  }
  std::vector<double> y;
  y.reserve(x.size() + 2);
  y.push_back(L);
  for (double v : x)
    if (v > L && v < H) y.push_back(v);
  y.push_back(H);
  const std::size_t n = y.size();
  z.assign(2 * n - 1, 0.0);
  w.assign(2 * n - 1, 0.0);
  for (std::size_t i = 0; i < n; ++i) z[2 * i] = y[i];
  for (std::size_t i = 0; i + 1 < n; ++i) {
    z[2 * i + 1] = 0.5 * (y[i] + y[i + 1]);
    const double d = (y[i + 1] - y[i]) / 6;
    w[2 * i] += d;
    w[2 * i + 1] += 4 * d;
    w[2 * i + 2] += d;
  }
}

// Sub-density of the z-statistic on the continuation region, carried from one
// analysis to the next.  At analysis k (starting at 0), upper(x) is
// P(Z_k ≥ x and no stop before k) and lower(x) is P(Z_k ≤ x and no stop
// before k).
class Track {
 public:
  Track(const std::vector<double>& t, double eta, int r) : t_(t), eta_(eta), r_(r) { set_stage(0); }

  std::size_t k() const { return k_; }
  double mean() const { return eta_ * stk_; }  // E(Z_k)
  double mass() const { return mass_; }        // P(no stop before k)

  double upper(double x) const {
    if (k_ == 0) return pnorm(x - mean(), false);
    if (x == -INFINITY) return mass_;
    if (x == INFINITY) return 0.0;
    const double X = (x * stk_ - eta_ * dt_) / sdt_;
    double s = 0;
    for (std::size_t i = 0; i < B_.size(); ++i) s += h_[i] * pnorm(X - B_[i], false);
    return s;
  }

  double lower(double x) const {
    if (k_ == 0) return pnorm(x - mean());
    if (x == -INFINITY) return 0.0;
    if (x == INFINITY) return mass_;
    const double X = (x * stk_ - eta_ * dt_) / sdt_;
    double s = 0;
    for (std::size_t i = 0; i < B_.size(); ++i) s += h_[i] * pnorm(X - B_[i]);
    return s;
  }

  // Keep the paths with lo < Z_k < hi and move on to analysis k + 1.
  void advance(double lo, double hi) {
    std::vector<double> z, w;
    make_grid(mean(), lo, hi, r_, z, w);
    std::vector<double> h(z.size(), 0.0);
    if (k_ == 0) {
      const double mu = mean();
      for (std::size_t j = 0; j < z.size(); ++j) {
        const double u = z[j] - mu;
        h[j] = w[j] * INV_SQRT_2PI * std::exp(-0.5 * u * u);
      }
    } else {
      // f_k(z) = ∫ f_{k−1}(u) φ((z√t_k − u√t_{k−1} − ηΔt)/√Δt) √t_k/√Δt du
      const double c = INV_SQRT_2PI * stk_ / sdt_, shift = eta_ * dt_;
      const double* B = B_.data();
      const double* hp = h_.data();
      const std::size_t m = B_.size();
      for (std::size_t j = 0; j < z.size(); ++j) {
        if (w[j] == 0) continue;
        const double A = (z[j] * stk_ - shift) / sdt_;
        double s = 0;
        for (std::size_t i = 0; i < m; ++i) {
          const double u = A - B[i];
          s += hp[i] * std::exp(-0.5 * u * u);
        }
        h[j] = w[j] * c * s;
      }
    }
    const double stk_old = stk_;
    set_stage(k_ + 1);
    B_.resize(z.size());
    for (std::size_t j = 0; j < z.size(); ++j) B_[j] = z[j] * stk_old / sdt_;
    h_ = std::move(h);
    mass_ = 0;
    for (double v : h_) mass_ += v;
  }

 private:
  void set_stage(std::size_t k) {
    k_ = k;
    stk_ = std::sqrt(t_[k]);
    dt_ = t_[k] - (k == 0 ? 0.0 : t_[k - 1]);
    sdt_ = std::sqrt(dt_);
  }

  const std::vector<double>& t_;
  double eta_;
  int r_;
  std::size_t k_ = 0;
  double stk_ = 0, dt_ = 0, sdt_ = 0, mass_ = 1;
  std::vector<double> B_, h_;  // previous-stage grid (scaled by √t_{k−1}/√Δt) and weighted sub-density
};

double log_ratio(double p, double target) { return std::log(std::max(p, 1e-320)) - std::log(target); }

// Efficacy boundary at the track's current analysis that spends `inc`.  Solved
// on the log scale so that very small increments are resolved.
double efficacy_bound(const Track& tr, double inc, bool two_sided) {
  if (!(inc > TINY)) return INFINITY;
  const double mu = tr.mean();
  if (tr.k() == 0 && mu == 0) return qnorm(two_sided ? inc / 2 : inc, false);
  auto g = [&](double x) { return log_ratio(two_sided ? tr.upper(x) + tr.lower(-x) : tr.upper(x), inc); };
  const double lo = two_sided ? 0.0 : mu - 12;
  const double glo = g(lo);
  if (glo <= 0) return lo;  // the increment exceeds what is left: every remaining path is rejected
  // the sub-density is dominated by the N(mu, 1) density, so g(hi) ≤ 0 up to rounding
  double hi = std::max(lo + 1e-3, std::fabs(mu) + qnorm(two_sided ? inc / 2 : inc, false));
  double ghi = g(hi);
  for (int i = 0; ghi > 0 && i < 100; ++i) ghi = g(hi += 0.25);
  return numerics::brent(g, lo, hi, glo, ghi, 1e-12);
}

// Futility boundary (at most b) at the track's current analysis that spends `inc`.
double futility_bound(const Track& tr, double inc, double b) {
  if (!(inc > TINY)) return -INFINITY;
  const double mu = tr.mean();
  auto g = [&](double x) { return log_ratio(tr.lower(x), inc); };
  const double hi = std::isfinite(b) ? b : mu + 12;
  const double ghi = g(hi);
  if (ghi <= 0) return hi;  // futility meets the efficacy boundary
  double lo = std::min(hi - 1e-3, mu - qnorm(inc, false));
  double glo = g(lo);
  for (int i = 0; glo > 0 && i < 100; ++i) glo = g(lo -= 0.25);
  return numerics::brent(g, lo, hi, glo, ghi, 1e-12);
}

void check_t(const std::vector<double>& t) {
  if (t.empty()) throw CalcError("at least one analysis is required");
  double prev = 0;
  for (double v : t) {
    if (!(v > prev)) throw CalcError("information fractions must be strictly increasing and positive");
    prev = v;
  }
  if (std::fabs(t.back() - 1) > 1e-9) throw CalcError("the last information fraction must be 1");
}

// Error to spend at each analysis: f(t_k) − f(t_{k−1}), with the total spent at the end.
std::vector<double> increments(const SpendFn& f, const std::vector<double>& t, double total) {
  std::vector<double> inc(t.size());
  double prev = 0;
  for (std::size_t k = 0; k < t.size(); ++k) {
    const double cum = k + 1 == t.size() ? total : std::min(total, f(t[k]));
    inc[k] = std::max(0.0, cum - prev);
    prev = std::max(prev, cum);
  }
  return inc;
}

// Alpha-spending efficacy boundaries under H0.  With binding futility
// boundaries `az` (analyses 1..K−1), paths below them leave the H0 track.
void spend_alpha(Design& d, const std::vector<double>& inc, const std::vector<double>* az, int r) {
  const std::size_t K = d.t.size();
  Track tr(d.t, 0.0, r);
  double cum = 0;
  for (std::size_t k = 0; k < K; ++k) {
    const double bk = efficacy_bound(tr, inc[k], d.two_sided);
    d.b[k] = bk;
    cum += d.two_sided ? tr.upper(bk) + tr.lower(-bk) : tr.upper(bk);
    d.alpha_spent[k] = cum;
    if (k + 1 == K) break;
    const double lo = d.two_sided ? -bk : (az ? std::min((*az)[k], bk) : -INFINITY);
    tr.advance(lo, bk);
  }
}

void cumulative_alpha(Design& d, int r) {
  const ExitProbs e = exit_probs(d, 0.0, 1.0, r);
  double c = 0;
  for (std::size_t k = 0; k < d.t.size(); ++k) {
    c += e.up[k] + e.lo[k];
    d.alpha_spent[k] = c;
  }
}

void set_fixed_futility(Design& d, const std::vector<double>& az) {
  const std::size_t K = d.t.size();
  d.a.assign(K, 0.0);
  for (std::size_t k = 0; k + 1 < K; ++k) d.a[k] = std::min(az[k], d.b[k]);
  d.a[K - 1] = d.b[K - 1];
}

// Wang–Tsiatis boundaries b_k = C t_k^(Δ − ½) with total level alpha
// (binding futility boundaries `az` included in the type I error).
void wang_tsiatis(Design& d, double alpha, double delta, const std::vector<double>* az,
                  const std::function<void()>& check, int r) {
  const std::size_t K = d.t.size();
  auto excess = [&](double C) {
    if (check) check();
    for (std::size_t k = 0; k < K; ++k) d.b[k] = C * std::pow(d.t[k], delta - 0.5);
    if (az) set_fixed_futility(d, *az);
    return exit_probs(d, 0.0, 1.0, r).reject - alpha;
  };
  const double z = qnorm(d.two_sided ? alpha / 2 : alpha, false);
  double lo = 0.5 * z, hi = 1.5 * z + 1;
  double flo = excess(lo), fhi = excess(hi);
  for (int i = 0; flo <= 0 && i < 60; ++i) flo = excess(lo *= 0.5);
  for (int i = 0; fhi >= 0 && i < 60; ++i) fhi = excess(hi *= 1.5);
  const double C = numerics::brent(excess, lo, hi, flo, fhi, 1e-12);
  excess(C);
  cumulative_alpha(d, r);
}

// One construction pass of beta-spending futility boundaries under drift eta:
// fills d.a (and d.b for binding designs, whose efficacy boundaries spend alpha
// among the paths that pass the futility boundaries).  Returns the power.
double beta_pass(Design& d, double eta, const std::vector<double>& ainc, const std::vector<double>& binc, int r) {
  const std::size_t K = d.t.size();
  Track t1(d.t, eta, r);
  std::optional<Track> t0;
  if (d.binding) t0.emplace(d.t, 0.0, r);
  double power = 0, cum_a = 0, cum_b = 0;
  for (std::size_t k = 0; k < K; ++k) {
    if (t0) {
      d.b[k] = efficacy_bound(*t0, ainc[k], false);
      cum_a += t0->upper(d.b[k]);
      d.alpha_spent[k] = cum_a;
    }
    const double bk = d.b[k];
    power += t1.upper(bk);
    const double ak = k + 1 < K ? futility_bound(t1, binc[k], bk) : bk;
    d.a[k] = ak;
    cum_b += t1.lower(ak);
    d.beta_spent[k] = cum_b;
    if (k + 1 == K) break;
    t1.advance(ak, bk);
    if (t0) t0->advance(ak, bk);
  }
  return power;
}

// Beta-spending futility: the drift at which the boundaries spend exactly beta.
void spend_beta(Design& d, double alpha, const std::vector<double>& ainc, const Futility& fut,
                const std::function<void()>& check, int r) {
  const auto binc = increments(fut.spend, d.t, fut.beta);
  const double target = 1 - fut.beta;
  auto f = [&](double eta) {
    if (check) check();
    return beta_pass(d, eta, ainc, binc, r) - target;
  };
  const double eta_f = qnorm(alpha, false) + qnorm(fut.beta, false);  // fixed-sample drift
  if (!(eta_f > 0)) throw CalcError("the power must exceed alpha");
  double lo = 0.95 * eta_f, hi = 1.3 * eta_f;
  double flo = f(lo), fhi = f(hi);
  for (int i = 0; flo > 0 && i < 60; ++i) {
    hi = lo;
    fhi = flo;
    flo = f(lo *= 0.8);
  }
  for (int i = 0; fhi < 0 && i < 60; ++i) {
    lo = hi;
    flo = fhi;
    fhi = f(hi *= 1.3);
  }
  if (!(flo <= 0 && fhi >= 0)) throw CalcError("the futility boundaries could not be determined");
  const double eta = numerics::brent(f, lo, hi, flo, fhi, 1e-11);
  beta_pass(d, eta, ainc, binc, r);
  d.eta_design = eta;
}

// Increments of alpha of at most TINY are not spent (their boundaries are infinite): a design that could spend
// none of its alpha would never reject.
void check_spent(const std::vector<double>& inc) {
  for (double v : inc)
    if (v > TINY) return;
  throw CalcError("alpha is too small to compute with");
}

}  // namespace

Design make_design(const std::vector<double>& t, double alpha, bool two_sided, const Efficacy& eff,
                   const Futility& fut, const std::function<void()>& check, int r) {
  check_t(t);
  if (!(alpha > 0 && alpha < (two_sided ? 1.0 : 0.5)))
    throw CalcError(two_sided ? "alpha must lie between 0 and 1" : "a one-sided alpha must be below 0.5");
  const std::size_t K = t.size();
  const bool futility = fut.kind != Futility::Kind::None && K > 1;
  if (futility && two_sided) throw CalcError("futility boundaries are available for one-sided designs only");
  if (futility && fut.kind == Futility::Kind::BetaSpending && eff.kind != Efficacy::Kind::Spending)
    throw CalcError("beta-spending futility boundaries need an alpha-spending efficacy boundary");
  if (futility && fut.kind == Futility::Kind::Fixed && fut.z.size() != K - 1)
    throw CalcError("enter one futility boundary per interim analysis (" + std::to_string(K - 1) + ")");

  Design d;
  d.t = t;
  d.two_sided = two_sided;
  d.binding = futility && fut.binding;
  d.b.assign(K, INFINITY);
  d.alpha_spent.assign(K, 0.0);
  const std::vector<double>* az = futility && fut.kind == Futility::Kind::Fixed && d.binding ? &fut.z : nullptr;

  if (eff.kind == Efficacy::Kind::Spending) {
    const auto ainc = increments(eff.spend, t, alpha);
    check_spent(ainc);
    if (futility && fut.kind == Futility::Kind::BetaSpending) {
      d.a.assign(K, -INFINITY);
      d.beta_spent.assign(K, 0.0);
      if (!d.binding) spend_alpha(d, ainc, nullptr, r);
      spend_beta(d, alpha, ainc, fut, check, r);
      return d;
    }
    spend_alpha(d, ainc, az, r);
  } else {
    check_spent({alpha});
    wang_tsiatis(d, alpha, eff.delta, az, check, r);
  }
  if (futility) set_fixed_futility(d, fut.z);
  return d;
}

ExitProbs exit_probs(const Design& d, double eta, double scale, int r) {
  const std::size_t K = d.t.size();
  ExitProbs e;
  e.up.assign(K, 0.0);
  e.lo.assign(K, 0.0);
  e.fut.assign(K, 0.0);
  Track tr(d.t, eta, r);
  const bool fut = d.has_futility();
  double stopped = 0;  // P(stop before the current analysis)
  for (std::size_t k = 0; k < K; ++k) {
    const double b = scale * d.b[k];
    e.up[k] = tr.upper(b);
    if (d.two_sided) e.lo[k] = tr.lower(-b);
    if (fut) e.fut[k] = tr.lower(scale * d.a[k]);
    e.reject += e.up[k] + e.lo[k];
    if (k + 1 == K) {
      e.info += std::max(0.0, 1 - stopped) * d.t[k];
      break;
    }
    const double stop = e.up[k] + e.lo[k] + e.fut[k];
    stopped += stop;
    e.info += stop * d.t[k];
    tr.advance(d.two_sided ? -b : (fut ? scale * d.a[k] : -INFINITY), b);
  }
  e.reject = std::min(1.0, e.reject);  // the quadrature error can push a sum near 1 just above it
  return e;
}

double drift_for_power(const Design& d, double power, double scale, int r) {
  auto f = [&](double eta) { return exit_probs(d, eta, scale, r).reject - power; };
  const double side = d.two_sided ? d.alpha_spent.back() / 2 : d.alpha_spent.back();
  const double zp = qnorm(power);
  // the fixed-sample test is most powerful at a given drift, so its drift is a lower bound
  double lo = std::max(1e-6, scale * qnorm(side, false) + zp - 0.25);
  double hi = std::max(lo + 0.5, scale * d.b.back() + zp + 0.25);
  double flo = f(lo), fhi = f(hi);
  for (int i = 0; flo > 0 && lo > 1e-9 && i < 60; ++i) flo = f(lo *= 0.5);
  if (flo > 0) return lo;
  for (int i = 0; fhi < 0; ++i) {
    if (i > 60 || hi > 1e3) throw CalcError("the target power cannot be reached");
    lo = hi;
    flo = fhi;
    fhi = f(hi = 1.25 * hi + 0.5);
  }
  return numerics::brent(f, lo, hi, flo, fhi, 1e-12);
}

}  // namespace ssc::gs
