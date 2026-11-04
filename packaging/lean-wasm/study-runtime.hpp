// Shared numerical backend for browser WASM and downloadable C++17 studies.
// Exact rational model arithmetic is checked separately by Lean. These normal
// tails and floating-point comparisons are numerical, not kernel certificates.
#pragma once
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace ssc_lean {
struct Group { std::string label; int64_t units, participants; };
struct Plan {
  double variance, effect, alpha, target;
  int64_t minimum, maximum, fixed;
  std::string mode;
  std::vector<Group> groups;
};
inline std::string quote(const std::string &text) {
  std::ostringstream out; out << '"';
  for (unsigned char c : text) {
    if (c == '"' || c == '\\') out << '\\' << c;
    else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
    else out << c;
  }
  out << '"'; return out.str();
}
inline double upper_tail(double z) { return 0.5 * std::erfc(z / std::sqrt(2.0)); }
inline double critical_value(double alpha) {
  double lo=0, hi=40;
  for(int i=0;i<120;++i) { double mid=(lo+hi)/2; if(upper_tail(mid)>alpha/2) lo=mid; else hi=mid; }
  return (lo+hi)/2;
}
inline double normal_quantile(double probability) {
  double lo=-40, hi=40;
  for(int i=0;i<120;++i) { double mid=(lo+hi)/2; if(upper_tail(-mid)<probability) lo=mid; else hi=mid; }
  return (lo+hi)/2;
}
inline double power_lambda(double lambda, double z) {
  const double p=upper_tail(z+lambda)+upper_tail(z-lambda);
  if(!std::isfinite(p)) throw std::runtime_error("Non-finite power; rescale planning assumptions");
  return p;
}
inline std::string evaluate(const Plan &p) {
  if(!std::isfinite(p.variance) || p.variance<=0 || !std::isfinite(p.effect))
    throw std::runtime_error("Planning variance and effect must be finite, with positive variance");
  if(!(p.alpha>=1e-12 && p.alpha<p.target && p.target<=1-1e-12))
    throw std::runtime_error("Numerical evaluation requires 1e-12 <= alpha < targetPower <= 1-1e-12");
  if(p.minimum<1 || p.minimum>p.maximum || p.maximum>10000000 || p.fixed<0 || p.fixed>10000000)
    throw std::runtime_error("Invalid integer search domain");
  if(p.groups.empty()) throw std::runtime_error("No allocation groups");
  int64_t multiplier=0;
  for(const auto &g:p.groups) {
    if(g.units<1 || g.units>1000000 || g.participants<1 || g.participants>1000000)
      throw std::runtime_error("Invalid allocation");
    multiplier+=g.units*g.participants;
  }
  if(multiplier>1000000) throw std::runtime_error("Allocation exceeds 1,000,000 participants per block");
  const double z=critical_value(p.alpha);
  const auto lambda=[&](int64_t n) {
    const double x=std::abs(p.effect)/std::sqrt(p.variance)*std::sqrt(double(n));
    if(!std::isfinite(x)) throw std::runtime_error("Non-finite noncentrality; rescale assumptions");
    return x;
  };
  const auto power=[&](int64_t n) { return power_lambda(lambda(n),z); };
  std::ostringstream out; out << std::setprecision(17);
  out << "{\"kernelCertified\":false,\"criticalValue\":" << z
      << ",\"varianceCoefficient\":" << p.variance << ",\"alpha\":" << p.alpha << ",\"targetPower\":" << p.target;
  const auto counts=[&](int64_t n) {
    out << ",\"totalParticipants\":" << n*multiplier << ",\"groups\":[";
    bool first=true;
    for(const auto &g:p.groups) {
      if(!first) out << ','; first=false;
      out << "{\"label\":" << quote(g.label) << ",\"count\":" << n*g.units
          << ",\"participants\":" << n*g.units*g.participants << '}';
    }
    out << ']';
  };
  if(p.mode=="sample_size") {
    if(p.effect!=0) {
      double approximate=p.variance*std::pow(z+normal_quantile(p.target),2)/(p.effect*p.effect);
      if(std::isfinite(approximate)) out << ",\"dominantTailApproximation\":" << approximate;
    }
    double previous=0;
    for(int64_t n=p.minimum;n<=p.maximum;++n) {
      double value=power(n);
      if(std::abs(value-p.target)<=1e-12) {
        out << ",\"status\":\"numerically_unresolved\",\"unresolvedN\":" << n << ",\"power\":" << value << '}'; return out.str();
      }
      if(value>p.target) {
        out << ",\"status\":\"found\",\"selectedN\":" << n << ",\"power\":" << value
            << ",\"standardError\":" << std::sqrt(p.variance)/std::sqrt(double(n)) << ",\"noncentrality\":" << lambda(n)
            << ",\"evaluated\":" << n-p.minimum+1;
        if(n>p.minimum) out << ",\"predecessorN\":" << n-1 << ",\"predecessorPower\":" << previous;
        counts(n); out << ",\"minimality\":\"all_smaller_admissible_counts_fail_in_binary64\"}"; return out.str();
      }
      previous=value; if(p.effect==0) break;
    }
    out << ",\"status\":\"not_found_in_domain\",\"powerAtLastEvaluated\":" << previous;
  } else if(p.mode=="power" && p.fixed>0) {
    out << ",\"status\":\"evaluated\",\"n\":" << p.fixed << ",\"power\":" << power(p.fixed)
        << ",\"standardError\":" << std::sqrt(p.variance)/std::sqrt(double(p.fixed)) << ",\"noncentrality\":" << lambda(p.fixed);
    counts(p.fixed);
  } else if(p.mode=="effect" && p.fixed>0) {
    double lo=0,hi=64;
    for(int i=0;i<120;++i) {
      double mid=(lo+hi)/2; if(mid==lo || mid==hi) break;
      if(power_lambda(mid,z)>=p.target) hi=mid; else lo=mid;
    }
    double scale=std::sqrt(p.variance)/std::sqrt(double(p.fixed));
    if(!std::isfinite(hi*scale)) throw std::runtime_error("Detectable effect overflowed");
    out << ",\"status\":\"approximated\",\"n\":" << p.fixed << ",\"detectableMagnitude\":" << hi*scale
        << ",\"effectLower\":" << lo*scale << ",\"effectUpper\":" << hi*scale
        << ",\"powerAtLower\":" << power_lambda(lo,z) << ",\"powerAtUpper\":" << power_lambda(hi,z);
    counts(p.fixed);
  } else throw std::runtime_error("Unsupported study query");
  out << '}'; return out.str();
}
}
