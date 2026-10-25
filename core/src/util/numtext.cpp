#include "numtext.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

namespace ssc::text {

namespace {
bool ieq(const std::string& s, std::size_t pos, const char* word) {
  std::size_t k = 0;
  for (; word[k]; ++k) {
    if (pos + k >= s.size()) return false;
    if (std::tolower(static_cast<unsigned char>(s[pos + k])) != word[k]) return false;
  }
  return true;
}
}  // namespace

bool parse_double(const std::string& s, std::size_t& pos, double& out) {
  std::size_t i = pos;
  bool neg = false;
  if (i < s.size() && (s[i] == '+' || s[i] == '-')) {
    neg = s[i] == '-';
    ++i;
  }
  if (ieq(s, i, "infinity") || ieq(s, i, "inf")) {
    out = neg ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
    pos = i + (ieq(s, i, "infinity") ? 8 : 3);
    return true;
  }
  if (ieq(s, i, "nan")) {
    out = std::numeric_limits<double>::quiet_NaN();
    pos = i + 3;
    return true;
  }
  // scan the extent of a decimal number: digits [. digits] [e[+-]digits]
  const std::size_t start = i;
  std::size_t digits = 0;
  while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) ++i, ++digits;
  if (i < s.size() && s[i] == '.') {
    ++i;
    while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) ++i, ++digits;
  }
  if (digits == 0) return false;
  if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
    std::size_t j = i + 1;
    if (j < s.size() && (s[j] == '+' || s[j] == '-')) ++j;
    if (j < s.size() && std::isdigit(static_cast<unsigned char>(s[j]))) {
      while (j < s.size() && std::isdigit(static_cast<unsigned char>(s[j]))) ++j;
      i = j;
    }
  }
  std::istringstream is(s.substr(start, i - start));
  is.imbue(std::locale::classic());
  double v = 0;
  is >> v;
  if (is.fail()) {
    // Out-of-range literals saturate like strtod.  Some C++ libraries (libc++)
    // also report subnormal results as out of range: those are parsed with the
    // exponent shifted into the normal range and scaled back (within one ulp).
    // Which way it is out of range is told by the place of the leading digit,
    // not by the sign of the exponent ("1000...0e-1" with 320 zeros is past
    // the largest double, "0.000...01" with 400 zeros below the smallest).
    const std::string tok = s.substr(start, i - start);
    const auto e = tok.find_first_of("eE");
    const std::string mantissa = tok.substr(0, e);
    long ex = 0;  // the exponent, up to a size no double comes near
    if (e != std::string::npos) {
      const bool down = tok[e + 1] == '-';
      for (std::size_t k = e + 1 + (tok[e + 1] == '+' || down); k < tok.size(); ++k)
        ex = std::min(10 * ex + (tok[k] - '0'), 100000L);
      if (down) ex = -ex;
    }
    // the power of ten of the leading digit other than 0
    const auto point = mantissa.find('.'), first = mantissa.find_first_of("123456789");
    const long whole = static_cast<long>(point == std::string::npos ? mantissa.size() : point);
    const long at = static_cast<long>(first == std::string::npos ? 0 : first);
    const long lead = ex + (at < whole ? whole - 1 - at : whole - at);
    v = lead > 0 && first != std::string::npos ? std::numeric_limits<double>::infinity() : 0.0;
    if (v == 0 && first != std::string::npos && lead > -400) {  // anything smaller is zero in double precision
      std::istringstream js(mantissa + "e" + std::to_string(ex + 300));
      js.imbue(std::locale::classic());
      double w = 0;
      js >> w;
      if (!js.fail()) v = w * 1e-300;
    }
  }
  out = neg ? -v : v;
  pos = i;
  return true;
}

bool parse_double(const std::string& s, double& out) {
  std::size_t i = 0;
  while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
  if (!parse_double(s, i, out)) return false;
  while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
  return i == s.size();
}

std::string format_g(double d, int precision) {
  std::ostringstream os;
  os.imbue(std::locale::classic());
  os << std::setprecision(precision) << d;
  return os.str();
}

std::string format_fixed(double d, int decimals) {
  std::ostringstream os;
  os.imbue(std::locale::classic());
  os << std::fixed << std::setprecision(decimals) << d;
  std::string r = os.str();
  if (r == "-0" || (r.size() > 1 && r[0] == '-' && r.find_first_not_of("-0.") == std::string::npos)) r.erase(0, 1);
  return r;
}

std::string format_roundtrip(double d) {
  if (!std::isfinite(d)) return std::isnan(d) ? "nan" : (d > 0 ? "inf" : "-inf");
  for (int p = 15; p <= 17; ++p) {
    const std::string s = format_g(d, p);
    double back;
    if (parse_double(s, back) && back == d) return s;
  }
  return format_g(d, 17);
}

}  // namespace ssc::text
