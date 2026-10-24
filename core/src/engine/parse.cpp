#include "engine/parse.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <locale>
#include <sstream>
#include <utility>

#include "util/numtext.hpp"

namespace ssc::parse {

namespace {

// Whether a character is shown in a message as it is: printable ASCII and the letters, digits and symbols of common
// scripts.  Any other (controls, the format characters that reorder or join text, such as U+202E, which turns the
// rest of the line around, combining marks, which land on the quote mark, and the scripts this leaves out) is
// written as \u{...}.
bool shown_as_is(char32_t c) {
  // (the digits of other scripts come with the Arabic comma, decimal and thousands separators and percent sign)
  static const char32_t ranges[][2] = {
      {0x20, 0x7E},     {0xA0, 0xAC},     {0xAE, 0x2FF},    {0x370, 0x482},   {0x48A, 0x52F},
      {0x60C, 0x60C},   {0x660, 0x66C},   {0x6F0, 0x6F9},   {0x966, 0x96F},   {0x9E6, 0x9EF},
      {0xA66, 0xA6F},   {0xAE6, 0xAEF},   {0xB66, 0xB6F},   {0xBE6, 0xBEF},   {0xC66, 0xC6F},
      {0xCE6, 0xCEF},   {0xD66, 0xD6F},   {0xE50, 0xE59},   {0xED0, 0xED9},   {0xF20, 0xF29},
      {0x1040, 0x1049}, {0x17E0, 0x17E9}, {0x1810, 0x1819}, {0x1E00, 0x1FFF}, {0x2010, 0x2027},
      {0x2030, 0x205E}, {0x2070, 0x20CF}, {0x2100, 0x27BF}, {0x3000, 0x3029}, {0x3030, 0x303F},
      {0x3041, 0x3096}, {0x309B, 0x30FF}, {0x4E00, 0x9FFF}, {0xAC00, 0xD7A3}, {0xFF01, 0xFF5E},
      {0xFF61, 0xFF9F}, {0xFFFD, 0xFFFD}, {0x1F300, 0x1FAFF}};  // U+FFFD: where text was not Unicode
  for (const auto& r : ranges)
    if (c >= r[0] && c <= r[1]) return true;
  return false;
}

// The character of the UTF-8 sequence at t[i] and its length; length 0 where the bytes are not UTF-8 (a lone
// continuation byte, an overlong form, a surrogate, a sequence cut short or past U+10FFFF).
std::pair<char32_t, std::size_t> decode(const std::string& t, std::size_t i) {
  const unsigned char b = static_cast<unsigned char>(t[i]);
  if (b < 0x80) return {b, 1};
  const std::size_t len = b >= 0xF0 ? 4 : b >= 0xE0 ? 3 : b >= 0xC2 ? 2 : 0;
  if (len == 0 || b > 0xF4 || i + len > t.size()) return {0, 0};
  char32_t c = b & (0x7F >> len);
  for (std::size_t k = 1; k < len; ++k) {
    const unsigned char x = static_cast<unsigned char>(t[i + k]);
    if ((x & 0xC0) != 0x80) return {0, 0};
    c = (c << 6) | (x & 0x3F);
  }
  static const char32_t least[] = {0, 0, 0x80, 0x800, 0x10000};
  if (c < least[len] || c > 0x10FFFF || (c >= 0xD800 && c < 0xE000)) return {0, 0};
  return {c, len};
}

}  // namespace

std::string quoted(const std::string& t) {
  std::string out;
  for (std::size_t i = 0; i < t.size();) {
    const auto [c, len] = decode(t, i);
    std::string piece;
    if (len == 0) {
      piece = t.substr(i, 1);  // not UTF-8: the byte, which the JSON writer shows as U+FFFD
    } else if (shown_as_is(c)) {
      piece = t.substr(i, len);
    } else {
      char hex[16];
      std::snprintf(hex, sizeof hex, "\\u{%X}", static_cast<unsigned>(c));
      piece = hex;
    }
    if (out.size() + piece.size() > 40) return out + "…";  // 40 bytes of the message, written out
    out += piece;
    i += len ? len : 1;
  }
  return out;
}

std::size_t characters(const std::string& t) {
  std::size_t n = 0;
  for (char ch : t) n += (static_cast<unsigned char>(ch) & 0xC0) != 0x80;
  return n;
}

namespace {

// Tokens: numbers and the keywords "to", "by"/"step".
struct Token {
  enum Kind { Num, To, By, Bad } kind;
  double value = 0;
  std::string text;
};

// Text typed or pasted with Unicode spacing and minus signs, written plainly: non-breaking, narrow and other
// spaces become spaces, and the minus sign (U+2212) a hyphen-minus.  (An en dash is left alone: "10–20" may mean
// a range.)
std::string plain(const std::string& s) {
  std::string out;
  out.reserve(s.size());
  for (std::size_t i = 0; i < s.size();) {
    const auto b = [&](std::size_t k) { return k < s.size() ? static_cast<unsigned char>(s[k]) : 0u; };
    if (b(i) == 0xC2 && b(i + 1) == 0xA0) {  // no-break space
      out += ' ', i += 2;
    } else if (b(i) == 0xE2 && b(i + 1) == 0x80 && ((b(i + 2) >= 0x80 && b(i + 2) <= 0x8A) || b(i + 2) == 0xAF)) {
      out += ' ', i += 3;  // en/em/thin/hair spaces, narrow no-break space
    } else if (b(i) == 0xE2 && ((b(i + 1) == 0x80 && b(i + 2) >= 0x8B && b(i + 2) <= 0x8D) ||
                                (b(i + 1) == 0x81 && b(i + 2) == 0xA0))) {
      i += 3;  // zero-width space, non-joiner, joiner, word joiner: unseen, so read as nothing ("0<ZWSP>.4" is 0.4)
    } else if (b(i) == 0xE2 && b(i + 1) == 0x81 && b(i + 2) == 0x9F) {  // medium mathematical space
      out += ' ', i += 3;
    } else if (b(i) == 0xE3 && b(i + 1) == 0x80 && b(i + 2) == 0x80) {  // ideographic space
      out += ' ', i += 3;
    } else if (b(i) == 0xEF && b(i + 1) == 0xBB && b(i + 2) == 0xBF) {  // byte-order mark
      i += 3;
    } else if (b(i) == 0xE2 && b(i + 1) == 0x88 && b(i + 2) == 0x92) {  // minus sign
      out += '-', i += 3;
    } else {
      out += s[i++];
    }
  }
  return out;
}

std::vector<Token> tokenize(const std::string& s, std::size_t max_tokens, std::string& err) {
  std::vector<Token> out;
  std::size_t i = 0;
  while (i < s.size()) {
    if (out.size() > max_tokens) {
      err = "too many values";
      return {};
    }
    const unsigned char c = static_cast<unsigned char>(s[i]);
    // a comma between two digits ("2,500", "1,5") may be a thousands separator, a decimal comma or a list:
    // rather than guess (a list would quietly make two scenarios), ask for the number written plainly
    if (c == ',' && i > 0 && i + 1 < s.size() && std::isdigit(static_cast<unsigned char>(s[i - 1])) &&
        std::isdigit(static_cast<unsigned char>(s[i + 1]))) {
      std::size_t a = i, b = i + 1;
      while (a > 0 && !std::isspace(static_cast<unsigned char>(s[a - 1]))) --a;
      while (b < s.size() && !std::isspace(static_cast<unsigned char>(s[b]))) ++b;
      err = "could not read '" + quoted(s.substr(a, b - a)) +
            "': write numbers without commas (2500, or 2.5 with a decimal point) and separate several values "
            "with spaces";
      return {};
    }
    if (std::isspace(c) || c == ',') {
      ++i;
      continue;
    }
    if (std::isalpha(c)) {
      std::size_t j = i;
      while (j < s.size() && std::isalpha(static_cast<unsigned char>(s[j]))) ++j;
      std::string w = s.substr(i, j - i);
      for (auto& ch : w) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
      if (w == "to") {
        out.push_back({Token::To, 0, w});
        i = j;
        continue;
      }
      if (w == "by" || w == "step") {
        out.push_back({Token::By, 0, w});
        i = j;
        continue;
      }
      // could be inf / nan
    }
    std::size_t j = i;
    double v = 0;
    if (text::parse_double(s, j, v) && (j == s.size() || !std::isalnum(static_cast<unsigned char>(s[j])) ||
                                        s[j] == ',')) {
      // "10-20" is two numbers, 10 and -20, but was most likely meant as a range: ask for it plainly
      if (j + 1 < s.size() && (s[j] == '-' || s[j] == '+') &&
          (std::isdigit(static_cast<unsigned char>(s[j + 1])) || s[j + 1] == '.')) {
        std::size_t k = j + 1;
        while (k < s.size() && !std::isspace(static_cast<unsigned char>(s[k])) && s[k] != ',') ++k;
        err = s[j] == '-' ? "could not read '" + quoted(s.substr(i, k - i)) + "': write a range as '" +
                                quoted(s.substr(i, j - i)) + " to " + quoted(s.substr(j + 1, k - j - 1)) +
                                "', and separate several values with spaces"
                          : "could not read '" + quoted(s.substr(i, k - i)) +
                                "' as a number (separate several values with spaces)";
        return {};
      }
      // a number too small for a double is not 0: it is read as the smallest double of its sign, so that it is
      // refused for its size as 1e-31 is, not taken for 0 ('1e-400' was)
      if (v == 0) {
        bool digit = false;  // other than 0, before any exponent
        for (std::size_t k = i; k < j && s[k] != 'e' && s[k] != 'E'; ++k) digit = digit || (s[k] >= '1' && s[k] <= '9');
        if (digit) v = std::copysign(std::numeric_limits<double>::denorm_min(), v);
      }
      if (j < s.size() && s[j] == '%') {
        err = "could not read '" + quoted(s.substr(i, j - i + 1)) +
              "': write the number without '%' (a probability as a decimal, such as 0.8)";
        return {};
      }
      out.push_back({Token::Num, v, s.substr(i, j - i)});
      i = j;
      continue;
    }
    // unreadable chunk
    std::size_t k = i;
    while (k < s.size() && !std::isspace(static_cast<unsigned char>(s[k])) && s[k] != ',') ++k;
    err = "could not read '" + quoted(s.substr(i, k - i)) + "' as a number";
    return {};
  }
  return out;
}

// A range's values are worked out from the decimals typed, exactly, and each is then rounded to a double once: in
// doubles, a + k * step leaves rounding noise (0.1 + 2 * 0.1 is 0.30000000000000004, -0.3 + 3 * 0.1 is 5.6e-17,
// -0.4999999999 + 2 * 0.25 is 1.00000008e-10), and rounding it to fewer digits loses digits typed.

// A number as typed: (-1)^neg * digits * 10^exp, its digits without leading or trailing zeros (none for 0)
struct Decimal {
  bool neg = false;
  std::string digits;
  long exp = 0;
};

// The decimal of a number token v (text::parse_double's syntax: a sign, digits with a point, an exponent).  One with
// more significant digits than a double holds, or an exponent no double has, is taken as v's 17 digits.
Decimal decimal(const std::string& t, double v) {
  Decimal d;
  std::size_t i = 0;
  if (i < t.size() && (t[i] == '+' || t[i] == '-')) d.neg = t[i++] == '-';
  for (bool point = false; i < t.size() && (std::isdigit(static_cast<unsigned char>(t[i])) || (t[i] == '.' && !point));
       ++i) {
    if (t[i] == '.') {
      point = true;
    } else {
      d.digits += t[i];
      if (point) --d.exp;
    }
  }
  if (i < t.size() && (t[i] == 'e' || t[i] == 'E')) {
    const bool down = ++i < t.size() && t[i] == '-';
    if (i < t.size() && (t[i] == '+' || t[i] == '-')) ++i;
    long e = 0;
    for (; i < t.size() && std::isdigit(static_cast<unsigned char>(t[i])); ++i)
      e = std::min(10 * e + (t[i] - '0'), 1000000L);  // far past any double's
    d.exp += down ? -e : e;
  }
  const std::size_t first = d.digits.find_first_not_of('0');
  if (first == std::string::npos) return {};
  const std::size_t last = d.digits.find_last_not_of('0');
  d.exp += static_cast<long>(d.digits.size() - 1 - last);
  d.digits = d.digits.substr(first, last + 1 - first);
  if (d.digits.size() > 17 || d.exp < -400 || d.exp > 400) return decimal(text::format_g(v, 17), v);
  return d;
}

// A whole number: its sign and its digits in base 10^9, least significant first, without leading zeros (none for 0,
// which is not negative).  Sums of numbers typed with exponents far apart (1e-300 and 1e300) have hundreds of digits.
struct Whole {
  bool neg = false;
  std::vector<std::uint32_t> limbs;
};

constexpr std::uint32_t BASE = 1000000000;

// The whole number written with the decimal digits t
std::vector<std::uint32_t> limbs_of(const std::string& t) {
  std::vector<std::uint32_t> out;
  for (std::size_t end = t.size(); end > 0;) {
    const std::size_t begin = end > 9 ? end - 9 : 0;
    std::uint32_t limb = 0;
    for (std::size_t k = begin; k < end; ++k) limb = 10 * limb + static_cast<std::uint32_t>(t[k] - '0');
    out.push_back(limb);
    end = begin;
  }
  while (!out.empty() && out.back() == 0) out.pop_back();
  return out;
}

int compare_limbs(const std::vector<std::uint32_t>& x, const std::vector<std::uint32_t>& y) {
  if (x.size() != y.size()) return x.size() < y.size() ? -1 : 1;
  for (std::size_t k = x.size(); k-- > 0;)
    if (x[k] != y[k]) return x[k] < y[k] ? -1 : 1;
  return 0;
}

// x += y, of magnitudes
void add_limbs(std::vector<std::uint32_t>& x, const std::vector<std::uint32_t>& y) {
  if (x.size() < y.size()) x.resize(y.size(), 0);
  std::uint32_t carry = 0;
  for (std::size_t k = 0; k < x.size() && (k < y.size() || carry); ++k) {
    const std::uint32_t d = x[k] + (k < y.size() ? y[k] : 0) + carry;  // below 2^31
    carry = d >= BASE;
    x[k] = carry ? d - BASE : d;
  }
  if (carry) x.push_back(1);
}

// x = x - y where x > y, or x = y - x where y > x (`reverse`), of magnitudes
void subtract_limbs(std::vector<std::uint32_t>& x, const std::vector<std::uint32_t>& y, bool reverse) {
  if (x.size() < y.size()) x.resize(y.size(), 0);
  std::uint32_t borrow = 0;
  for (std::size_t k = 0; k < x.size() && (reverse || k < y.size() || borrow); ++k) {
    const std::uint32_t yk = k < y.size() ? y[k] : 0, from = reverse ? yk : x[k], take = (reverse ? x[k] : yk) + borrow;
    borrow = from < take;
    x[k] = borrow ? from + BASE - take : from - take;
  }
  while (!x.empty() && x.back() == 0) x.pop_back();
}

// x += y
void plus(Whole& x, const Whole& y) {
  if (y.limbs.empty()) return;
  if (x.limbs.empty() || x.neg == y.neg) {
    x.neg = y.neg;
    add_limbs(x.limbs, y.limbs);
    return;
  }
  const int c = compare_limbs(x.limbs, y.limbs);
  if (c == 0) {
    x = Whole{};
  } else {
    subtract_limbs(x.limbs, y.limbs, c < 0);
    if (c < 0) x.neg = y.neg;
  }
}

Whole negative(Whole w) {
  if (!w.limbs.empty()) w.neg = !w.neg;
  return w;
}

int compare(const Whole& x, const Whole& y) {
  if (x.neg != y.neg) return x.neg ? -1 : 1;
  const int c = compare_limbs(x.limbs, y.limbs);
  return x.neg ? -c : c;
}

// x * m, of a magnitude and a factor below 2^32
std::vector<std::uint32_t> times(const std::vector<std::uint32_t>& x, std::uint64_t m) {
  std::vector<std::uint32_t> out;
  std::uint64_t carry = 0;
  for (std::size_t k = 0; k < x.size() || carry; ++k) {
    carry += (k < x.size() ? x[k] : 0) * m;  // below 2^62 + 2^33
    out.push_back(static_cast<std::uint32_t>(carry % BASE));
    carry /= BASE;
  }
  while (!out.empty() && out.back() == 0) out.pop_back();
  return out;
}

// The decimal digits of the limbs of w from its most significant down to limb `last`
std::string digits_of(const std::vector<std::uint32_t>& w, std::size_t last) {
  std::string out = std::to_string(w.back());
  for (std::size_t k = w.size() - 1; k-- > last;) {
    char d[9];
    for (std::uint32_t j = 9, x = w[k]; j-- > 0; x /= 10) d[j] = static_cast<char>('0' + x % 10);
    out.append(d, 9);
  }
  return out;
}

// Decimal text read as the double nearest it, as text::parse_double reads it, from one stream for all the values of a
// range: making a stream takes longer than reading a number from it
class Reader {
 public:
  Reader() { in_.imbue(std::locale::classic()); }
  double operator()(const std::string& t) {
    in_.clear();
    in_.str(t);
    double v = 0;
    if (!(in_ >> v)) text::parse_double(t, v);  // out of the doubles' range, which it handles
    return v;
  }

 private:
  std::istringstream in_;
};

// The double nearest w * 10^e.  One too small for a double is the smallest double of its sign, not 0, as a number
// typed so is.
double nearest(const Whole& w, long e, Reader& read) {
  if (w.limbs.empty()) return 0;
  const std::size_t n = w.limbs.size();
  double v = 0;
  if ((n == 1 || (n == 2 && w.limbs[1] < 1000000)) && e >= -22 && e <= 22) {
    // a whole number below 2^53 and a power of ten up to 1e22 are doubles exactly: one product or quotient rounds once
    static const double ten[] = {1e0,  1e1,  1e2,  1e3,  1e4,  1e5,  1e6,  1e7,  1e8,  1e9,  1e10, 1e11,
                                 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22};
    v = (n == 2 ? w.limbs[1] * 1e9 : 0) + w.limbs[0];
    v = e >= 0 ? v * ten[e] : v / ten[-e];
  } else if (n < 12) {
    v = read(digits_of(w.limbs, 0) + "e" + std::to_string(e));
  } else {
    // 100 digits or more, which take microseconds to read: its first 19 digits t, from the three leading limbs.
    // Where every digit after them is 0, w is t times the power of ten of t's last digit; where one is not, w lies
    // between that and t + 1 times it, and where those two round to the same double, so does w.  Only where they do
    // not is every digit read.
    const std::string lead = digits_of(w.limbs, n - 3), t = lead.substr(0, 19);
    const long last = e + 9 * static_cast<long>(n - 3) + static_cast<long>(lead.size()) - 19;
    const std::string exp = "e" + std::to_string(last);
    bool rest = lead.find_first_not_of('0', 19) != std::string::npos;
    for (std::size_t k = 0; k < n - 3 && !rest; ++k) rest = w.limbs[k] != 0;
    v = read(t + exp);
    if (rest && read(std::to_string(std::stoull(t) + 1) + exp) != v)
      v = read(digits_of(w.limbs, 0) + "e" + std::to_string(e));
  }
  if (v == 0) v = std::numeric_limits<double>::denorm_min();
  return w.neg ? -v : v;
}

}  // namespace

Numbers numbers(const std::string& raw, std::size_t max_values) {
  Numbers res;
  // read before the value count could be checked, it took memory and seconds
  if (raw.size() > MAX_TEXT && (raw.size() > 4 * MAX_TEXT || characters(raw) > MAX_TEXT)) {
    res.error = "the text is too long (limit 300,000 characters)";
    return res;
  }
  const std::string text = plain(raw);
  std::string err;
  // a value takes at most 5 tokens ("a to b by c"); more than that many per value cannot be within the limit
  const auto toks = tokenize(text, 5 * max_values + 5, err);
  if (!err.empty()) {
    res.error = err == "too many values" ? "too many values (limit " + std::to_string(max_values) + ")" : err;
    return res;
  }
  for (std::size_t i = 0; i < toks.size();) {
    if (toks[i].kind != Token::Num) {
      res.error = "unexpected '" + quoted(toks[i].text) + "'";
      return res;
    }
    if (i + 1 < toks.size() && toks[i + 1].kind == Token::To) {
      if (i + 2 >= toks.size() || toks[i + 2].kind != Token::Num) {
        res.error = "a range needs the form 'a to b by c'";
        return res;
      }
      const double a = toks[i].value, b = toks[i + 2].value;
      double step = 1;
      std::size_t used = 3;
      if (i + 3 < toks.size() && toks[i + 3].kind == Token::By) {
        if (i + 4 >= toks.size() || toks[i + 4].kind != Token::Num) {
          res.error = "missing step after 'by'";
          return res;
        }
        step = toks[i + 4].value;
        used = 5;
      }
      if (!(std::isfinite(a) && std::isfinite(b) && std::isfinite(step))) {
        res.error = "range limits must be finite";
        return res;
      }
      if (step == 0) {
        res.error = "range step must not be zero";
        return res;
      }
      // a + k * step from the decimals typed, so the start is the value typed, so is the end where the steps reach
      // it (each is checked as it would be alone), and a range that crosses 0 has 0 among its values.  The steps go
      // toward the end whatever the sign of the step ("10 to 1 by 1" descends), and the count allows a step past the
      // end by a billionth of a step ('0 to 0.6666666666 by 0.11111111111' has 7 values): that value is the end.
      // The count is worked out from the decimals too: in doubles, the rounding of the ends shifts (b - a) / step
      // by more than that billionth once the start is a million steps or more ('1.96 to 1.9600003 by 0.00000001'
      // lost its end).
      const Decimal da = decimal(toks[i].text, a), db = decimal(toks[i + 2].text, b);
      const Decimal ds = used == 5 ? decimal(toks[i + 4].text, toks[i + 4].value) : Decimal{false, "1", 0};
      long e = ds.exp;  // the exponent of the least significant digit typed (the step is not 0)
      for (const Decimal* d : {&da, &db})
        if (!d->digits.empty()) e = std::min(e, d->exp);
      const auto whole = [&](const Decimal& d) {
        if (d.digits.empty()) return Whole{};
        return Whole{d.neg, limbs_of(d.digits + std::string(static_cast<std::size_t>(d.exp - e), '0'))};
      };
      const Whole start = whole(da), end = whole(db);
      Whole length = end, s = whole(ds);
      plus(length, negative(start));
      s.neg = length.neg;
      // the steps k with k * step <= (end - start) + step / 10^9, in whole numbers k * step * 10^9 <= (end - start)
      // * 10^9 + step, are those with k * step <= (end - start) + step without its last limb (10^9 is a limb)
      std::vector<std::uint32_t> reach = length.limbs;
      if (!s.limbs.empty()) add_limbs(reach, std::vector<std::uint32_t>(s.limbs.begin() + 1, s.limbs.end()));
      const std::size_t room = max_values - res.values.size();
      if (compare_limbs(times(s.limbs, room), reach) <= 0) {  // room + 1 values or more
        res.error = "too many values (limit " + std::to_string(max_values) + ")";
        return res;
      }
      std::size_t steps = 0, over = room;  // steps reach no further than the end allows, over does
      while (over - steps > 1) {
        const std::size_t mid = steps + (over - steps) / 2;
        if (compare_limbs(times(s.limbs, mid), reach) <= 0) {
          steps = mid;
        } else {
          over = mid;
        }
      }
      res.ends.push_back(a);
      res.ends.push_back(b);
      Whole v = start;
      Reader read;
      for (std::size_t k = 0; k <= steps; ++k) {
        if (k > 0) plus(v, s);
        res.values.push_back(k == 0 ? a : compare(v, end) * (length.neg ? -1 : 1) >= 0 ? b : nearest(v, e, read));
      }
      i += used;
      continue;
    }
    res.values.push_back(toks[i].value);
    if (res.values.size() > max_values) {
      res.error = "too many values (limit " + std::to_string(max_values) + ")";
      return res;
    }
    ++i;
  }
  if (res.values.empty()) res.error = "a value is required";
  return res;
}

Vectors vectors(const std::string& text, std::size_t max_vectors) {
  Vectors res;
  if (text.size() > MAX_TEXT && (text.size() > 4 * MAX_TEXT || characters(text) > MAX_TEXT)) {
    res.error = "the text is too long (limit 300,000 characters)";
    return res;
  }
  std::size_t start = 0, total = 0;
  while (start <= text.size()) {
    std::size_t end = text.find_first_of(";|", start);
    if (end == std::string::npos) end = text.size();
    const std::string part = text.substr(start, end - start);
    bool blank = true;
    for (char c : part)
      if (!std::isspace(static_cast<unsigned char>(c))) blank = false;
    if (!blank) {
      Numbers n = numbers(part);
      if (!n.error.empty()) {
        res.error = n.error;
        return res;
      }
      // the limit of a list of numbers holds for all the alternatives together: each becomes scenarios of its own
      total += n.values.size();
      if (total > 10000) {
        res.error = "too many values in all the lists together (limit 10000)";
        return res;
      }
      res.vectors.push_back(std::move(n.values));
      res.ends.insert(res.ends.end(), n.ends.begin(), n.ends.end());
      if (res.vectors.size() > max_vectors) {
        res.error = "too many alternatives";
        return res;
      }
    }
    if (end == text.size()) break;
    start = end + 1;
  }
  if (res.vectors.empty()) res.error = "a list of values is required";
  return res;
}

}  // namespace ssc::parse
