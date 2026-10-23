#include "engine/parse.hpp"

#include <cctype>
#include <cmath>
#include <cstdio>
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

double tidy(double v) {
  // remove floating-point noise from generated range values (e.g. 0.30000000000000004)
  if (v == 0 || !std::isfinite(v)) return v;
  double r;
  return text::parse_double(text::format_g(v, 12), r) ? r : v;
}

}  // namespace

Numbers numbers(const std::string& raw, std::size_t max_values) {
  Numbers res;
  // read before the value count could be checked, it took memory and seconds
  if (raw.size() > MAX_TEXT && (raw.size() > 4 * MAX_TEXT || characters(raw) > MAX_TEXT)) {
    res.error = "the text is too long (limit 200,000 characters)";
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
      if ((b - a) * step < 0) step = -step;  // allow "10 to 1 by 1" (descending)
      const double count = std::floor((b - a) / step + 1e-9) + 1;
      if (count > static_cast<double>(max_values) || res.values.size() + static_cast<std::size_t>(count) > max_values) {
        res.error = "too many values (limit " + std::to_string(max_values) + ")";
        return res;
      }
      res.ends.push_back(a);
      res.ends.push_back(b);
      for (long k = 0; k < static_cast<long>(count); ++k) {
        // where a range crosses 0, the sum leaves rounding noise (-0.3 + 3 * 0.1 is 5.6e-17, -3e-15 + 3 * 1e-15 is
        // 3.9e-31): a value within a billionth of a step of 0 is 0 (the count above allows that much for rounding)
        const double v = a + static_cast<double>(k) * step;
        res.values.push_back(k > 0 && std::fabs(v) < 1e-9 * std::fabs(step) ? 0 : tidy(v));
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
    res.error = "the text is too long (limit 200,000 characters)";
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
