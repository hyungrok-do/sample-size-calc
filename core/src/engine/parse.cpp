#include "engine/parse.hpp"

#include <cctype>
#include <cmath>
#include <cstdio>

#include "util/numtext.hpp"

namespace ssc::parse {

namespace {

// A character that would act on a message rather than show in it: controls, the format characters that reorder or
// join text (U+202E turns the rest of the line around) and combining marks (they land on the quote mark).
bool acts_on_text(char32_t c) {
  return c < 0x20 || (c >= 0x7F && c < 0xA0) || c == 0xAD || c == 0x61C || (c >= 0x300 && c < 0x370) ||
         (c >= 0x1AB0 && c < 0x1B00) || (c >= 0x1DC0 && c < 0x1E00) || (c >= 0x200B && c < 0x2010) ||
         (c >= 0x2028 && c < 0x202F) || (c >= 0x2060 && c < 0x2070) || (c >= 0x20D0 && c < 0x2100) ||
         (c >= 0xFE00 && c < 0xFE10) || (c >= 0xFE20 && c < 0xFE30) || c == 0xFEFF || (c >= 0xFFF9 && c < 0xFFFC) ||
         (c >= 0xE0000 && c < 0xE1000);
}

}  // namespace

std::string quoted(const std::string& t) {
  std::size_t end = t.size();
  if (end > 40) {
    end = 40;
    while (end > 0 && (static_cast<unsigned char>(t[end]) & 0xC0) == 0x80) --end;
  }
  // such characters are written as \u{...}
  std::string out;
  for (std::size_t i = 0; i < end;) {
    const unsigned char b = static_cast<unsigned char>(t[i]);
    const std::size_t len = b < 0x80 ? 1 : b >= 0xF0 ? 4 : b >= 0xE0 ? 3 : b >= 0xC0 ? 2 : 1;
    char32_t c = len == 1 ? b : b & (0x7F >> len);
    bool whole = i + len <= end;
    for (std::size_t k = 1; whole && k < len; ++k) {
      const unsigned char x = static_cast<unsigned char>(t[i + k]);
      whole = (x & 0xC0) == 0x80;
      c = (c << 6) | (x & 0x3F);
    }
    if (whole && acts_on_text(c)) {
      char hex[16];
      std::snprintf(hex, sizeof hex, "\\u{%X}", static_cast<unsigned>(c));
      out += hex;
      i += len;
    } else if (whole) {
      out.append(t, i, len);
      i += len;
    } else {
      out += t[i++];  // not UTF-8: the byte as it is
    }
  }
  return end < t.size() ? out + "…" : out;
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
    } else if (b(i) == 0xE2 && b(i + 1) == 0x80 && ((b(i + 2) >= 0x80 && b(i + 2) <= 0x8B) || b(i + 2) == 0xAF)) {
      out += ' ', i += 3;  // en/em/thin/hair/zero-width spaces, narrow no-break space
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
    res.error = "the text is too long (limit 1,000,000 characters)";
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
      for (long k = 0; k < static_cast<long>(count); ++k) res.values.push_back(tidy(a + static_cast<double>(k) * step));
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
    res.error = "the text is too long (limit 1,000,000 characters)";
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
