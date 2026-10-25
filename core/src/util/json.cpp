#include "json.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <unordered_map>

#include "numtext.hpp"

namespace ssc::json {

const std::string& Value::as_string() const {
  static const std::string empty;
  return t_ == Type::String ? s_ : empty;
}

std::string Value::to_text() const {
  switch (t_) {
    case Type::String:
      return s_;
    case Type::Number: {
      if (!std::isfinite(n_)) return text::format_roundtrip(n_);  // "inf", which an input refuses as not finite
      std::string s;
      append_number(n_, s);
      return s;
    }
    case Type::Bool:
      return b_ ? "true" : "false";
    default:
      return std::string();
  }
}

Value& Value::push(Value v) {
  if (t_ != Type::Array) {
    t_ = Type::Array;
    a_.clear();
  }
  a_.push_back(std::move(v));
  return a_.back();
}

Value& Value::set(const std::string& key, Value v) {
  if (t_ != Type::Object) {
    t_ = Type::Object;
    o_.clear();
  }
  for (auto& kv : o_)
    if (kv.first == key) {
      kv.second = std::move(v);
      return kv.second;
    }
  o_.emplace_back(key, std::move(v));
  return o_.back().second;
}

Value& Value::append(std::string key, Value v) {
  if (t_ != Type::Object) {
    t_ = Type::Object;
    o_.clear();
  }
  o_.emplace_back(std::move(key), std::move(v));
  return o_.back().second;
}

const Value* Value::find(const std::string& key) const {
  if (t_ != Type::Object) return nullptr;
  for (const auto& kv : o_)
    if (kv.first == key) return &kv.second;
  return nullptr;
}

const Value& Value::get(const std::string& key) const {
  static const Value null_value;
  const Value* v = find(key);
  return v ? *v : null_value;
}

namespace {
// The length of the valid UTF-8 sequence starting at s[i], or 0 if it is not one (a stray byte, an overlong
// form, a surrogate, a code point beyond U+10FFFF or a sequence cut short).
std::size_t utf8_length(const std::string& s, std::size_t i) {
  const auto b = [&](std::size_t k) { return static_cast<unsigned char>(s[k]); };
  const unsigned char c = b(i);
  std::size_t n;
  unsigned cp;
  if (c >= 0xC2 && c <= 0xDF) { n = 2; cp = c & 0x1Fu; }
  else if (c >= 0xE0 && c <= 0xEF) { n = 3; cp = c & 0x0Fu; }
  else if (c >= 0xF0 && c <= 0xF4) { n = 4; cp = c & 0x07u; }
  else return 0;
  if (i + n > s.size()) return 0;
  for (std::size_t k = 1; k < n; ++k) {
    if ((b(i + k) & 0xC0) != 0x80) return 0;
    cp = (cp << 6) | (b(i + k) & 0x3Fu);
  }
  if ((n == 3 && (cp < 0x800 || (cp >= 0xD800 && cp < 0xE000))) || (n == 4 && (cp < 0x10000 || cp > 0x10FFFF)))
    return 0;
  return n;
}
}  // namespace

void escape_string(const std::string& s, std::string& out) {
  out.push_back('"');
  for (std::size_t i = 0; i < s.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(s[i]);
    if (c >= 0x80) {  // the output is always valid UTF-8: text that is not becomes U+FFFD
      const std::size_t n = utf8_length(s, i);
      if (n == 0) {
        out += "\\ufffd";
      } else {
        out.append(s, i, n);
        i += n - 1;
      }
      continue;
    }
    switch (c) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      case '\b': out += "\\b"; break;
      case '\f': out += "\\f"; break;
      default:
        if (c < 0x20) {
          char buf[8];
          std::snprintf(buf, sizeof buf, "\\u%04x", c);
          out += buf;
        } else {
          out.push_back(static_cast<char>(c));
        }
    }
  }
  out.push_back('"');
}

void append_number(double d, std::string& out) {
  if (!std::isfinite(d)) {
    out += "null";
    return;
  }
  if (d == std::floor(d) && std::fabs(d) < 1e15) {
    out += text::format_fixed(d, 0);
    return;
  }
  out += text::format_roundtrip(d);
}

void Value::dump_to(std::string& out) const {
  switch (t_) {
    case Type::Null: out += "null"; break;
    case Type::Bool: out += b_ ? "true" : "false"; break;
    case Type::Number: append_number(n_, out); break;
    case Type::String: escape_string(s_, out); break;
    case Type::Array: {
      out.push_back('[');
      bool first = true;
      for (const auto& v : a_) {
        if (!first) out.push_back(',');
        first = false;
        v.dump_to(out);
      }
      out.push_back(']');
      break;
    }
    case Type::Object: {
      out.push_back('{');
      bool first = true;
      for (const auto& kv : o_) {
        if (!first) out.push_back(',');
        first = false;
        escape_string(kv.first, out);
        out.push_back(':');
        kv.second.dump_to(out);
      }
      out.push_back('}');
      break;
    }
  }
}

std::string Value::dump() const {
  std::string s;
  dump_to(s);
  return s;
}

// ------------------------------------------------------------------ parser
namespace {
struct Parser {
  // Objects and arrays nested deeper than this are refused: a request is three levels deep, and each level is a
  // recursive call (100,000 levels overflowed the stack).
  static constexpr int MAX_DEPTH = 64;

  const std::string& s;
  std::size_t i = 0;
  int depth = 0;

  struct Nest {
    int& depth;
    ~Nest() { --depth; }
  };
  Nest nest() {
    if (++depth > MAX_DEPTH) fail("nested too deeply");
    return Nest{depth};
  }

  [[noreturn]] void fail(const char* what) const {
    throw std::runtime_error(std::string("JSON parse error at offset ") + std::to_string(i) + ": " + what);
  }
  void ws() {
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\n' || s[i] == '\r')) ++i;
  }
  bool lit(const char* w) {
    const std::size_t n = std::strlen(w);
    if (s.compare(i, n, w) == 0) {
      i += n;
      return true;
    }
    return false;
  }
  static void utf8(unsigned cp, std::string& out) {
    if (cp < 0x80) {
      out.push_back(static_cast<char>(cp));
    } else if (cp < 0x800) {
      out.push_back(static_cast<char>(0xC0 | (cp >> 6)));
      out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    } else if (cp < 0x10000) {
      out.push_back(static_cast<char>(0xE0 | (cp >> 12)));
      out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
      out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    } else {
      out.push_back(static_cast<char>(0xF0 | (cp >> 18)));
      out.push_back(static_cast<char>(0x80 | ((cp >> 12) & 0x3F)));
      out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
      out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    }
  }
  unsigned hex4() {
    if (i + 4 > s.size()) fail("bad \\u escape");
    unsigned v = 0;
    for (int k = 0; k < 4; ++k) {
      const char c = s[i++];
      v <<= 4;
      if (c >= '0' && c <= '9') v |= static_cast<unsigned>(c - '0');
      else if (c >= 'a' && c <= 'f') v |= static_cast<unsigned>(c - 'a' + 10);
      else if (c >= 'A' && c <= 'F') v |= static_cast<unsigned>(c - 'A' + 10);
      else fail("bad hex digit");
    }
    return v;
  }
  std::string str() {
    if (s[i] != '"') fail("expected string");
    ++i;
    std::string out;
    while (true) {
      if (i >= s.size()) fail("unterminated string");
      const char c = s[i++];
      if (c == '"') break;
      if (c != '\\') {
        out.push_back(c);
        continue;
      }
      if (i >= s.size()) fail("bad escape");
      const char e = s[i++];
      switch (e) {
        case '"': out.push_back('"'); break;
        case '\\': out.push_back('\\'); break;
        case '/': out.push_back('/'); break;
        case 'b': out.push_back('\b'); break;
        case 'f': out.push_back('\f'); break;
        case 'n': out.push_back('\n'); break;
        case 'r': out.push_back('\r'); break;
        case 't': out.push_back('\t'); break;
        case 'u': {
          unsigned cp = hex4();
          if (cp >= 0xD800 && cp < 0xDC00) {  // a surrogate pair; half of one alone is U+FFFD
            const std::size_t at = i;
            unsigned lo = 0;
            if (i + 1 < s.size() && s[i] == '\\' && s[i + 1] == 'u') {
              i += 2;
              lo = hex4();
            }
            if (lo >= 0xDC00 && lo < 0xE000) {
              cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
            } else {
              i = at;
              cp = 0xFFFD;
            }
          } else if (cp >= 0xDC00 && cp < 0xE000) {
            cp = 0xFFFD;
          }
          utf8(cp, out);
          break;
        }
        default: fail("bad escape");
      }
    }
    return out;
  }
  Value value() {
    ws();
    if (i >= s.size()) fail("unexpected end");
    const char c = s[i];
    if (c == '{') {
      ++i;
      const Nest level = nest();
      Value o = Value::object();
      ws();
      if (i < s.size() && s[i] == '}') {
        ++i;
        return o;
      }
      // the keys seen, so that a long object is read in linear time (a repeated key keeps its first place and
      // takes the last value, as set() does)
      std::unordered_map<std::string, std::size_t> seen;
      while (true) {
        ws();
        std::string k = str();
        ws();
        if (i >= s.size() || s[i] != ':') fail("expected ':'");
        ++i;
        Value v = value();
        const auto at = seen.find(k);
        if (at != seen.end()) {
          o.member(at->second) = std::move(v);
        } else {
          seen.emplace(k, o.size());
          o.append(std::move(k), std::move(v));
        }
        ws();
        if (i < s.size() && s[i] == ',') { ++i; continue; }
        if (i < s.size() && s[i] == '}') { ++i; break; }
        fail("expected ',' or '}'");
      }
      return o;
    }
    if (c == '[') {
      ++i;
      const Nest level = nest();
      Value a = Value::array();
      ws();
      if (i < s.size() && s[i] == ']') {
        ++i;
        return a;
      }
      while (true) {
        a.push(value());
        ws();
        if (i < s.size() && s[i] == ',') { ++i; continue; }
        if (i < s.size() && s[i] == ']') { ++i; break; }
        fail("expected ',' or ']'");
      }
      return a;
    }
    if (c == '"') return Value(str());
    if (lit("true")) return Value(true);
    if (lit("false")) return Value(false);
    if (lit("null")) return Value();
    // number (locale independent).  One too small for a double whose digits are not all 0 is the smallest double of
    // its sign, as in typed text (parse::numbers), so that an input sent as 1e-400 is refused for its size, not
    // taken for 0
    const std::size_t start = i;
    double d = 0;
    if (!text::parse_double(s, i, d)) fail("unexpected character");
    for (std::size_t k = start; d == 0 && k < i && s[k] != 'e' && s[k] != 'E'; ++k)
      if (s[k] >= '1' && s[k] <= '9') d = std::copysign(std::numeric_limits<double>::denorm_min(), d);
    return Value(d);
  }
};
}  // namespace

Value parse(const std::string& text) {
  Parser p{text};
  Value v = p.value();
  p.ws();
  if (p.i != text.size()) p.fail("trailing characters");
  return v;
}

}  // namespace ssc::json
