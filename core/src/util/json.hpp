// Minimal JSON value, parser and serializer (UTF-8, insertion-ordered objects).
#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ssc::json {

class Value {
 public:
  enum class Type { Null, Bool, Number, String, Array, Object };

  Value() = default;
  Value(std::nullptr_t) {}
  Value(bool b) : t_(Type::Bool), b_(b) {}
  Value(double d) : t_(Type::Number), n_(d) {}
  Value(int i) : t_(Type::Number), n_(i) {}
  Value(long i) : t_(Type::Number), n_(static_cast<double>(i)) {}
  Value(long long i) : t_(Type::Number), n_(static_cast<double>(i)) {}
  Value(unsigned long i) : t_(Type::Number), n_(static_cast<double>(i)) {}
  Value(unsigned long long i) : t_(Type::Number), n_(static_cast<double>(i)) {}
  Value(const char* s) : t_(Type::String), s_(s) {}
  Value(std::string s) : t_(Type::String), s_(std::move(s)) {}

  static Value array() {
    Value v;
    v.t_ = Type::Array;
    return v;
  }
  static Value object() {
    Value v;
    v.t_ = Type::Object;
    return v;
  }

  Type type() const { return t_; }
  bool is_null() const { return t_ == Type::Null; }
  bool is_bool() const { return t_ == Type::Bool; }
  bool is_number() const { return t_ == Type::Number; }
  bool is_string() const { return t_ == Type::String; }
  bool is_array() const { return t_ == Type::Array; }
  bool is_object() const { return t_ == Type::Object; }

  bool as_bool(bool def = false) const { return t_ == Type::Bool ? b_ : def; }
  double as_number(double def = 0) const { return t_ == Type::Number ? n_ : def; }
  const std::string& as_string() const;
  /// String content, or the JSON text of numbers/booleans (inputs may be either).
  std::string to_text() const;

  // arrays
  Value& push(Value v);
  std::size_t size() const { return t_ == Type::Array ? a_.size() : (t_ == Type::Object ? o_.size() : 0); }
  const Value& at(std::size_t i) const { return a_.at(i); }
  const std::vector<Value>& elements() const { return a_; }

  // objects (insertion ordered; set() replaces an existing key)
  Value& set(const std::string& key, Value v);
  const Value* find(const std::string& key) const;
  const Value& get(const std::string& key) const;  // null value if absent
  const std::vector<std::pair<std::string, Value>>& items() const { return o_; }

  std::string dump() const;
  void dump_to(std::string& out) const;

 private:
  Type t_ = Type::Null;
  bool b_ = false;
  double n_ = 0;
  std::string s_;
  std::vector<Value> a_;
  std::vector<std::pair<std::string, Value>> o_;
};

/// Parse JSON text; throws std::runtime_error with a position on failure.
Value parse(const std::string& text);

/// Append a JSON-escaped string literal (with quotes).
void escape_string(const std::string& s, std::string& out);
/// Shortest round-trip representation of a double (null for non-finite).
void append_number(double d, std::string& out);

}  // namespace ssc::json
