#include "program.hpp"
#include "program_algebra.hpp"
#include "program_runtime.hpp"
#include "program_sources.hpp"
#include "util/numtext.hpp"
#include <array>
#include <cctype>
#include <memory>
#include <set>
#include <variant>

namespace ssc::design {
namespace {
namespace rt = ssc_study;
using J = json::Value;
J obj(std::initializer_list<std::pair<std::string, J>> fields) {
  J j = J::object();
  for (const auto &f : fields)
    j.set(f.first, f.second);
  return j;
}
enum class Type { real, boolean, array, lambda };
const char *typeName(Type t) {
  return t == Type::real      ? "Real"
         : t == Type::boolean ? "Bool"
         : t == Type::array   ? "Array Real"
                              : "function";
}
struct Expr;
using E = std::shared_ptr<Expr>;
struct Expr {
  std::string op;
  std::string literal;
  E exactSource;
  Type type = Type::real;
  double value = 0;
  int id = 0, serial = 0;
  std::vector<E> a;
  bool random = false;
};
struct TokenText {
  std::string text;
  std::size_t at = 0;
  bool number = false;
};
struct Function {
  std::vector<int> ids;
  std::vector<Type> types;
  std::vector<bool> naturals;
  E body;
  bool naturalResult = false;
};
using V = std::variant<double, bool, rt::Vec>;
using Env = std::map<int, V>;
double real(const V &v) {
  if (!std::holds_alternative<double>(v))
    throw std::runtime_error("Expected a real value.");
  return std::get<double>(v);
}
bool boolean(const V &v) {
  if (!std::holds_alternative<bool>(v))
    throw std::runtime_error("Expected a Boolean value.");
  return std::get<bool>(v);
}
const rt::Vec &array(const V &v) {
  if (!std::holds_alternative<rt::Vec>(v))
    throw std::runtime_error("Expected an array.");
  return std::get<rt::Vec>(v);
}
V eval(const E &, Env &, rt::Context &);
struct Program {
  E expression, event;
  E proofLeft, proofRight, proofLo, proofHi;
  bool probabilityQuantity = false;
  std::vector<std::pair<std::string, E>> analyticStages;
  rt::Query query;
  std::string source, reason = "A deterministic user-defined expression.",
                      condition;
  int nodes = 0;
};
class Compiler {
  std::string source;
  std::size_t pos = 0;
  TokenText tok;
  int slot = 1, serial = 0, depth = 0, cloneDepth = 0;
  std::map<std::string, std::pair<int, Type>> locals;
  std::map<std::string, E> constants;
  std::map<std::string, Function> functions;
  std::string requestedMethod = "auto";
  std::vector<std::pair<std::string, E>> analyticStages;
  const std::map<std::string, int> builtins{{"sqrt", 1},
                                            {"log", 1},
                                            {"exp", 1},
                                            {"abs", 1},
                                            {"sin", 1},
                                            {"cos", 1},
                                            {"floor", 1},
                                            {"ceil", 1},
                                            {"min", 2},
                                            {"max", 2},
                                            {"pow", 2},
                                            {"logistic", 1},
                                            {"normalCDF", 1},
                                            {"normalSF", 1},
                                            {"normalQuantile", 1},
                                            {"normal", 2},
                                            {"uniform", 2},
                                            {"bernoulli", 1},
                                            {"exponential", 1},
                                            {"gamma", 2},
                                            {"binomial", 2},
                                            {"normalSample", 3},
                                            {"lognormalSample", 3},
                                            {"bernoulliSample", 2},
                                            {"exponentialSample", 2},
                                            {"sum", 1},
                                            {"mean", 1},
                                            {"variance", 1},
                                            {"sd", 1},
                                            {"length", 1},
                                            {"dot", 2},
                                            {"quantile", 2},
                                            {"probability", 1},
                                            {"sample", 2},
                                            {"map", 2},
                                            {"integrate", 4}};
  [[noreturn]] void error(const std::string &message) const {
    int line = 1, column = 1;
    for (std::size_t i = 0; i < tok.at && i < source.size(); ++i) {
      unsigned char c = source[i];
      if (c == '\n') {
        ++line;
        column = 1;
      } else if ((c & 0xc0) != 0x80)
        ++column;
    }
    throw std::invalid_argument("Study source, line " + std::to_string(line) +
                                ", column " + std::to_string(column) + ": " +
                                message);
  }
  void next() {
    while (pos < source.size()) {
      if (std::isspace(static_cast<unsigned char>(source[pos]))) {
        ++pos;
        continue;
      }
      if (source.compare(pos, 2, "--") == 0) {
        while (pos < source.size() && source[pos] != '\n')
          ++pos;
        continue;
      }
      if (source.compare(pos, 2, "/-") == 0) {
        pos += 2;
        int d = 1;
        while (pos < source.size() && d) {
          if (source.compare(pos, 2, "/-") == 0) {
            ++d;
            pos += 2;
          } else if (source.compare(pos, 2, "-/") == 0) {
            --d;
            pos += 2;
          } else
            ++pos;
        }
        if (d)
          error("Unterminated comment.");
        continue;
      }
      break;
    }
    tok = {"", pos, false};
    if (pos == source.size())
      return;
    for (const auto &op :
         std::vector<std::pair<std::string, std::string>>{{"≥", ">="},
                                                          {"≤", "<="},
                                                          {"≠", "!="},
                                                          {"∧", "&&"},
                                                          {"∨", "||"},
                                                          {"¬", "!"},
                                                          {"ℝ", "Real"},
                                                          {"ℕ", "Nat"},
                                                          {"∈", "in"},
                                                          {":=", ":="},
                                                          {"=>", "=>"},
                                                          {">=", ">="},
                                                          {"<=", "<="},
                                                          {"==", "="},
                                                          {"!=", "!="},
                                                          {"&&", "&&"},
                                                          {"||", "||"}})
      if (source.compare(pos, op.first.size(), op.first) == 0) {
        tok.text = op.second;
        pos += op.first.size();
        return;
      }
    unsigned char c = source[pos];
    std::size_t start = pos++;
    if (std::isdigit(c) ||
        (c == '.' && pos < source.size() &&
         std::isdigit(static_cast<unsigned char>(source[pos])))) {
      while (pos < source.size() &&
             (std::isdigit(static_cast<unsigned char>(source[pos])) ||
              source[pos] == '.'))
        ++pos;
      if (pos < source.size() && (source[pos] == 'e' || source[pos] == 'E')) {
        ++pos;
        if (pos < source.size() && (source[pos] == '+' || source[pos] == '-'))
          ++pos;
        while (pos < source.size() &&
               std::isdigit(static_cast<unsigned char>(source[pos])))
          ++pos;
      }
      tok.text = source.substr(start, pos - start);
      tok.number = true;
      return;
    }
    if (std::isalpha(c) || c == '_' || c >= 128) {
      while (pos < source.size()) {
        unsigned char z = source[pos];
        bool operatorStart = false;
        for (const auto *symbol : {"≥", "≤", "≠", "∧", "∨", "¬", "∈"})
          operatorStart |=
              source.compare(pos, std::char_traits<char>::length(symbol),
                             symbol) == 0;
        if (operatorStart)
          break;
        if (!(std::isalnum(z) || z == '_' || z == '.' || z >= 128))
          break;
        ++pos;
      }
      tok.text = source.substr(start, pos - start);
      return;
    }
    if (std::string("()[]{}:,+-*/^<>=!;").find(char(c)) != std::string::npos) {
      tok.text = char(c);
      return;
    }
    error("Unsupported token. This is the Lean-style study language, not "
          "unrestricted Lean.");
  }
  bool take(const std::string &s) {
    if (tok.text == s) {
      next();
      return true;
    }
    return false;
  }
  void expect(const std::string &s) {
    if (!take(s))
      error("Expected '" + s + "'.");
  }
  std::string name() {
    if (tok.text.empty() || tok.number ||
        !(std::isalpha(static_cast<unsigned char>(tok.text[0])) ||
          tok.text[0] == '_' || static_cast<unsigned char>(tok.text[0]) >= 128))
      error("Expected a name.");
    std::string n = tok.text;
    next();
    return n;
  }
  E node(const std::string &op, Type t, std::vector<E> a = {}) {
    if (++serial > 30000)
      error("Expanded program exceeds 30000 expression nodes.");
    auto e = std::make_shared<Expr>();
    e->op = op;
    e->type = t;
    e->a = std::move(a);
    e->serial = serial;
    for (auto x : e->a)
      e->random |= x->random;
    return e;
  }
  E num(double v) {
    auto e = node("number", Type::real);
    e->value = rt::finite(v);
    return e;
  }
  E truth(bool b) {
    auto e = node("boolean", Type::boolean);
    e->value = b;
    return e;
  }
  E var(int id, Type t) {
    auto e = node("variable", t);
    e->id = id;
    return e;
  }
  void need(const E &e, Type t) {
    if (e->type != t)
      error("Expected " + std::string(typeName(t)) + ", got " +
            typeName(e->type) + ".");
  }
  std::pair<Type, bool> type() {
    bool nat = tok.text == "Nat";
    Type t;
    if (take("Real") || take("Float") || take("Nat"))
      t = Type::real;
    else if (take("Bool"))
      t = Type::boolean;
    else if (take("Array")) {
      expect("Real");
      t = Type::array;
    } else {
      error("Use Real, Nat, Bool, or Array Real.");
    }
    return {t, nat};
  }
  E clone(const E &e, std::map<int, int> &ids) {
    if (++cloneDepth > 120)
      error("Expanded expression nesting exceeds 120 levels.");
    auto r = node(e->op, e->type);
    r->value = e->value;
    r->literal = e->literal;
    r->exactSource = e->exactSource;
    r->random = e->random;
    if (e->op == "let" || e->op == "lambda")
      ids[e->id] = slot++;
    r->id = ids.count(e->id) ? ids[e->id] : e->id;
    for (auto a : e->a)
      r->a.push_back(clone(a, ids));
    --cloneDepth;
    return r;
  }
  E call(const std::string &name, std::vector<E> a) {
    if (functions.count(name)) {
      const auto &f = functions.at(name);
      if (a.size() != f.ids.size())
        error("Wrong argument count for " + name + ".");
      std::map<int, int> ids;
      for (std::size_t i = 0; i < a.size(); ++i) {
        need(a[i], f.types[i]);
        ids[f.ids[i]] = slot++;
        if (f.naturals[i])
          a[i] = node("natural", Type::real, {a[i]});
      }
      E body = clone(f.body, ids);
      for (std::size_t i = a.size(); i-- > 0;) {
        auto e = node("let", body->type, {a[i], body});
        e->id = ids[f.ids[i]];
        body = e;
      }
      return body;
    }
    if (!builtins.count(name) || a.size() != std::size_t(builtins.at(name)))
      error("Unknown function or wrong argument count: " + name + ".");
    Type t = Type::real;
    if (name == "probability") {
      need(a[0], Type::boolean);
      t = Type::real;
    } else if (name == "sample" || name == "map") {
      need(a[0], name == "sample" ? Type::real : Type::array);
      need(a[1], Type::lambda);
      t = Type::array;
    } else if (name == "integrate") {
      need(a[0], Type::lambda);
      for (std::size_t i = 1; i < a.size(); ++i)
        need(a[i], Type::real);
      if (a[0]->random)
        error("integrate requires a deterministic integrand.");
    } else if (name == "sum" || name == "mean" || name == "variance" ||
               name == "sd" || name == "length" || name == "quantile") {
      need(a[0], Type::array);
      if (a.size() > 1)
        need(a[1], Type::real);
    } else if (name == "dot") {
      need(a[0], Type::array);
      need(a[1], Type::array);
    } else {
      for (auto x : a)
        need(x, Type::real);
      if (name.find("Sample") != std::string::npos)
        t = Type::array;
    }
    auto e = node(name, t, std::move(a));
    if (name == "normal" || name == "uniform" || name == "bernoulli" ||
        name == "exponential" || name == "gamma" || name == "binomial" ||
        name.find("Sample") != std::string::npos)
      e->random = true;
    if (name == "probability")
      e->random = false;
    return e;
  }
  int precedence(const std::string &op) {
    if (op == "||")
      return 1;
    if (op == "&&")
      return 2;
    if (op == "=" || op == "!=" || op == "<" || op == ">" || op == "<=" ||
        op == ">=")
      return 3;
    if (op == "+" || op == "-")
      return 4;
    if (op == "*" || op == "/")
      return 5;
    if (op == "^")
      return 7;
    return -1;
  }
  E binary(const std::string &op, E a, E b) {
    bool logic = op == "&&" || op == "||", compare = precedence(op) == 3;
    need(a, logic ? Type::boolean : Type::real);
    need(b, logic ? Type::boolean : Type::real);
    return node(op, logic || compare ? Type::boolean : Type::real, {a, b});
  }
  E expr(int minimum = 0) {
    if (++depth > 120)
      error("Expression nesting exceeds 120 levels.");
    E a = prefix();
    while (true) {
      if (take("[")) {
        need(a, Type::array);
        E i = expr();
        need(i, Type::real);
        expect("]");
        a = node("at", Type::real, {a, i});
        continue;
      }
      int p = precedence(tok.text);
      if (p < minimum)
        break;
      std::string op = tok.text;
      next();
      E b = expr(p + (op == "^" ? 0 : 1));
      a = binary(op, a, b);
    }
    --depth;
    return a;
  }
  E prefix() {
    if (take("-")) {
      auto a = expr(6);
      need(a, Type::real);
      return node("neg", Type::real, {a});
    }
    if (take("+"))
      return expr(6);
    if (take("!")) {
      auto a = expr(6);
      need(a, Type::boolean);
      return node("not", Type::boolean, {a});
    }
    if (take("(")) {
      auto e = expr();
      expect(")");
      return e;
    }
    if (take("[")) {
      std::vector<E> a;
      if (!take("]")) {
        do {
          auto x = expr();
          need(x, Type::real);
          a.push_back(x);
        } while (take(","));
        expect("]");
      }
      return node("array", Type::array, a);
    }
    if (take("let")) {
      auto n = name();
      std::pair<Type, bool> annotation{Type::real, false};
      bool annotated = take(":");
      if (annotated)
        annotation = type();
      expect(":=");
      E value = expr();
      if (annotated) {
        need(value, annotation.first);
        if (annotation.second)
          value = node("natural", Type::real, {value});
      }
      take(";");
      take("in");
      auto before = locals;
      int id = slot++;
      locals[n] = {id, value->type};
      E body = expr();
      locals = before;
      auto e = node("let", body->type, {value, body});
      e->id = id;
      return e;
    }
    if (take("if")) {
      auto c = expr();
      need(c, Type::boolean);
      expect("then");
      auto a = expr();
      expect("else");
      auto b = expr();
      need(b, a->type);
      return node("if", a->type, {c, a, b});
    }
    if (take("fun")) {
      auto n = name();
      if (take(":")) {
        auto t = type();
        if (t.first != Type::real)
          error("Lambda arguments must be numeric.");
      }
      expect("=>");
      auto before = locals;
      int id = slot++;
      locals[n] = {id, Type::real};
      auto body = expr();
      locals = before;
      need(body, Type::real);
      auto e = node("lambda", Type::lambda, {body});
      e->id = id;
      return e;
    }
    if (tok.number) {
      double d;
      if (!text::parse_double(tok.text, d) || !std::isfinite(d))
        error("Invalid finite numeric literal.");
      auto e = num(d);
      e->literal = tok.text;
      next();
      return e;
    }
    if (take("true"))
      return truth(true);
    if (take("false"))
      return truth(false);
    auto n = name();
    if (n == "pi")
      return num(3.1415926535897932384626433832795);
    if (n.rfind("Real.", 0) == 0)
      n = n.substr(5);
    if (locals.count(n))
      return var(locals[n].first, locals[n].second);
    if (constants.count(n)) {
      std::map<int, int> ids;
      return clone(constants[n], ids);
    }
    if (!functions.count(n) && !builtins.count(n))
      error("Unknown name '" + n + "'. Define functions before using them.");
    int arity =
        functions.count(n) ? int(functions[n].ids.size()) : builtins.at(n);
    std::vector<E> a;
    if (take("(")) {
      if (!take(")")) {
        do {
          a.push_back(expr());
        } while (take(","));
        expect(")");
      }
      // Accept Lean's curried unary-parenthesized syntax as well as f(a,b).
      while (int(a.size()) < arity)
        a.push_back(expr(8));
    } else
      for (int i = 0; i < arity; ++i)
        a.push_back(expr(8));
    return call(n, a);
  }
  double constant(E e) {
    if (e->random)
      error("This setting must be deterministic.");
    Env env;
    rt::Context c;
    try {
      return real(eval(e, env, c));
    } catch (const std::exception &ex) {
      error(std::string("Expected a constant expression: ") + ex.what());
    }
  }
  bool contains(const E &e, const std::string &op) {
    if (e->op == op)
      return true;
    for (auto a : e->a)
      if (contains(a, op))
        return true;
    return false;
  }
  E inlineLets(E e, std::map<int, E> env = {}, int level = 0) {
    if (level > 120)
      error("Expanded expression nesting exceeds 120 levels.");
    if (e->op == "variable" && env.count(e->id))
      return env[e->id];
    if (e->op == "let") {
      env[e->id] = inlineLets(e->a[0], env, level + 1);
      return inlineLets(e->a[1], env, level + 1);
    }
    auto r = std::make_shared<Expr>(*e);
    r->a.clear();
    for (auto a : e->a)
      r->a.push_back(inlineLets(a, env, level + 1));
    r->random = e->op == "normal" || e->op == "uniform" ||
                e->op == "bernoulli" || e->op == "exponential" ||
                e->op == "gamma" || e->op == "binomial" ||
                e->op.find("Sample") != std::string::npos;
    if (e->op != "probability")
      for (auto a : r->a)
        r->random |= a->random;
    return r;
  }
  E probabilityEvent(E e) {
    if (e->op == "probability") {
      if (contains(e->a[0], "probability"))
        error("Nested probability expressions are not supported.");
      return e->a[0];
    }
    if (e->op == "let") {
      auto r = std::make_shared<Expr>(*e);
      r->type = Type::boolean;
      r->a[1] = probabilityEvent(e->a[1]);
      return r;
    }
    error(
        "Simulation fallback needs probability(event) as the query quantity.");
  }
  void sizeCheck(const E &e, int &count, int level = 0) {
    if (++count > 30000 || level > 120)
      error("Expanded expression is too large or deeply nested.");
    for (auto a : e->a)
      sizeCheck(a, count, level + 1);
  }
  struct Gaussian {
    E mean;
    std::map<int, E> noise;
  };
  struct Unsupported {};
  Gaussian gaussian(E e) {
    if (!e->random)
      return {e, {}};
    if (e->op == "normal") {
      if (e->a[1]->random)
        throw Unsupported{};
      auto g = gaussian(e->a[0]);
      g.noise[e->serial] = node("checkedSD", Type::real, {e->a[1]});
      return g;
    }
    if (e->op == "neg") {
      auto g = gaussian(e->a[0]);
      g.mean = node("neg", Type::real, {g.mean});
      for (auto &c : g.noise)
        c.second = node("neg", Type::real, {c.second});
      return g;
    }
    if (e->op == "+" || e->op == "-") {
      auto g = gaussian(e->a[0]), h = gaussian(e->a[1]);
      g.mean = binary(e->op, g.mean, h.mean);
      for (auto c : h.noise) {
        E value = e->op == "-" ? node("neg", Type::real, {c.second}) : c.second;
        if (g.noise.count(c.first))
          g.noise[c.first] = binary("+", g.noise[c.first], value);
        else
          g.noise[c.first] = value;
      }
      return g;
    }
    if (e->op == "*" || e->op == "/") {
      E x = e->a[0], scale = e->a[1];
      if (e->op == "*" && !x->random)
        std::swap(x, scale);
      if (scale->random)
        throw Unsupported{};
      auto g = gaussian(x);
      g.mean = binary(e->op, g.mean, scale);
      for (auto &c : g.noise)
        c.second = binary(e->op, c.second, scale);
      return g;
    }
    throw Unsupported{};
  }
  E analytic(E event) {
    event = inlineLets(event);
    if (!event->random)
      return node("if", Type::real, {event, num(1), num(0)});
    std::string op = event->op;
    if (op != "<" && op != "<=" && op != ">" && op != ">=")
      throw Unsupported{};
    E left = event->a[0], right = event->a[1];
    bool absolute = left->op == "abs" && !right->random;
    Gaussian g =
        absolute ? gaussian(left->a[0]) : gaussian(binary("-", left, right));
    E variance = num(0);
    for (auto c : g.noise)
      variance = binary("+", variance, binary("*", c.second, c.second));
    E sd = node("sqrt", Type::real, {variance});
    if (analyticStages.size() < 15) {
      analyticStages.push_back(
          {"Mean of the affine Gaussian statistic", g.mean});
      analyticStages.push_back(
          {"Variance, including covariance from shared draws", variance});
    }
    bool greater = op == ">" || op == ">=";
    E p, degenerate;
    if (absolute) {
      E low = binary("/", binary("-", node("neg", Type::real, {right}), g.mean),
                     sd),
        high = binary("/", binary("-", right, g.mean), sd);
      E inside = binary("-", node("normalCDF", Type::real, {high}),
                        node("normalCDF", Type::real, {low}));
      E outside = binary("+", node("normalCDF", Type::real, {low}),
                         node("normalSF", Type::real, {high}));
      p = node("if", Type::real,
               {binary("<", right, num(0)), num(greater ? 1 : 0),
                greater ? outside : inside});
      degenerate = binary(op, node("abs", Type::real, {g.mean}), right);
    } else {
      E z = binary("/", g.mean, sd);
      p = node(greater ? "normalCDF" : "normalSF", Type::real, {z});
      degenerate = binary(op, g.mean, num(0));
    }
    return node("if", Type::real,
                {binary("=", variance, num(0)),
                 node("if", Type::real, {degenerate, num(1), num(0)}), p});
  }
  E lowerProbabilities(E e, bool &lowered) {
    if (e->op == "probability") {
      if (contains(e->a[0], "probability"))
        error("Nested probability expressions are not supported.");
      auto r = analytic(e->a[0]);
      lowered = true;
      return r;
    }
    auto r = std::make_shared<Expr>(*e);
    r->a.clear();
    for (auto a : e->a)
      r->a.push_back(lowerProbabilities(a, lowered));
    r->random = e->op == "normal" || e->op == "uniform" ||
                e->op == "bernoulli" || e->op == "exponential" ||
                e->op == "gamma" || e->op == "binomial" ||
                e->op.find("Sample") != std::string::npos;
    for (auto a : r->a)
      r->random |= a->random;
    return r;
  }
  std::vector<double> polynomial(E e) {
    e = inlineLets(e);
    if (e->op == "number")
      return {e->value};
    if (e->op == "variable" && e->id == 0)
      return {0, 1};
    if (e->op == "neg") {
      auto a = polynomial(e->a[0]);
      for (double &x : a)
        x = -x;
      return a;
    }
    if (e->op == "+" || e->op == "-") {
      auto a = polynomial(e->a[0]), b = polynomial(e->a[1]);
      if (a.empty() || b.empty())
        return {};
      a.resize(std::max(a.size(), b.size()));
      for (std::size_t i = 0; i < b.size(); ++i)
        a[i] += e->op == "+" ? b[i] : -b[i];
      return a;
    }
    if (e->op == "*" || e->op == "/") {
      auto a = polynomial(e->a[0]), b = polynomial(e->a[1]);
      if (a.empty() || b.empty())
        return {};
      if (e->op == "/") {
        if (b.size() != 1 || b[0] == 0)
          return {};
        for (double &x : a)
          x /= b[0];
        return a;
      }
      if (a.size() + b.size() > 4)
        return {};
      std::vector<double> c(a.size() + b.size() - 1, 0);
      for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j)
          c[i + j] += a[i] * b[j];
      return c;
    }
    if (e->op == "^" && e->a[1]->op == "number") {
      double n = e->a[1]->value;
      if (n == 0)
        return {1};
      if (n == 1)
        return polynomial(e->a[0]);
      if (n == 2)
        return polynomial(binary("*", e->a[0], e->a[0]));
    }
    if (!contains(e, "variable") && !contains(e, "probability") && !e->random) {
      try {
        return {constant(e)};
      } catch (const std::exception &) {
        return {};
      }
    }
    return {};
  }

public:
  explicit Compiler(std::string s) : source(std::move(s)) {
    if (source.empty() || source.size() > 131072)
      throw std::invalid_argument("Use 1–131072 bytes of study source.");
    next();
  }
  Program compile() {
    Program p;
    p.source = source;
    if (take("import")) {
      auto m = name();
      if (m != "DesignCompiler.Input" && m != "DesignCompiler.Study")
        error("The study runtime is built in; external imports are not "
              "supported.");
    }
    if (take("open"))
      expect("DesignCompiler");
    while (take("def")) {
      auto n = name();
      if (constants.count(n) || functions.count(n) || builtins.count(n) ||
          n == "pi")
        error("Duplicate or reserved name '" + n + "'.");
      Function f;
      auto before = locals;
      while (take("(")) {
        auto param = name();
        expect(":");
        auto t = type();
        expect(")");
        if (locals.count(param))
          error("Duplicate parameter.");
        int id = slot++;
        locals[param] = {id, t.first};
        f.ids.push_back(id);
        f.types.push_back(t.first);
        f.naturals.push_back(t.second);
      }
      bool annotated = take(":");
      std::pair<Type, bool> t{Type::real, false};
      if (annotated)
        t = type();
      expect(":=");
      f.body = expr();
      locals = before;
      if (annotated) {
        need(f.body, t.first);
        if (t.second)
          f.body = node("natural", Type::real, {f.body});
      }
      if (f.ids.empty()) {
        if (f.body->random)
          error("Put random draws inside a function or a local let binding.");
        constants[n] =
            f.body->type == Type::real && !contains(f.body, "probability")
                ? num(constant(f.body))
                : f.body;
        if (constants[n] != f.body)
          constants[n]->exactSource = f.body;
      } else
        functions[n] = f;
      take(";");
    }
    if (take("simulation")) {
      expect("{");
      std::set<std::string> seen;
      while (!take("}")) {
        auto k = name();
        if (!seen.insert(k).second)
          error("Duplicate simulation setting.");
        expect(":=");
        if (k == "method") {
          requestedMethod = name();
          if (requestedMethod != "auto" && requestedMethod != "analytic" &&
              requestedMethod != "simulation")
            error("method must be auto, analytic, or simulation.");
        } else {
          double x = constant(expr());
          if (k == "replicates")
            p.query.replicates = rt::count(x, 100, 1000000);
          else if (k == "seed") {
            if (x < 0 || x > 9007199254740991.0 || x != std::floor(x))
              error("seed must be a nonnegative exact integer.");
            p.query.seed = std::uint64_t(x);
          } else if (k == "confidence") {
            if (!(x > 0.5 && x < 1))
              error("confidence must be between 0.5 and 1.");
            p.query.confidence = x;
          } else
            error("Unknown simulation setting '" + k + "'.");
        }
        take(",");
      }
    }
    if (take("evaluate")) {
      p.query.evaluate = true;
      p.expression = expr();
      need(p.expression, Type::real);
    } else {
      expect("solve");
      expect("for");
      p.query.variable = name();
      if (constants.count(p.query.variable) ||
          functions.count(p.query.variable) ||
          builtins.count(p.query.variable) || p.query.variable == "pi")
        error("The solve variable must be a fresh name.");
      expect(":");
      auto domain = name();
      if (domain != "Nat" && domain != "Real" && domain != "Grid")
        error("Solve over Nat, Real, or Grid.");
      p.query.integer = domain != "Real";
      expect("in");
      expect("[");
      p.proofLo = expr();
      p.query.lo = constant(p.proofLo);
      expect(",");
      p.proofHi = expr();
      p.query.hi = constant(p.proofHi);
      expect("]");
      if (take("step"))
        p.query.step = constant(expr());
      if (take("tolerance"))
        p.query.tolerance = constant(expr());
      if (!(p.query.lo <= p.query.hi && p.query.step > 0 &&
            p.query.tolerance > 0 && p.query.tolerance <= 0.1))
        error("Invalid search bounds, step, or tolerance.");
      if (domain == "Nat") {
        rt::count(p.query.lo, 0, 10000000);
        rt::count(p.query.hi, 0, 10000000);
        rt::count(p.query.step, 1, 10000000);
      }
      locals[p.query.variable] = {0, Type::real};
      expect("where");
      auto start = tok.at;
      auto condition = expr();
      p.condition = source.substr(start);
      need(condition, Type::boolean);
      if (precedence(condition->op) != 3 || condition->op == "!=")
        error("Use a comparison =, <, <=, >, or >= as the solve condition.");
      p.query.comparison = condition->op;
      p.proofLeft = condition->a[0];
      p.proofRight = condition->a[1];
      p.expression = condition->a[0];
      try {
        p.query.target = constant(condition->a[1]);
      } catch (const std::exception &) {
        p.expression = binary("-", condition->a[0], condition->a[1]);
        p.query.target = 0;
      }
      if (!p.query.integer && p.query.comparison != "=")
        error("Real solving finds an equation root. Use =, or use Grid with an "
              "explicit step for an inequality.");
    }
    take(";");
    if (!tok.text.empty())
      error("Unexpected code after the query. Use one solve or evaluate "
            "command.");
    auto flat = inlineLets(p.expression);
    p.probabilityQuantity = flat->op == "probability";
    int expanded = 0;
    sizeCheck(flat, expanded);
    bool hasProbability = contains(p.expression, "probability");
    bool lowered = false;
    if (hasProbability) {
      try {
        if (requestedMethod == "simulation")
          throw Unsupported{};
        p.expression = lowerProbabilities(flat, lowered);
        p.query.evaluation = "analytic_gaussian_probability";
        p.reason = "The event reduces to a comparison of an affine Gaussian "
                   "expression. Shared draws and independence are tracked "
                   "symbolically.";
      } catch (Unsupported &) {
        if (requestedMethod == "analytic")
          error("No analytic probability lowering is available for this event. "
                "Use auto or simulation.");
        if (flat->op != "probability")
          error("Simulation fallback requires probability(event) as the query "
                "quantity, possibly through a user function.");
        p.event = probabilityEvent(p.expression);
        p.expression = nullptr;
        p.query.monteCarlo = true;
        p.query.evaluation = "monte_carlo";
        p.reason =
            "The probability event is outside the analytic Gaussian rules; "
            "execute the supplied generative model and count the event.";
        if (!p.query.evaluate && p.query.comparison == "=")
          error("Simulation cannot certify exact equality. Use >= or <= on a "
                "finite grid.");
      }
    }
    if (p.expression && p.expression->random)
      error("Wrap a random Boolean event in probability(...). Raw random "
            "values are not deterministic solve objectives.");
    if (p.expression && !p.query.evaluate && !p.query.integer) {
      p.query.polynomial = polynomial(p.expression);
      for (double x : p.query.polynomial)
        if (!std::isfinite(x)) {
          p.query.polynomial.clear();
          break;
        }
    }
    int finalNodes = 0;
    sizeCheck(p.expression ? p.expression : p.event, finalNodes);
    p.nodes = finalNodes;
    p.analyticStages = std::move(analyticStages);
    return p;
  }
};
V eval(const E &e, Env &env, rt::Context &c) {
  c.tick();
  const auto &op = e->op;
  if (op == "number")
    return e->value;
  if (op == "boolean")
    return bool(e->value);
  if (op == "variable") {
    auto i = env.find(e->id);
    if (i == env.end())
      throw std::runtime_error("Expression depends on an unbound variable.");
    return i->second;
  }
  if (op == "let") {
    V value = eval(e->a[0], env, c);
    auto old = env.find(e->id);
    bool present = old != env.end();
    V saved = present ? old->second : V(0.0);
    env[e->id] = std::move(value);
    V result = eval(e->a[1], env, c);
    if (present)
      env[e->id] = std::move(saved);
    else
      env.erase(e->id);
    return result;
  }
  if (op == "if")
    return eval(e->a[boolean(eval(e->a[0], env, c)) ? 1 : 2], env, c);
  if (op == "&&") {
    bool a = boolean(eval(e->a[0], env, c));
    return a && boolean(eval(e->a[1], env, c));
  }
  if (op == "||") {
    bool a = boolean(eval(e->a[0], env, c));
    return a || boolean(eval(e->a[1], env, c));
  }
  if (op == "sample" || op == "map" || op == "integrate") {
    std::size_t index = op == "integrate" ? 0 : 1;
    E lambda = e->a[index];
    auto apply = [&](double x) {
      auto old = env.find(lambda->id);
      bool present = old != env.end();
      V saved = present ? old->second : V(0.0);
      env[lambda->id] = x;
      double value = real(eval(lambda->a[0], env, c));
      if (present)
        env[lambda->id] = std::move(saved);
      else
        env.erase(lambda->id);
      return value;
    };
    if (op == "sample")
      return rt::sample(c, real(eval(e->a[0], env, c)), apply);
    if (op == "map") {
      auto xs = array(eval(e->a[0], env, c));
      rt::Vec out;
      out.reserve(xs.size());
      for (double x : xs) {
        c.tick();
        out.push_back(rt::finite(apply(x)));
      }
      return out;
    }
    double lo = real(eval(e->a[1], env, c)), hi = real(eval(e->a[2], env, c)),
           n = real(eval(e->a[3], env, c));
    return rt::integrate(c, apply, lo, hi, n);
  }
  std::vector<V> a;
  for (auto x : e->a)
    a.push_back(eval(x, env, c));
  if (op == "array") {
    rt::Vec out;
    for (const auto &x : a)
      out.push_back(real(x));
    return out;
  }
  if (op == "not")
    return !boolean(a[0]);
  if (op == "sum")
    return rt::sum(array(a[0]));
  if (op == "mean")
    return rt::mean(array(a[0]));
  if (op == "variance")
    return rt::variance(array(a[0]));
  if (op == "sd")
    return std::sqrt(rt::variance(array(a[0])));
  if (op == "length")
    return double(array(a[0]).size());
  if (op == "dot")
    return rt::dot(array(a[0]), array(a[1]));
  if (op == "quantile")
    return rt::quantile(array(a[0]), real(a[1]));
  if (op == "at")
    return rt::at(array(a[0]), real(a[1]));
  double x = a.empty() ? 0 : real(a[0]), y = a.size() > 1 ? real(a[1]) : 0,
         z = a.size() > 2 ? real(a[2]) : 0;
  if (op == "+")
    return rt::finite(x + y);
  if (op == "-")
    return rt::finite(x - y);
  if (op == "*")
    return rt::finite(x * y);
  if (op == "/")
    return rt::divide(x, y);
  if (op == "^" || op == "pow")
    return rt::power(x, y);
  if (op == "=")
    return x == y;
  if (op == "!=")
    return x != y;
  if (op == "<")
    return x < y;
  if (op == "<=")
    return x <= y;
  if (op == ">")
    return x > y;
  if (op == ">=")
    return x >= y;
  if (op == "neg")
    return -x;
  if (op == "natural")
    return rt::natural(x);
  if (op == "checkedSD")
    return rt::checkedSD(x);
  if (op == "sqrt")
    return rt::squareRoot(x);
  if (op == "log")
    return rt::logarithm(x);
  if (op == "exp")
    return rt::finite(std::exp(x));
  if (op == "abs")
    return std::abs(x);
  if (op == "sin")
    return std::sin(x);
  if (op == "cos")
    return std::cos(x);
  if (op == "floor")
    return std::floor(x);
  if (op == "ceil")
    return std::ceil(x);
  if (op == "min")
    return std::min(x, y);
  if (op == "max")
    return std::max(x, y);
  if (op == "logistic")
    return rt::logistic(x);
  if (op == "normalCDF")
    return rt::normalCDF(x);
  if (op == "normalSF")
    return rt::normalSF(x);
  if (op == "normalQuantile")
    return rt::normalQuantile(x);
  if (op == "normal")
    return rt::normal(c, x, y);
  if (op == "uniform")
    return rt::uniform(c, x, y);
  if (op == "bernoulli")
    return rt::bernoulli(c, x);
  if (op == "exponential")
    return rt::exponential(c, x);
  if (op == "gamma")
    return rt::gamma(c, x, y);
  if (op == "binomial")
    return rt::binomial(c, x, y);
  if (op == "normalSample")
    return rt::normalSample(c, x, y, z);
  if (op == "lognormalSample")
    return rt::lognormalSample(c, x, y, z);
  if (op == "bernoulliSample")
    return rt::bernoulliSample(c, x, y);
  if (op == "exponentialSample")
    return rt::exponentialSample(c, x, y);
  throw std::runtime_error("Unsupported compiled operation: " + op);
}
J wire(const E &e) {
  if (!e)
    return J();
  J a = J::array();
  for (auto x : e->a)
    a.push(wire(x));
  return obj({{"op", e->op},
              {"type", typeName(e->type)},
              {"value", e->value},
              {"id", e->id},
              {"args", a}});
}
J queryWire(const rt::Query &q) {
  J polynomial = J::array();
  for (double x : q.polynomial)
    polynomial.push(x);
  return obj({{"evaluate", q.evaluate},
              {"integer", q.integer},
              {"monteCarlo", q.monteCarlo},
              {"variable", q.variable},
              {"comparison", q.comparison},
              {"evaluation", q.evaluation},
              {"lo", q.lo},
              {"hi", q.hi},
              {"step", q.step},
              {"target", q.target},
              {"tolerance", q.tolerance},
              {"replicates", q.replicates},
              {"seed", static_cast<unsigned long long>(q.seed)},
              {"confidence", q.confidence},
              {"polynomial", polynomial}});
}
J ir(const Program &p) {
  return obj({{"expression", wire(p.expression)},
              {"event", wire(p.event)},
              {"query", queryWire(p.query)}});
}
std::string ctype(Type t) {
  return t == Type::real ? "double" : t == Type::boolean ? "bool" : "rt::Vec";
}
std::string emit(const E &e) {
  auto v = [](int id) { return "v" + std::to_string(id); };
  if (e->op == "number") {
    auto value = rt::number(e->value);
    // Source numbers have Real semantics even when written without a dot.
    // Integer C++ temporaries would overflow before rt::finite sees a product.
    if (value.find_first_of(".eE") == std::string::npos)
      value += ".0";
    return value;
  }
  if (e->op == "boolean")
    return e->value ? "true" : "false";
  if (e->op == "variable")
    return v(e->id);
  if (e->op == "let")
    return "([&]() -> " + ctype(e->type) + " { const auto " + v(e->id) + " = " +
           emit(e->a[0]) + "; return " + emit(e->a[1]) + "; }())";
  if (e->op == "lambda")
    return "[&](double " + v(e->id) + ") -> double { return " + emit(e->a[0]) +
           "; }";
  if (e->op == "if")
    return "(" + emit(e->a[0]) + " ? " + emit(e->a[1]) + " : " + emit(e->a[2]) +
           ")";
  if (e->op == "&&" || e->op == "||")
    return "(" + emit(e->a[0]) + " " + e->op + " " + emit(e->a[1]) + ")";
  std::string body = "([&]() -> " + ctype(e->type) + " { ";
  std::vector<std::string> a;
  for (std::size_t i = 0; i < e->a.size(); ++i) {
    auto n = "a" + std::to_string(e->serial) + "_" + std::to_string(i);
    body += "const auto " + n + " = " + emit(e->a[i]) + "; ";
    a.push_back(n);
  }
  std::string x = a.empty() ? "" : a[0], y = a.size() > 1 ? a[1] : "",
              z = a.size() > 2 ? a[2] : "", s, op = e->op;
  if (op == "+" || op == "-" || op == "*")
    s = "rt::finite(" + x + op + y + ")";
  else if (op == "=" || op == "!=" || op == "<" || op == ">" || op == "<=" ||
           op == ">=")
    s = "(" + x + (op == "=" ? "==" : op) + y + ")";
  else if (op == "neg")
    s = "(-" + x + ")";
  else if (op == "not")
    s = "(!" + x + ")";
  else if (op == "/")
    s = "rt::divide(" + x + "," + y + ")";
  else if (op == "^" || op == "pow")
    s = "rt::power(" + x + "," + y + ")";
  else if (op == "sqrt")
    s = "rt::squareRoot(" + x + ")";
  else if (op == "log")
    s = "rt::logarithm(" + x + ")";
  else if (op == "exp")
    s = "rt::finite(std::exp(" + x + "))";
  else if (op == "abs" || op == "sin" || op == "cos" || op == "floor" ||
           op == "ceil")
    s = "std::" + op + "(" + x + ")";
  else if (op == "min" || op == "max")
    s = "std::" + op + "(" + x + "," + y + ")";
  else if (op == "array") {
    s = "rt::Vec{";
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (i)
        s += ",";
      s += a[i];
    }
    s += "}";
  } else if (op == "length")
    s = "double(" + x + ".size())";
  else if (op == "sd")
    s = "std::sqrt(rt::variance(" + x + "))";
  else {
    bool context = op == "normal" || op == "uniform" || op == "bernoulli" ||
                   op == "exponential" || op == "gamma" || op == "binomial" ||
                   op.find("Sample") != std::string::npos || op == "sample" ||
                   op == "map" || op == "integrate";
    s = "rt::" + op + "(" + (context ? "ctx" : "");
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (i || context)
        s += ",";
      s += a[i];
    }
    s += ")";
  }
  return body + "return " + s + "; }())";
}
std::string readable(const E &e, int depth = 0) {
  if (!e)
    return "";
  if (depth > 16)
    return "…";
  // Simplify presentation only. The executable expression keeps its domain
  // checks and evaluation semantics. Limit constant evaluation work here.
  if (!e->random && e->op != "variable" && e->op != "number" &&
      e->op != "boolean") {
    try {
      Env env;
      rt::Context context;
      context.limit = 1000;
      if (e->op == "if")
        return readable(e->a[boolean(eval(e->a[0], env, context)) ? 1 : 2],
                        depth + 1);
      if (e->type == Type::real)
        return rt::number(real(eval(e, env, context)));
    } catch (const std::exception &) {
    }
  }
  if (e->op == "number")
    return rt::number(e->value);
  if (e->op == "boolean")
    return e->value ? "true" : "false";
  if (e->op == "variable")
    return e->id == 0 ? "unknown" : "local" + std::to_string(e->id);
  if (e->op == "if")
    return "(if " + readable(e->a[0], depth + 1) + " then " +
           readable(e->a[1], depth + 1) + " else " +
           readable(e->a[2], depth + 1) + ")";
  if (e->op == "let")
    return "(let local" + std::to_string(e->id) +
           " := " + readable(e->a[0], depth + 1) + "; " +
           readable(e->a[1], depth + 1) + ")";
  if (e->op == "neg")
    return "(-" + readable(e->a[0], depth + 1) + ")";
  if (e->a.size() == 2 &&
      std::string("+ - * / ^ < <= > >= = != && ||").find(e->op) !=
          std::string::npos)
    return "(" + readable(e->a[0], depth + 1) + " " + e->op + " " +
           readable(e->a[1], depth + 1) + ")";
  std::string s = e->op + "(";
  for (std::size_t i = 0; i < e->a.size(); ++i) {
    if (i)
      s += ", ";
    s += readable(e->a[i], depth + 1);
  }
  return s + ")";
}
std::string readableQuery(const E &e, const std::string &variable) {
  auto out = readable(e);
  if (variable.empty())
    return out;
  std::size_t at = 0;
  while ((at = out.find("unknown", at)) != std::string::npos) {
    out.replace(at, 7, variable);
    at += variable.size();
  }
  return out;
}
std::string generate(const Program &p) {
  const auto &q = p.query;
  std::string s = "// Generated from the Lean-style study expression compiler. "
                  "C++17.\n// Build: c++ -std=c++17 -O2 study.cpp -o study\n// "
                  "The expression below is compiled directly; no parser or "
                  "expression interpreter is included.\n// This is not a Lean "
                  "proof certificate. Arithmetic is ordinary binary64.\n";
  s += program_runtime_source;
  s += "\n#include <iostream>\nnamespace rt=ssc_study;\nint main() {\n try {\n "
       " rt::Context ctx; rt::Query q;\n";
  auto set = [&](const std::string &k, const std::string &value) {
    s += "  q." + k + " = " + value + ";\n";
  };
  set("evaluate", q.evaluate ? "true" : "false");
  set("integer", q.integer ? "true" : "false");
  set("monteCarlo", q.monteCarlo ? "true" : "false");
  set("variable", rt::quote(q.variable));
  set("comparison", rt::quote(q.comparison));
  set("evaluation", rt::quote(q.evaluation));
  for (auto kv : std::vector<std::pair<std::string, double>>{
           {"lo", q.lo},
           {"hi", q.hi},
           {"step", q.step},
           {"target", q.target},
           {"tolerance", q.tolerance},
           {"replicates", double(q.replicates)},
           {"seed", double(q.seed)},
           {"confidence", q.confidence}})
    set(kv.first, rt::number(kv.second));
  for (double x : q.polynomial)
    s += "  q.polynomial.push_back(" + rt::number(x) + ");\n";
  s += "  auto expression = [&](double v0, rt::Context& ctx) -> double { "
       "return " +
       (p.expression ? emit(p.expression) : "0.0") + "; };\n";
  s += "  auto event = [&](double v0, rt::Context& ctx) -> bool { return " +
       (p.event ? emit(p.event) : "false") + "; };\n";
  s += "  std::cout << rt::run(q, expression, event, ctx).json() << '\\n';\n  "
       "return 0;\n } catch (const std::exception& error) { std::cerr << "
       "error.what() << '\\n'; return 1; }\n}\n";
  if (s.size() > 2000000)
    throw std::invalid_argument(
        "Generated C++ exceeds 2 MB. Factor or simplify the expressions.");
  return s;
}

algebra::Expression exactExpression(const E &e,
                                    std::map<int, algebra::Expression> env = {},
                                    int depth = 0) {
  using namespace algebra;
  if (!e || depth > 120)
    throw Unsupported{};
  if (e->exactSource)
    return exactExpression(e->exactSource, env, depth + 1);
  if (e->op == "number")
    return constant(decimal(e->literal));
  if (e->op == "variable") {
    if (env.count(e->id))
      return env.at(e->id);
    if (e->id == 0)
      return {Polynomial{Rational(0), Rational(1), Rational(0)}, "x"};
    throw Unsupported{};
  }
  if (e->op == "let") {
    env[e->id] = exactExpression(e->a[0], env, depth + 1);
    return exactExpression(e->a[1], env, depth + 1);
  }
  if (e->op == "neg")
    return binary("-", constant(0), exactExpression(e->a[0], env, depth + 1));
  if (e->op == "+" || e->op == "-" || e->op == "*" || e->op == "/")
    return binary(e->op, exactExpression(e->a[0], env, depth + 1),
                  exactExpression(e->a[1], env, depth + 1));
  if (e->op == "^" || e->op == "pow") {
    auto a = exactExpression(e->a[0], env, depth + 1),
         b = exactExpression(e->a[1], env, depth + 1);
    if (!b.p[1].zero() || !b.p[2].zero() || b.p[0].d != 1 || b.p[0].n < 0 ||
        b.p[0].n > 2)
      throw Unsupported{};
    Expression out = constant(1);
    for (int i = 0; i < b.p[0].n; ++i)
      out = binary("*", out, a);
    return out;
  }
  throw Unsupported{};
}

J derive(const Program &p) {
  J steps = J::array();
  auto step = [&](const std::string &rule, const std::string &expression) {
    steps.push(obj({{"rule", rule}, {"expression", expression}}));
  };
  J result =
      obj({{"kind", "numerical"},
           {"leanKernelChecked", false},
           {"leanSource", ""},
           {"proofScope", "No Lean proof is attached to this calculation."}});
  if (!p.query.evaluate)
    step("Requested condition", p.condition);
  try {
    using namespace algebra;
    if (p.query.evaluate || p.query.integer || p.query.comparison != "=" ||
        p.query.monteCarlo)
      throw Unsupported{};
    auto equation = binary("-", exactExpression(p.proofLeft),
                           exactExpression(p.proofRight));
    auto lo = exactExpression(p.proofLo).p[0],
         hi = exactExpression(p.proofHi).p[0];
    Rational c = equation.p[0], b = equation.p[1], a = equation.p[2];
    std::string poly = "(" + a.text() + ") * x^2 + (" + b.text() + ") * x + (" +
                       c.text() + ")";
    step("Move both sides to the left and collect exact rational coefficients",
         poly + " = 0");
    std::string proof =
        "import Mathlib.Tactic\nimport Mathlib.Data.Real.Sqrt\n\n";
    proof +=
        "-- Exact real-number interpretation of the source arithmetic.\n-- "
        "This proves the algebra, not binary64 execution or a numerical "
        "tolerance.\nnoncomputable section\nnamespace StudyDerivation\n\n";
    proof += "def residual (x : ℝ) : ℝ := " + equation.lean + "\n\n";
    proof += "theorem normalized (x : ℝ) : residual x = " + poly +
             " := by\n  norm_num [residual] <;> ring\n\n";
    std::vector<Rational> rationalRoots;
    bool exactRoots = true, quadratic = !a.zero();
    std::string rootCondition;
    if (!quadratic && !b.zero()) {
      Rational root = (Rational(0) - c) / b;
      rationalRoots.push_back(root);
      rootCondition = "x = " + root.text();
      step("Isolate the unknown", rootCondition);
      proof += "theorem roots (x : ℝ) : residual x = 0 ↔ " + rootCondition +
               " := by\n  rw [normalized]\n  constructor <;> intro h <;> "
               "nlinarith\n\n";
    } else if (quadratic) {
      Rational d = b * b - Rational(4) * a * c;
      std::string square =
          "(2 * (" + a.text() + ") * x + (" + b.text() + "))^2";
      step("Complete the square", square + " = " + d.text());
      step("Compute the discriminant", "D = b² − 4ac = " + d.text());
      proof += "theorem completed_square (x : ℝ) : residual x = 0 ↔ " + square +
               " = " + d.text() +
               " := by\n  rw [normalized]\n  constructor <;> intro h <;> "
               "nlinarith\n\n";
      if (d.n < 0) {
        rootCondition = "False";
        step("A real square cannot be negative", "No real solutions.");
        proof += "theorem roots (x : ℝ) : residual x = 0 ↔ False := by\n  rw "
                 "[completed_square]\n  constructor\n  · intro h\n    "
                 "nlinarith [sq_nonneg (2 * (" +
                 a.text() + ") * x + (" + b.text() +
                 "))]\n  · intro h\n    exact False.elim h\n\n";
      } else {
        std::string s = "Real.sqrt (" + d.text() + ")", r1, r2;
        try {
          auto sqrtD = squareRoot(d);
          s = sqrtD.text();
          rationalRoots = {(Rational(0) - b + sqrtD) / (Rational(2) * a),
                           (Rational(0) - b - sqrtD) / (Rational(2) * a)};
          r1 = rationalRoots[0].text();
          r2 = rationalRoots[1].text();
        } catch (Unsupported &) {
          exactRoots = false;
          r1 = "((-(" + b.text() + ") + " + s + ") / (2 * (" + a.text() + ")))";
          r2 = "((-(" + b.text() + ") - " + s + ") / (2 * (" + a.text() + ")))";
        }
        rootCondition = "x = " + r1 + " ∨ x = " + r2;
        step("Take both square-root branches and divide by 2a", rootCondition);
        proof += "theorem roots (x : ℝ) : residual x = 0 ↔ " + rootCondition +
                 " := by\n  rw [normalized]\n";
        if (!exactRoots)
          proof += "  have hs : (Real.sqrt (" + d.text() +
                   "))^2 = " + d.text() + " := Real.sq_sqrt (by norm_num)\n";
        proof += "  constructor\n  · intro h\n    have hf : (x - (" + r1 +
                 ")) * (x - (" + r2 +
                 ")) = 0 := by\n      ring_nf at h ⊢ <;> nlinarith" +
                 std::string(exactRoots ? "" : " [hs]") +
                 "\n    rcases mul_eq_zero.mp hf with h1 | h2\n    · left; "
                 "linarith\n    · right; linarith\n  · rintro (rfl | rfl) <;> "
                 "norm_num <;> nlinarith" +
                 std::string(exactRoots ? "" : " [hs]") + "\n\n";
      }
    } else {
      rootCondition = c.zero() ? "True" : "False";
      step("Resolve the constant equation",
           c.zero() ? "Every value in the domain satisfies the equation."
                    : "No values satisfy the equation.");
      proof += "theorem roots (x : ℝ) : residual x = 0 ↔ " + rootCondition +
               " := by\n  rw [normalized]\n  norm_num\n\n";
    }
    if (!rationalRoots.empty()) {
      std::sort(rationalRoots.begin(), rationalRoots.end(),
                [](Rational x, Rational y) { return x < y; });
      std::vector<Rational> kept;
      for (auto r : rationalRoots)
        if (!(r < lo) && !(hi < r) &&
            (kept.empty() || kept.back().text() != r.text()))
          kept.push_back(r);
      std::string allowed = kept.empty() ? "False" : "";
      for (auto r : kept) {
        if (!allowed.empty())
          allowed += " ∨ ";
        allowed += "x = " + r.text();
      }
      step("Apply the supplied domain [" + lo.text() + ", " + hi.text() + "]",
           kept.empty() ? "No roots lie in this domain." : allowed);
      proof += "theorem domain_solution (x : ℝ) : (" + lo.text() +
               " ≤ x ∧ x ≤ " + hi.text() + " ∧ residual x = 0) ↔ (" + allowed +
               ") := by\n  rw [roots]\n  constructor\n  · rintro ⟨hl, hu, h⟩\n";
      if (quadratic)
        proof +=
            "    rcases h with h | h\n    all_goals subst x; norm_num at *\n";
      else
        proof += "    subst x; norm_num at *\n";
      proof += "  · intro h\n";
      if (kept.empty())
        proof += "    exact False.elim h\n";
      else if (kept.size() == 1)
        proof += "    subst x; norm_num\n";
      else
        proof += "    rcases h with h | h\n    all_goals subst x; norm_num\n";
      proof += "\n#print axioms domain_solution\n";
    } else if (!exactRoots)
      step(
          "Restrict the roots to the supplied interval",
          "Keep each exact root only if " + lo.text() + " ≤ x ≤ " + hi.text() +
              ". No floating-point approximation is used in the root formula.");
    proof +=
        "#print axioms normalized\n#print axioms roots\nend StudyDerivation\n";
    result.set("kind", "exact_algebraic");
    result.set("leanSource", proof);
    result.set("proofScope",
               "Generated Lean proof for the exact real-number equation and "
               "its algebraic roots. Browser status: not kernel checked. Use "
               "--check-lean to check the exported proof.");
  } catch (algebra::Unsupported &) {
    if (p.query.monteCarlo) {
      result.set("kind", "simulation");
      step("Expand the supplied random model and event", readable(p.event));
      step("Estimate the event probability",
           "p̂ = hits / R; radius = sqrt(log(2K / (1 − confidence)) / (2R)).");
      step("Apply simultaneous bounds",
           "Compare [max(0,p̂−radius), min(1,p̂+radius)] with the target at each "
           "grid point. Report a minimum, a range, or unresolved.");
    } else if (p.query.evaluation == "analytic_gaussian_probability") {
      result.set("kind", "analytic_probability");
      step("Reduce the Gaussian event",
           "An affine Gaussian statistic is characterized by its mean and "
           "variance, retaining covariance from shared draws. Standardize the "
           "rejection boundary and evaluate the appropriate normal tails.");
      for (const auto &stage : p.analyticStages)
        step(stage.first, readableQuery(stage.second, p.query.variable));
      step("Derived probability expression (numeric constants are approximate)",
           readableQuery(p.expression, p.query.variable));
      step("Solve the resulting condition",
           p.query.evaluate ? "Evaluate the derived probability."
           : p.query.integer
               ? "Test each grid point in ascending order against the target; "
                 "the probability formula is analytic, while the inversion is "
                 "numerical."
               : "Bisect the stated interval until the equation residual meets "
                 "the tolerance. This is numerical inversion, not a "
                 "closed-form solution.");
    } else {
      step("Expand the user-defined objective", readable(p.expression));
      step("Available solution method",
           p.query.evaluate ? "Evaluate this expression."
           : p.query.integer
               ? "Evaluate grid points in ascending order and retain the first "
                 "qualifying value."
               : "No supported exact algebraic derivation is available. Use "
                 "bracketed bisection and report the residual; quadrature "
                 "errors require a separate convergence check.");
    }
  }
  result.set("steps", steps);
  return result;
}
} // namespace

J compile_program(const J &request) {
  for (const auto &field : request.items())
    if (field.first != "op" && field.first != "source" &&
        field.first != "solveFor")
      throw std::invalid_argument(
          "Program solve settings belong in the source: " + field.first);
  if (request.find("solveFor") && request.get("solveFor").as_string() != "code")
    throw std::invalid_argument(
        "Choose 'Query written in code' for a study program.");
  if (!request.get("source").is_string())
    throw std::invalid_argument("Supply study source text.");
  Program p = Compiler(request.get("source").as_string()).compile();
  J plan = obj({{"schema", "ssc.expression-plan"},
                {"version", 1},
                {"source", p.source},
                {"ir", ir(p)}});
  J steps = J::array();
  steps.push("Parse and type-check the user definitions; expand function calls "
             "while retaining local bindings.");
  steps.push(p.reason);
  std::string solver =
      p.query.evaluate ? "Evaluate the expression."
      : p.query.monteCarlo
          ? "Scan the specified grid using simultaneous Monte Carlo bounds. "
            "Return a minimum or an uncertainty range."
      : p.query.integer ? "Scan the grid in ascending order; no monotonicity "
                          "assumption is needed."
      : !p.query.polynomial.empty()
          ? "Try the algebraic linear/quadratic roots and check them against "
            "the original expression."
          : "Bisect the supplied sign-changing interval and check the equation "
            "residual. Continuity is required for the bracketing "
            "interpretation.";
  steps.push(solver);
  steps.push("Emit the expression as standalone C++17. Browser Solve evaluates "
             "the same typed expression plan in C++/WASM.");
  std::string expression =
      readableQuery(p.expression ? p.expression : p.event, p.query.variable);
  if (expression.size() > 12000)
    expression = expression.substr(0, 12000) + "…";
  return obj(
      {{"plan", plan},
       {"cpp", generate(p)},
       {"compiler", "ssc-study-expression-v2"},
       {"leanKernelChecked", false},
       {"language", "Lean-style study language with compiler-defined "
                    "solve/probability syntax; headers optional"},
       {"execution",
        "Standalone output compiles to a native C++ executable. Browser Solve "
        "executes the typed plan in the precompiled C++/WASM runtime."},
       {"query", queryWire(p.query)},
       {"derivation", derive(p)},
       {"explanation", obj({{"steps", steps},
                            {"probabilityQuantity", p.probabilityQuantity},
                            {"condition", p.condition},
                            {"loweredExpression", expression},
                            {"reason", p.reason}})}});
}
J execute_program(const J &plan, Token *token) {
  if (plan.get("schema").as_string() != "ssc.expression-plan" ||
      plan.get("version").as_number() != 1 || !plan.get("source").is_string())
    throw std::invalid_argument("Expected ssc.expression-plan version 1.");
  for (const auto &field : plan.items())
    if (field.first != "schema" && field.first != "version" &&
        field.first != "source" && field.first != "ir")
      throw std::invalid_argument("Unexpected expression-plan field: " +
                                  field.first);
  Program p = Compiler(plan.get("source").as_string()).compile();
  if (ir(p).dump() != plan.get("ir").dump())
    throw std::invalid_argument(
        "Compiled plan and source differ; compile the study again.");
  rt::Context ctx;
  if (token) {
    token->check();
    ctx.check = [&]() { token->check(); };
  }
  Env env;
  auto expression = [&](double x, rt::Context &c) {
    env[0] = x;
    return real(eval(p.expression, env, c));
  };
  auto event = [&](double x, rt::Context &c) {
    env[0] = x;
    return boolean(eval(p.event, env, c));
  };
  auto report = rt::run(p.query, expression, event, ctx);
  report.put("operations", double(ctx.operations));
  if (token)
    token->check();
  return json::parse(report.json());
}
} // namespace ssc::design
