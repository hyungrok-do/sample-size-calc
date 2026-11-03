#include "model.hpp"
#include <algorithm>
#include <cmath>
#include <set>

namespace ssc::dag {
Transform parse_transform(const J &v) {
  Transform t;
  if (v.is_null()) return t;
  if (!v.is_object()) fail("A transformation must be an object.");
  t.kind = v.get("kind").as_string();
  if (t.kind != "identity" && t.kind != "square" && t.kind != "threshold" &&
      t.kind != "hinge" && t.kind != "saturation" && t.kind != "spline" &&
      t.kind != "log")
    fail("Unknown response transformation.");
  if (v.find("center")) t.center = number(v, "center");
  if (v.find("scale")) t.scale = number(v, "scale");
  if (!(t.scale > 0)) fail("Transformation scale must be positive.");
  if (v.find("knot")) t.knot = number(v, "knot");
  if (t.kind == "spline") {
    const auto &k = v.get("knots");
    if (!k.is_array() || k.size() < 3 || k.size() > 7)
      fail("A restricted cubic spline needs 3–7 ordered, distinct knots.");
    for (const auto &x : k.elements()) {
      if (!x.is_number() || !std::isfinite(x.as_number()) ||
          (!t.knots.empty() && x.as_number() <= t.knots.back()))
        fail("Spline knots must be finite and strictly increasing.");
      t.knots.push_back(x.as_number());
    }
    t.index = integer(v, "index", 0, static_cast<int>(t.knots.size()) - 3, 0);
  }
  return t;
}
J transform_json(const Transform &t) {
  J v = object({{"kind", t.kind}});
  if (t.kind == "square" || t.kind == "saturation") {
    v.set("center", t.center); v.set("scale", t.scale);
  }
  if (t.kind == "threshold" || t.kind == "hinge") v.set("knot", t.knot);
  if (t.kind == "spline") { v.set("knots", array(t.knots)); v.set("index", t.index); }
  return v;
}
double transform_value(double x, const Transform &t) {
  if (t.kind == "identity") return x;
  if (t.kind == "square") { const double z = (x - t.center) / t.scale; return z * z; }
  if (t.kind == "saturation") return std::tanh((x - t.center) / t.scale);
  if (t.kind == "threshold") return x > t.knot ? 1 : 0;
  if (t.kind == "hinge") return std::max(0.0, x - t.knot);
  if (t.kind == "log") {
    if (!(x > 0))
      throw Error("generation_error", "log(x) needs positive values; a parent of a log term took a value of zero or less.");
    return std::log(x);
  }
  const auto &k = t.knots;
  const double a = k[t.index], b = k[k.size() - 2], c = k.back(),
               width = c - k.front(), z = std::min(x, c),
               u = std::max(0.0, z - a) / width,
               v = std::max(0.0, z - b) / width,
               w = std::max(0.0, z - c) / width,
               r = (c - a) / (c - b), s = (b - a) / (c - b);
  // Restricted cubic basis, divided by the squared boundary-knot span.
  // Linear tails are evaluated directly to avoid cancellation for large x.
  double value = width * (u*u*u - r*v*v*v + s*w*w*w);
  if (x > c) value += 3 * (u*u - r*v*v) * (x - c);
  return value;
}
double term_value(const Term &t, double left, double right) {
  const double value = transform_value(left, t.transform) *
      (t.right < 0 ? 1 : transform_value(right, t.rightTransform));
  if (!std::isfinite(value)) throw Error("numerical", "A transformed term overflowed.");
  return value;
}
std::string term_key(const Term &t) {
  std::string a=std::to_string(t.left)+":"+transform_json(t.transform).dump();
  if(t.right<0) return a;
  std::string b=std::to_string(t.right)+":"+transform_json(t.rightTransform).dump();
  return a<b ? a+"*"+b : b+"*"+a;
}
bool nonlinear_in(const Term &t, int v) {
  return (t.left == v && t.transform.kind != "identity") ||
         (t.right == v && t.rightTransform.kind != "identity") ||
         (t.left == v && t.right == v);
}
std::string transform_text(const std::string &id, const Transform &t) {
  if (t.kind == "identity") return id;
  if (t.kind == "square" || t.kind == "saturation") {
    const auto z = "(" + id + " - " + J(t.center).to_text() + ") / " + J(t.scale).to_text();
    return t.kind == "square" ? "square(" + z + ")" : "tanh(" + z + ")";
  }
  if (t.kind == "threshold") return "I(" + id + " > " + J(t.knot).to_text() + ")";
  if (t.kind == "hinge") return "max(0, " + id + " - " + J(t.knot).to_text() + ")";
  if (t.kind == "log") return "log(" + id + ")";
  return "rcs" + std::to_string(t.index + 1) + "(" + id + "; " + array(t.knots).dump() + ")";
}
Term parse_term(const Graph &g, const J &v, bool coefficient) {
  Term t{g.index(v.get("left").as_string()), -1, coefficient ? number(v,"coefficient") : 0};
  if(v.find("right") && !v.get("right").is_null() && !v.get("right").is_string()) fail("A response term needs parent identifiers.");
  if (v.find("right") && !v.get("right").is_null() && !v.get("right").as_string().empty())
    t.right = g.index(v.get("right").as_string());
  t.transform = parse_transform(v.get("transform"));
  t.rightTransform = parse_transform(v.get("rightTransform"));
  return t;
}
J term_json(const Graph &g, const Term &t, bool coefficient) {
  J v = object({{"left", g.nodes[t.left].id}});
  if (t.right >= 0) v.set("right", g.nodes[t.right].id);
  if (t.transform.kind != "identity") v.set("transform", transform_json(t.transform));
  if (t.right >= 0 && t.rightTransform.kind != "identity") v.set("rightTransform", transform_json(t.rightTransform));
  if (coefficient) v.set("coefficient", t.coefficient);
  return v;
}
static J basis(const J &raw) {
  J out = raw.is_null() ? J::object() : raw;
  if (!out.is_object()) fail("An analysis model must be an object.");
  const auto b = (out.find("basis") ? out.get("basis").as_string() : "match");
  if (b != "match" && b != "linear" && b != "quadratic" && b != "spline" && b != "custom")
    fail("Choose matching, linear, quadratic, spline or custom analysis terms.");
  out.set("basis", b);
  if (out.find("interactions") && !out.get("interactions").is_bool()) fail("Analysis interactions must be true or false.");
  if (b == "spline") {
    const auto knots = out.find("knots") ? out.get("knots") : array({-1,0,1});
    parse_transform(object({{"kind","spline"},{"knots",knots}}));
    out.set("knots",knots);
  }
  if (b == "custom" && (!out.get("terms").is_array() || out.get("terms").size() > 120))
    fail("A custom analysis needs an explicit list of at most 120 terms.");
  return out;
}
J parse_analysis(const J &raw) {
  J out = basis(raw);
  const auto method = (out.find("method") ? out.get("method").as_string() : "auto");
  if (method != "auto" && method != "gformula" && method != "aipw" && method != "tmle" && method != "mmrm" && method != "lmm" && method != "trajectory")
    fail("Choose automatic analysis, g-computation, AIPW or binary-outcome TMLE.");
  out.set("method",method);
  const int folds = integer(out,"folds",1,10,1);
  out.set("folds",folds);
  const double bound = out.find("propensityBound") ? number(out,"propensityBound") : 0;
  if (!(bound>=0 && bound<.5)) fail("The propensity bound must be in [0, 0.5).");
  out.set("propensityBound",bound);
  if (out.find("outcomeBounds")) {
    const auto &b=out.get("outcomeBounds");
    if (!b.is_array() || b.size()!=2 || !b.at(0).is_number() || !b.at(1).is_number() ||
        !std::isfinite(b.at(0).as_number()) || !std::isfinite(b.at(1).as_number()) ||
        !(b.at(0).as_number()<b.at(1).as_number()))
      fail("Supply finite, increasing prespecified outcome bounds.");
  }
  if (out.find("learner")) {
    const auto learner = out.get("learner").as_string();
    if (learner != "glm" && learner != "superlearner") fail("Choose parametric or Super Learner nuisance models.");
  }
  if (out.find("propensity")) out.set("propensity",basis(out.get("propensity")));
  if (out.find("models")) {
    if (!out.get("models").is_object() || out.get("models").size()>40) fail("Invalid analysis model overrides.");
    J models=J::object();
    for(const auto &kv:out.get("models").items()) models.set(kv.first,basis(kv.second));
    out.set("models",models);
  }
  return out;
}
std::vector<Term> analysis_terms(const Graph &g, const std::vector<Term> &terms,
                                 const J &analysis, int node, bool propensity) {
  const auto &models=analysis.get("models");
  const J &cfg=models.find(g.nodes[node].id) ? models.get(g.nodes[node].id) :
      propensity && analysis.find("propensity") ? analysis.get("propensity") : analysis;
  const auto mode=(cfg.find("basis") ? cfg.get("basis").as_string() : "match");
  if(mode=="match") return terms;
  std::vector<Term> out;
  if(mode=="custom") {
    for(const auto &t:cfg.get("terms").elements()) {
      out.push_back(parse_term(g,t,false));
      for(int v:{out.back().left,out.back().right})
        if(v>=0 && !g.nodes[v].measured) fail("Analysis terms cannot use "+g.nodes[v].id+": it is unmeasured.");
    }
  } else {
    std::set<int> variables;
    for(const auto &t:terms) { variables.insert(t.left); if(t.right>=0) variables.insert(t.right); }
    for(int v:variables) {
      out.push_back({v,-1,0});
      if(g.nodes[v].kind=="binary") continue;
      if(mode=="quadratic") { Term t{v,-1,0}; t.transform.kind="square"; out.push_back(t); }
      if(mode=="spline") {
        const J knots=cfg.find("knots") ? cfg.get("knots") : array({-1,0,1});
        for(int j=0;j<static_cast<int>(knots.size())-2;++j) {
          Term t{v,-1,0}; t.transform=parse_transform(object({{"kind","spline"},{"knots",knots},{"index",j}})); out.push_back(t);
        }
      }
    }
    if(cfg.get("interactions").as_bool(false))
      for(int a:variables) for(int b:variables) if(a<b) out.push_back({a,b,0});
  }
  if(out.size()>120) fail("An analysis model supports at most 120 terms.");
  return out;
}
} // namespace ssc::dag
