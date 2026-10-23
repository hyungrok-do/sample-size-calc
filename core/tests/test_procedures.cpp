// Procedure regression tests: run each case of fixtures/procedures.json
// through the engine and compare with reference values computed independently
// of it (the fixture records the software that computed them).  Also checks
// engine behavior that does not depend on reference values: thread-count
// invariance, cancellation and input validation.
#include <fstream>
#include <new>
#include <sstream>

#include "engine/parse.hpp"
#include "engine/runner.hpp"
#include "testing.hpp"
#include "util/json.hpp"
#include "util/parallel.hpp"

using namespace ssc;

namespace {

// The reference cases.
json::Value reference_cases() {
  std::ifstream in(std::string(SSC_FIXTURE_DIR) + "/procedures.json");
  std::stringstream ss;
  ss << in.rdbuf();
  return json::parse(ss.str()).get("cases");
}

json::Value run(const std::string& request, Token* token = nullptr) {
  Token local;
  return run_request(json::parse(request), token ? token : &local);
}

json::Value request_of(const json::Value& c, int threads, bool curve) {
  json::Value req = json::Value::object();
  req.set("procedure", c.get("procedure").as_string());
  req.set("solve_for", c.get("solve_for").as_string());
  req.set("inputs", c.get("inputs"));
  json::Value opts = json::Value::object();
  opts.set("threads", static_cast<double>(threads));
  opts.set("power_curve", curve);
  req.set("options", opts);
  return req;
}

}  // namespace

TEST_CASE("procedures: independent reference values") {
  const json::Value cases = reference_cases();
  CHECK(cases.size() > 40);
  for (const auto& c : cases.elements()) {
    const std::string name = c.get("name").as_string();
    Token token;
    const json::Value res = run_request(request_of(c, 0, false), &token);
    if (!res.get("ok").as_bool()) {
      ::ssc::test::count_check();
      ::ssc::test::report_failure(name + ": " + res.dump(), __FILE__, __LINE__);
      continue;
    }
    const double rtol = c.get("rtol").as_number(1e-8);
    for (const auto& e : c.get("expect").elements()) {
      const auto& row = res.get("rows").at(static_cast<std::size_t>(e.get("row").as_number()));
      const std::string key = e.get("key").as_string();
      const double got = row.get(key).as_number(NAN), want = e.get("value").as_number();
      ::ssc::test::count_check();
      if (!::ssc::test::close(got, want, rtol, 0)) {
        std::ostringstream os;
        os.precision(12);
        os << name << ": " << key << " = " << got << ", expected " << want << " ("
           << c.get("source").as_string() << ")";
        if (row.find("_error")) os << " error: " << row.get("_error").as_string();
        ::ssc::test::report_failure(os.str(), __FILE__, __LINE__);
      }
    }
  }
}

TEST_CASE("procedures: results do not depend on the number of threads") {
  const json::Value cases = reference_cases();
  int n = 0;
  for (const auto& c : cases.elements()) {
    if (++n % 3 != 0) continue;  // a sample keeps the test quick
    Token t1, t4;
    const json::Value a = run_request(request_of(c, 1, true), &t1);
    const json::Value b = run_request(request_of(c, 4, true), &t4);
    // identical rows and power curves (bitwise: no reductions depend on scheduling)
    CHECK_EQ(a.get("rows").dump(), b.get("rows").dump());
    CHECK_EQ(a.get("plots").dump(), b.get("plots").dump());
  }
}

TEST_CASE("procedures: input validation and errors") {
  const json::Value bad = run(R"({"procedure":"means.two.t","solve_for":"n","inputs":{"power":"1.5"}})");
  CHECK(!bad.get("ok").as_bool());
  CHECK(bad.get("input_errors").find("power") != nullptr);

  const json::Value unknown = run(R"({"procedure":"no.such.procedure","solve_for":"n","inputs":{}})");
  CHECK(!unknown.get("ok").as_bool());

  // a scenario that cannot be computed reports a per-row error, the others still run
  const json::Value mixed = run(R"({"procedure":"means.two.t","solve_for":"n","inputs":{"mu1":"0 1","mu2":"0"}})");
  CHECK(mixed.get("ok").as_bool());
  CHECK_EQ(mixed.get("rows").size(), std::size_t{2});
  CHECK(mixed.get("rows").at(0).find("_error") != nullptr);
  CHECK(mixed.get("rows").at(1).find("_error") == nullptr);

  // lists and ranges expand into scenarios
  const json::Value grid = run(
      R"({"procedure":"means.two.t","solve_for":"power","inputs":{"n1":"10 to 50 by 10","mu1":"0.5 1"}})");
  CHECK_EQ(grid.get("rows").size(), std::size_t{10});
}

TEST_CASE("procedures: a scenario without the memory it needs says so") {
  // a procedure of the test's own whose scenario at x = 2 cannot have the memory it asks for, as one past the 2 GB of
  // the page's engine cannot
  Procedure memory;
  memory.id = "test.memory";
  memory.name = "Memory";
  memory.solve_for = {{"power", "Power", "", ""}};
  memory.params = {num("x", "X", "1"), num("n", "N", "10")};
  memory.curve_param = "n";  // and so are the points of its power curve past n = 15
  memory.run = [](const Scenario& s, Row& r, RunContext&) {
    if (s("x") == 2 || s("n") > 15) throw std::bad_alloc();
    r.set("power", 0.5);
  };
  registry().push_back(&memory);
  const json::Value res =
      run(R"({"procedure":"test.memory","solve_for":"power","inputs":{"x":"1 2"},"options":{"power_curve":false}})");
  const json::Value curve =
      run(R"({"procedure":"test.memory","solve_for":"power","inputs":{"x":"1"},"options":{"power_curve":true}})");
  registry().pop_back();
  CHECK(res.get("ok").as_bool(false));
  CHECK_EQ(res.get("rows").size(), std::size_t{2});
  CHECK(res.get("rows").at(0).find("_error") == nullptr);
  CHECK_EQ(res.get("rows").at(1).get("_error").as_string(),
           std::string("there is not enough memory to compute the answer for these inputs"));
  CHECK(res.get("notes").size() == 0);
  CHECK(curve.get("ok").as_bool(false));
  CHECK(curve.get("rows").at(0).find("_error") == nullptr);
  CHECK_EQ(curve.get("notes").at(0).as_string(),
           std::string("Some points of the power curve are left out: there was not enough memory to compute them."));
  const json::Value& xs = curve.get("plots").at(0).get("series").at(0).get("x");
  CHECK(xs.size() > 0);
  for (const auto& x : xs.elements()) CHECK(x.as_number(NAN) <= 15);
}

TEST_CASE("procedures: typed and pasted numbers") {
  auto message = [](const std::string& mu1) {
    const json::Value r = run(R"({"procedure":"means.two.t","solve_for":"n","inputs":{"mu1":")" + mu1 + R"("}})");
    return r.get("ok").as_bool() ? std::string() : r.get("input_errors").get("mu1").as_string();
  };
  auto first_mu1 = [](const std::string& mu1) {
    const json::Value r = run(R"({"procedure":"means.two.t","solve_for":"n","inputs":{"mu1":")" + mu1 + R"("}})");
    return r.get("rows").at(0).get("mu1").as_number();
  };
  // "10-20" read as 10 and -20 would quietly make two scenarios; a range is written "10 to 20"
  CHECK(message("10-20").find("'10 to 20'") != std::string::npos);
  CHECK(message("0.5-1").find("'0.5 to 1'") != std::string::npos);
  CHECK(message("3+4").find("as a number") != std::string::npos);
  CHECK(message("80%").find("without '%'") != std::string::npos);
  CHECK_EQ(message("1e-5 -1 -2"), std::string());  // exponents and negative values in a list
  // the minus sign, and the spaces of pasted text (no-break, thin, narrow no-break)
  CHECK_EQ(first_mu1("\u22120.5"), -0.5);
  CHECK_EQ(first_mu1("0.75\u00a0"), 0.75);
  CHECK_EQ(first_mu1("\u2009 0.25\u202f"), 0.25);
  // a long text is refused before it is read; the limit counts characters, not bytes
  CHECK(message(std::string(200001, '1')).find("too long") != std::string::npos);
  std::string accents;
  for (int i = 0; i < 200000; ++i) accents += "\u00e9";
  CHECK(message(accents).find("too long") == std::string::npos);
  // a range that crosses or reaches 0 has 0 among its values, not the rounding noise of a + k * step, at any scale
  // (0.9 - 3 * 0.3 is 1.1e-16, 0.3 - 3 * 0.1 is -5.6e-17)
  CHECK_EQ(parse::numbers("-0.3 to 0.3 by 0.1").values.at(3), 0.0);
  CHECK_EQ(parse::numbers("-3e-15 to 3e-15 by 1e-15").values.at(3), 0.0);
  CHECK_EQ(parse::numbers("0.3e-29 to -0.3e-29 by 1e-30").values.at(3), 0.0);
  CHECK_EQ(parse::numbers("0.9 to 0 by 0.3").values.at(3), 0.0);
  CHECK_EQ(parse::numbers("0.3 to 0 by 0.1").values.at(3), 0.0);
  // its start is kept as typed, and so is its end where the steps reach it, while the values between are rounded to
  // 12 digits; both ends as typed are kept besides, reached or not
  CHECK_EQ(parse::numbers("1e-12 to 1 by 0.5").values.at(0), 1e-12);
  CHECK(parse::numbers("0.99999999999999 to 0.99999999999999").values == std::vector<double>({0.99999999999999}));
  CHECK(parse::numbers("0.1234567890123456 to 1 by 0.5").values ==
        std::vector<double>({0.1234567890123456, 0.623456789012}));
  CHECK(parse::numbers("0.9 to 0.99999999999999 by 0.09999999999999").values ==
        std::vector<double>({0.9, 0.99999999999999}));
  CHECK(parse::numbers("0 to 0.3 by 0.1").values == std::vector<double>({0, 0.1, 0.2, 0.3}));
  CHECK(parse::numbers("0 to 1.5 by 1").ends == std::vector<double>({0, 1.5}));
  CHECK(parse::vectors("1 2; 0 to 1.5 by 1").ends == std::vector<double>({0, 1.5}));
  // characters that would act on the message are written out: a right-to-left override, a joiner, NUL, a
  // combining mark
  const std::string rlo = message("\u202e0.5");
  CHECK(rlo.find("\\u{202E}0.5") != std::string::npos);
  CHECK(rlo.find("\u202e") == std::string::npos);
  CHECK(message("1\\u0000x").find("\\u{0}") != std::string::npos);  // NUL, escaped in the request
  CHECK(message("0.0\u0301").find("\\u{301}") != std::string::npos);
  CHECK(message("1\u2066x").find("\\u{2066}") != std::string::npos);
  // the zero-width characters show nothing and are read as nothing: "0<ZWSP>.4" is the 0.4 the field shows
  CHECK_EQ(first_mu1("0\u200b.4"), 0.4);
  CHECK_EQ(first_mu1("1\u200d2\u2060"), 12.0);
  CHECK_EQ(parse::quoted("abc\u00e9"), std::string("abc\u00e9"));
  // marks of any script, and enclosing ones, are written out; letters and symbols of common scripts are not
  CHECK_EQ(parse::quoted("\u05b0abc"), std::string("\\u{5B0}abc"));
  CHECK_EQ(parse::quoted("\u0489\u064e\u093c"), std::string("\\u{489}\\u{64E}\\u{93C}"));
  CHECK_EQ(parse::quoted("\u03bc \u2264 \u00bd \u4e2d"), std::string("\u03bc \u2264 \u00bd \u4e2d"));
  // the digits of other scripts, with the Arabic separators
  CHECK_EQ(parse::quoted("\u06f0\u066b\u06f0\u06f5 \u0665\u066a \u0662\u066c\u0665 \u0e50.\u0e55 \u09e6.\u09eb"),
           std::string("\u06f0\u066b\u06f0\u06f5 \u0665\u066a \u0662\u066c\u0665 \u0e50.\u0e55 \u09e6.\u09eb"));
  // bytes that are not UTF-8 stay bytes (the JSON writer shows U+FFFD), never a character they could spell
  CHECK_EQ(parse::quoted("a\x85" "b"), std::string("a\x85" "b"));
  CHECK_EQ(parse::quoted("a\xc0\x80" "b"), std::string("a\xc0\x80" "b"));
  CHECK_EQ(parse::quoted("a\xf0\x82\x80\xae" "b"), std::string("a\xf0\x82\x80\xae" "b"));
  CHECK_EQ(parse::quoted("a\xe2\x80"), std::string("a\xe2\x80"));
  // the cut counts the text as written: 37 controls give at most 40 bytes and the ellipsis
  CHECK(parse::quoted(std::string(37, '\x1f')).size() <= 40 + std::string("\u2026").size());
  CHECK_EQ(parse::quoted(std::string(40, 'x')), std::string(40, 'x'));
  CHECK_EQ(parse::quoted(std::string(41, 'x')), std::string(40, 'x') + "\u2026");
  // objects and arrays nested beyond any request are refused, not followed into a stack overflow
  std::string deep = R"({"procedure":"means.two.t","inputs":{"x":)";
  for (int i = 0; i < 100000; ++i) deep += '[';
  bool refused = false;
  try {
    json::parse(deep);
  } catch (const std::exception& e) {
    refused = std::string(e.what()).find("nested too deeply") != std::string::npos;
  }
  CHECK(refused);
}

TEST_CASE("procedures: cancellation stops a running calculation") {
  Token token;
  token.cancel();
  const json::Value res = run(
      R"({"procedure":"prop.two","solve_for":"power",)"
      R"("inputs":{"n1":"500 to 5000 by 5","test":"fisher","p1":"0.5","p2":"0.45"}})",
      &token);
  CHECK(!res.get("ok").as_bool());
  CHECK(res.get("cancelled").as_bool());
}

TEST_CASE("procedures: every procedure runs its example") {
  const json::Value cat = catalog_json();
  for (const auto& p : cat.get("procedures").elements()) {
    const std::string id = p.get("id").as_string();
    const json::Value schema = json::parse(procedure_json(*find_procedure(id)).dump());
    const json::Value& example = schema.get("example");
    json::Value req = json::Value::object();
    req.set("procedure", id);
    req.set("solve_for", example.get("solve_for").as_string());
    req.set("inputs", example.get("inputs"));
    Token token;
    const json::Value res = run_request(req, &token);
    ::ssc::test::count_check();
    bool ok = res.get("ok").as_bool();
    for (const auto& r : res.get("rows").elements())
      if (r.find("_error")) ok = false;
    if (!ok) ::ssc::test::report_failure(id + " example failed: " + res.dump().substr(0, 400), __FILE__, __LINE__);
  }
}

TEST_CASE("procedures: every input has a help text") {
  // The app explains each input with its help text: one or two short, complete sentences.
  const json::Value cat = catalog_json();
  int params = 0;
  for (const auto& p : cat.get("procedures").elements()) {
    const std::string id = p.get("id").as_string();
    const json::Value schema = procedure_json(*find_procedure(id));
    for (const auto& q : schema.get("params").elements()) {
      ++params;
      const std::string where = id + " " + q.get("id").as_string() + ": ";
      const std::string help = q.get("help").as_string();
      std::size_t chars = 0;  // Unicode code points of the UTF-8 text
      for (char c : help)
        if ((static_cast<unsigned char>(c) & 0xC0) != 0x80) ++chars;
      ::ssc::test::count_check();
      if (help.empty()) {
        ::ssc::test::report_failure(where + "no help text", __FILE__, __LINE__);
        continue;
      }
      if (help.back() != '.')
        ::ssc::test::report_failure(where + "the help text does not end with a period: " + help, __FILE__, __LINE__);
      if (chars > 300)
        ::ssc::test::report_failure(where + "the help text has " + std::to_string(chars) + " characters (at most 300)",
                                    __FILE__, __LINE__);
      if (help.front() == ' ' || help.find("  ") != std::string::npos)
        ::ssc::test::report_failure(where + "stray spaces in the help text: " + help, __FILE__, __LINE__);
    }
  }
  CHECK(params > 700);
}
