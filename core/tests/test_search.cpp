// Procedure search: text normalization, the bounded edit distance, ranking on
// the relevance benchmark of fixtures/search.json (queries a user might type,
// with the procedure that should come first or among the first few), spelling
// corrections, and that every procedure is found by its own name and id.
#include <algorithm>
#include <chrono>
#include <fstream>
#include <sstream>

#include "engine/procedure.hpp"
#include "engine/search.hpp"
#include "testing.hpp"
#include "util/json.hpp"
#include "util/parallel.hpp"

using namespace ssc;

namespace {

bool has(const std::vector<std::string>& v, const std::string& x) {
  return std::find(v.begin(), v.end(), x) != v.end();
}

std::vector<std::string> ids(const search::Result& r) {
  std::vector<std::string> out;
  for (const auto& h : r.hits) out.push_back(h.id);
  return out;
}

int rank_of(const search::Result& r, const std::string& id) {
  for (std::size_t i = 0; i < r.hits.size(); ++i)
    if (r.hits[i].id == id) return static_cast<int>(i) + 1;
  return 0;
}

json::Value benchmark() {
  std::ifstream in(std::string(SSC_FIXTURE_DIR) + "/search.json");
  std::stringstream ss;
  ss << in.rdbuf();
  return json::parse(ss.str());
}

}  // namespace

TEST_CASE("search: text normalization") {
  const auto t = search::terms("Log-Rank Tests");
  CHECK(has(t, "log") && has(t, "rank") && has(t, "logrank") && has(t, "test"));
  CHECK(has(search::terms("Cohen's κ"), "cohen"));  // apostrophe dropped, then the plural-like s
  CHECK(has(search::terms("Cohen's κ"), "kappa"));
  CHECK(has(search::terms("χ² test"), "chi2"));
  CHECK(has(search::terms("ρ0"), "rho0"));
  CHECK(has(search::terms("Analyses"), "analysis"));
  CHECK(has(search::terms("proportions"), "proportion"));
  CHECK(has(search::terms("studies"), "study"));
  CHECK(has(search::terms("classes"), "class"));
  CHECK(has(search::terms("Équivalence naïve"), "equivalence"));
  CHECK(has(search::terms("Équivalence naïve"), "naive"));
  CHECK(search::terms("the power of a test") == (std::vector<std::string>{"power", "test"}));
  CHECK(search::terms("").empty());
  CHECK(search::terms(" ,;!? ").empty());
  CHECK(search::terms("\xff\xfe broken \xc3").size() == 1);  // invalid UTF-8 is skipped
}

TEST_CASE("search: bounded edit distance") {
  CHECK_EQ(search::edit_distance("kitten", "sitting", 5), 3);
  CHECK_EQ(search::edit_distance("mcnemar", "mcnemra", 2), 1);  // a transposition counts once
  CHECK_EQ(search::edit_distance("poisson", "poison", 2), 1);
  CHECK_EQ(search::edit_distance("same", "same", 1), 0);
  CHECK_EQ(search::edit_distance("abc", "xyz", 1), 2);  // limit + 1 once the limit is exceeded
  CHECK_EQ(search::edit_distance("short", "much longer", 2), 3);
}

TEST_CASE("search: relevance benchmark") {
  const json::Value cases = benchmark();
  CHECK(cases.size() >= 100);
  double reciprocal = 0;
  int ranked = 0;
  for (const auto& c : cases.elements()) {
    const std::string q = c.get("query").as_string();
    const search::Result r = search::query(q, 10);
    if (c.get("none").as_bool()) {
      if (!r.hits.empty()) ::ssc::test::report_failure("'" + q + "' should find nothing", __FILE__, __LINE__);
      ::ssc::test::count_check();
      continue;
    }
    const std::string want = c.get("expect").as_string();
    const int within = static_cast<int>(c.get("within").as_number(1));
    const int rank = rank_of(r, want);
    ::ssc::test::count_check();
    if (rank == 0 || rank > within) {
      std::string got;
      for (const auto& id : ids(r)) got += (got.empty() ? "" : ", ") + id;
      ::ssc::test::report_failure("'" + q + "': " + want + " at rank " + std::to_string(rank) + " (want ≤ " +
                                      std::to_string(within) + "); got " + got,
                                  __FILE__, __LINE__);
    }
    if (const json::Value* fix = c.find("corrected")) {
      ::ssc::test::count_check();
      if (r.corrected != fix->as_string())
        ::ssc::test::report_failure("'" + q + "' corrected to '" + r.corrected + "', want '" + fix->as_string() + "'",
                                    __FILE__, __LINE__);
    }
    reciprocal += rank ? 1.0 / rank : 0.0;
    ++ranked;
  }
  CHECK(reciprocal / ranked > 0.9);  // mean reciprocal rank of the expected procedure
}

TEST_CASE("search: every procedure is found by its name and id") {
  for (const Procedure* p : catalog_procedures()) {
    for (const std::string& q : {p->name, p->id}) {
      const search::Result r = search::query(q, 3);
      ::ssc::test::count_check();
      if (r.hits.empty() || r.hits[0].id != p->id)
        ::ssc::test::report_failure("'" + q + "' does not find " + p->id + " first", __FILE__, __LINE__);
    }
  }
}

TEST_CASE("search: corrections, prefixes, partial matches and limits") {
  // a misspelled word is matched and reported as corrected; a correctly spelled one is not
  CHECK_EQ(search::query("mcnemmar").corrected, std::string("mcnemar"));
  CHECK(search::query("mcnemar").corrected.empty());
  // a word still being typed matches as a prefix, without being "corrected"
  const search::Result typing = search::query("bioequ");
  CHECK(!typing.hits.empty() && typing.hits[0].id == "means.be.ratio");
  CHECK(typing.corrected.empty());
  // words that no procedure has together: the procedures with most of them, flagged
  const search::Result partial = search::query("survival phase ii");
  CHECK(partial.partial);
  CHECK(!partial.hits.empty());
  CHECK(!search::query("logrank survival").partial);
  // nothing to search for
  CHECK(search::query("").hits.empty());
  CHECK(search::query("   ").hits.empty());
  CHECK(search::query("the of and").hits.empty());
  CHECK(search::query("qqqqzzzz").hits.empty());
  // limits: 0 returns every match, ranked
  const search::Result all = search::query("test", 0);
  CHECK(all.hits.size() > 20);
  CHECK_EQ(search::query("test", 3).hits.size(), std::size_t{3});
  for (std::size_t i = 1; i < all.hits.size(); ++i) CHECK(all.hits[i - 1].score >= all.hits[i].score);
  // matched words are reported for each hit
  const search::Result km = search::query("kaplan meier");
  CHECK(!km.hits.empty() && !km.hits[0].terms.empty());
}

TEST_CASE("search: concurrent queries and speed") {
  const std::vector<std::string> queries = {"two sample t test", "logrnk", "NI proportions", "κ",
                                            "margin of error",   "icc",    "group sequential", "poison"};
  std::vector<std::vector<std::string>> expected;
  for (const auto& q : queries) expected.push_back(ids(search::query(q)));
  std::vector<int> same(64, 0);
  parallel_for(same.size(), 8, [&](std::size_t i) {
    same[i] = ids(search::query(queries[i % queries.size()])) == expected[i % queries.size()];
  });
  CHECK(std::all_of(same.begin(), same.end(), [](int s) { return s == 1; }));

  const auto start = std::chrono::steady_clock::now();
  int n = 0;
  for (int rep = 0; rep < 50; ++rep)
    for (const auto& q : queries) n += static_cast<int>(search::query(q).hits.size());
  const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
  CHECK(n > 0);
  CHECK(ms / (50.0 * queries.size()) < 5.0);  // typically a few hundredths of a millisecond
}
