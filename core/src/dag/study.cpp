#include "study.hpp"
#include "procedures/sim_mmrm.hpp"
#include "sim/lmm.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>

namespace ssc::dag {
namespace {
double optional(const J &v, const char *key, double fallback = 0) {
  return v.find(key) ? number(v,key) : fallback;
}
int level(const Node &n) { return n.scope=="cluster" ? 0 : n.scope=="subject" ? 1 : 2; }
double logistic(double v) { return v>=0 ? 1/(1+std::exp(-v)) : std::exp(v)/(1+std::exp(v)); }
int cluster_size(const Graph &g, int k) {
  const auto &sizes=g.study.clusterSizes;
  return sizes.empty() ? 1 : sizes[k % sizes.size()];
}
}
void parse_study(Graph &g) {
  const auto &raw=g.json.get("study");
  if(raw.is_null()) return;
  if(!raw.is_object()) fail("Study structure must be an object.");
  auto &s=g.study; s.enabled=true;
  const auto &times=raw.get("times");
  if(!times.is_array() || times.size()<1 || times.size()>12) fail("Specify 1–12 increasing visit times.");
  for(const auto &t:times.elements()) {
    if(!t.is_number() || !std::isfinite(t.as_number()) || (!s.times.empty() && t.as_number()<=s.times.back()))
      fail("Visit times must be finite and strictly increasing.");
    s.times.push_back(t.as_number());
  }
  if(raw.find("clusterSizes")) {
    if(!raw.get("clusterSizes").is_array() || raw.get("clusterSizes").size()>100) fail("Specify at most 100 cluster sizes.");
    for(const auto &v:raw.get("clusterSizes").elements()) {
      const double n=v.as_number(-1);
      if(n<1 || n>1000 || n!=std::floor(n)) fail("Cluster sizes must be integers from 1 to 1,000.");
      s.clusterSizes.push_back(static_cast<int>(n));
    }
  }
  std::set<std::string> outputs={"subject_id","cluster_id","visit_index","visit_time"};
  for(const auto &n:g.nodes) if(!outputs.insert(n.id).second) fail("A variable ID conflicts with a study identifier column.");
  for(int i=0;i<static_cast<int>(g.nodes.size());++i) {
    auto &n=g.nodes[i]; const auto &v=g.json.get("nodes").at(i);
    n.scope=v.find("scope") ? v.get("scope").as_string() : "subject";
    if(n.scope!="subject" && n.scope!="cluster" && n.scope!="visit" && n.scope!="time") fail("Choose cluster, subject, visit or time scope.");
    if(n.scope=="cluster" && s.clusterSizes.empty()) fail("Cluster variables require cluster sizes in the study structure.");
    if(n.scope=="time" && (n.family!="gaussian" || !n.terms.empty())) fail("A visit-time variable must be a continuous Normal root; its values come from the schedule.");
    if(n.scope=="time" && !n.measured) fail("The visit time is part of the schedule and is always measured.");
    if(n.scope=="visit" && n.kind=="survival") fail("Repeated event processes require a separate event-history model. Use a subject-level survival variable.");
    const auto &r=v.get("repeated");
    if(!r.is_null() && !r.is_object()) fail("Repeated-measure settings must be an object.");
    n.subjectSD=optional(r,"subjectSD"); n.slopeSD=optional(r,"slopeSD");
    n.slopeCorrelation=optional(r,"slopeCorrelation"); n.clusterSD=optional(r,"clusterSD");
    n.residualCorrelation=optional(r,"residualCorrelation");
    if(n.subjectSD<0 || n.slopeSD<0 || n.clusterSD<0 || std::fabs(n.slopeCorrelation)>1 || n.residualCorrelation<0 || n.residualCorrelation>=1)
      fail("Random-effect SDs must be nonnegative; slope correlation in [-1,1]; residual correlation in [0,1).");
    if((n.subjectSD || n.slopeSD || n.clusterSD || n.residualCorrelation) && n.scope!="visit") fail("Repeated random effects require a visit-level variable.");
    if(n.clusterSD && s.clusterSizes.empty()) fail("A cluster random effect requires cluster sizes.");
    if(n.residualCorrelation && n.family!="gaussian" && n.family!="lognormal") fail("Residual correlation is available for Normal and log-normal repeated errors.");
    if(v.find("lags")) {
      if(!v.get("lags").is_array() || v.get("lags").size()>40) fail("Use at most 40 lagged relationships per variable.");
      std::set<std::pair<int,int>> used;
      for(const auto &l:v.get("lags").elements()) {
        Node::Lag lag; lag.source=g.index(l.get("source").as_string());
        lag.lag=integer(l,"lag",1,12,1); lag.coefficient=number(l,"coefficient"); lag.initial=optional(l,"initial");
        if(n.scope!="visit" || !used.insert({lag.source,lag.lag}).second) fail("Lagged relationships need a visit-level child and distinct source/lag pairs.");
        n.lags.push_back(lag);
      }
    }
    if(n.scope=="visit") for(const char *suffix:{"_observed","_recorded"})
      if(!outputs.insert(n.id+suffix).second) fail("A repeated observation column conflicts with a variable ID.");
  }
  // Scopes are now known for every parent (node order need not be topological).
  for(const auto &n:g.nodes) for(const auto &t:n.terms) for(int p:{t.left,t.right})
    if(p>=0 && level(g.nodes[p])>level(n)) fail("A variable cannot depend on a finer observation level.");
  const auto &o=raw.get("observation");
  if(!o.is_null() && !o.is_object()) fail("Observation settings must be an object.");
  s.missing=o.find("mode") ? o.get("mode").as_string() : "none";
  if(s.missing!="none" && s.missing!="dropout" && s.missing!="intermittent") fail("Choose no missing visits, dropout or intermittent missing visits.");
  s.missingIntercept=optional(o,"intercept",-3);
  s.historyCoefficient=optional(o,"historyCoefficient"); s.currentCoefficient=optional(o,"currentCoefficient");
  s.treatmentCoefficient=optional(o,"treatmentCoefficient");
  if(o.find("historyVariable") && !o.get("historyVariable").as_string().empty()) {
    s.history=g.index(o.get("historyVariable").as_string());
    if(g.nodes[s.history].scope!="visit") fail("Missing visits must reference a visit-level history variable.");
  }
  if(o.find("treatmentVariable") && !o.get("treatmentVariable").as_string().empty()) {
    s.treatment=g.index(o.get("treatmentVariable").as_string());
    if(level(g.nodes[s.treatment])>1) fail("Missing-visit treatment effects require a subject- or cluster-level variable.");
  }
  if((s.historyCoefficient || s.currentCoefficient) && s.history<0) fail("Choose the variable determining missing visits.");
  if(s.treatmentCoefficient && s.treatment<0) fail("Choose the treatment variable determining missing visits.");
}
int study_rows(const Graph &g, int units) {
  long long subjects=0;
  for(int c=0;c<units;++c) subjects+=cluster_size(g,c);
  const long long rows=subjects*g.study.times.size();
  if(rows>50000) fail("This study exceeds 50,000 measurement rows. Reduce subjects, clusters, cluster sizes or visits.");
  return static_cast<int>(rows);
}
Data simulate_study(const Graph &g, int units, std::uint64_t seed,
                    const std::map<int,double> &interventions, bool observe, const Token *token) {
  if(g.nodes.empty()) fail("Add at least one variable before simulating data.");
  if(units<1 || units>50000) fail("Study size must be an integer from 1 to 50,000.");
  const int visits=static_cast<int>(g.study.times.size());
  Data d; d.n=study_rows(g,units); d.units=units; d.subjects=d.n/visits;
  const int p=static_cast<int>(g.nodes.size());
  d.values.assign(p,Vec(d.n)); d.times.resize(p); d.events.resize(p);
  d.observed.assign(p,Vec(d.n,1));
  d.subject.resize(d.n); d.cluster.resize(d.n); d.visit.resize(d.n);
  for(const auto &kv:interventions) {
    const auto &v=g.nodes.at(kv.first); const double x=kv.second;
    if(!std::isfinite(x) || v.scope=="time" || (v.kind=="binary" && x!=0 && x!=1) ||
       (v.kind=="count" && (x<0 || x!=std::floor(x))) || (v.kind=="survival" && x<=0)) fail("Invalid intervention value or scope.");
  }
  int subject=0;
  for(int c=0;c<units;++c) {
    const int firstCluster=subject*visits;
    for(int member=0;member<cluster_size(g,c);++member,++subject) {
      const int start=subject*visits;
      Vec innovation(p); bool dropped=false; int lastObserved=0;
      for(int v=0;v<visits;++v) {
        const int r=start+v;
        if(token) token->check();
        d.subject[r]=subject; d.cluster[r]=c; d.visit[r]=v;
        for(int j:g.order) {
          const auto &node=g.nodes[j]; double value=0;
          if(node.scope=="time") value=g.study.times[v];
          else if(node.scope=="cluster" && r!=firstCluster) value=d.values[j][firstCluster];
          else if(node.scope=="subject" && v>0) value=d.values[j][start];
          else if(interventions.count(j)) value=interventions.at(j);
          else {
            double eta=node.intercept;
            for(const auto &t:node.terms) eta+=t.coefficient*term_value(t,d.values[t.left][r],t.right<0?1:d.values[t.right][r]);
            for(const auto &l:node.lags) eta+=l.coefficient*(v>=l.lag ? d.values[l.source][r-l.lag] : l.initial);
            const std::uint64_t identity=node.scope=="cluster" ? c : subject;
            sim::Rng rng(seed,node.stream,sim::stream_seed(identity,static_cast<std::uint64_t>(v),0));
            double z=NAN;
            if(node.scope=="visit") {
              sim::Rng person(seed,node.stream ^ 0x6375727665ULL,subject), group(seed,node.stream ^ 0x67726f7570ULL,c);
              const double b0=person.normal(), b1=node.slopeCorrelation*b0+std::sqrt(1-node.slopeCorrelation*node.slopeCorrelation)*person.normal();
              eta+=node.subjectSD*b0+node.slopeSD*b1*(g.study.times[v]-g.study.times.front())+node.clusterSD*group.normal();
              if(node.family=="gaussian" || node.family=="lognormal") {
                const double rho=v && node.residualCorrelation ? std::pow(node.residualCorrelation,g.study.times[v]-g.study.times[v-1]) : 0;
                innovation[j]=rho*innovation[j]+std::sqrt(1-rho*rho)*rng.normal(); z=innovation[j];
              }
            }
            value=draw_value(node,eta,rng,z);
          }
          d.values[j][r]=value;
        }
        if(observe && v>0 && g.study.missing!="none") {
          const auto &s=g.study;
          double eta=s.missingIntercept;
          if(s.history>=0) eta+=s.historyCoefficient*d.values[s.history][start+lastObserved]+s.currentCoefficient*d.values[s.history][r];
          if(s.treatment>=0) eta+=s.treatmentCoefficient*d.values[s.treatment][r];
          if(!std::isfinite(eta)) throw Error("generation_error","The missing-visit predictor overflowed.");
          sim::Rng rng(seed,0x6d697373696e67ULL,sim::stream_seed(subject,v,0));
          const bool missing=(s.missing=="dropout" && dropped) || rng.uniform()<logistic(eta);
          if(missing) { dropped=true; for(int j=0;j<p;++j) if(g.nodes[j].scope=="visit") d.observed[j][r]=0; }
          else lastObserved=v;
        }
      }
    }
  }
  if(observe) for(int j=0;j<p;++j) if(g.nodes[j].kind=="survival") {
    const auto &node=g.nodes[j]; d.times[j].resize(d.n); d.events[j].resize(d.n);
    for(int r=0;r<d.n;++r) {
      const double c=node.censor=="none" ? INFINITY : node.censor=="administrative" ? node.cutoff : d.values[node.censor_index][r];
      if(!(c>=0)) throw Error("generation_error","Censoring times must be nonnegative.");
      d.times[j][r]=std::min(d.values[j][r],c); d.events[j][r]=d.values[j][r]<=c;
    }
  }
  return d;
}
void read_study_data(const Graph &g, const J &raw, Data &d) {
  const int visits=static_cast<int>(g.study.times.size());
  if(d.n%visits) fail("A study dataset must retain a row for every scheduled visit, with observation indicators.");
  d.subjects=d.n/visits; d.observed.assign(g.nodes.size(),Vec(d.n,1));
  auto ids=[&](const char *name,std::vector<int> &out) {
    const auto &a=raw.get("columns").get(name);
    if(!a.is_array() || a.size()!=static_cast<std::size_t>(d.n)) fail(std::string("Missing study identifiers: ")+name);
    for(const auto &v:a.elements()) { const double x=v.as_number(-1); if(x<0 || x!=std::floor(x) || x>d.n) fail("Invalid study identifier."); out.push_back(static_cast<int>(x)); }
  };
  ids("subject_id",d.subject); ids("cluster_id",d.cluster); ids("visit_index",d.visit);
  int subject=0,cluster=0,left=cluster_size(g,0);
  for(int r=0;r<d.n;++r) {
    if(d.subject[r]!=subject || d.cluster[r]!=cluster || d.visit[r]!=r%visits) fail("Study rows must follow their saved cluster, subject and visit structure.");
    if(r%visits==visits-1) { ++subject; if(--left==0) { ++cluster; left=cluster_size(g,cluster); } }
  }
  if(left!=cluster_size(g,cluster)) fail("A dataset contains an incomplete cluster.");
  d.units=cluster;
  for(int j=0;j<static_cast<int>(g.nodes.size());++j) if(g.nodes[j].scope=="visit") {
    const auto &a=raw.get("columns").get(g.nodes[j].id+"_observed");
    if(!a.is_array() || a.size()!=static_cast<std::size_t>(d.n)) fail("Repeated datasets need observation indicators.");
    for(int r=0;r<d.n;++r) { const double x=a.at(r).as_number(-1); if(x!=0 && x!=1) fail("Observation indicators must be 0 or 1."); d.observed[j][r]=x; }
  }
  for(int j=0;j<static_cast<int>(g.nodes.size());++j) if(g.nodes[j].scope=="subject" || g.nodes[j].scope=="cluster") {
    int first=0;
    for(int r=0;r<d.n;++r) {
      if(!r || (g.nodes[j].scope=="subject" ? d.subject[r]!=d.subject[r-1] : d.cluster[r]!=d.cluster[r-1])) first=r;
      if(d.values[j][r]!=d.values[j][first]) fail("A subject- or cluster-level variable changed within its unit.");
    }
  }
}
Vec study_weights(const Graph &g, const J &e) {
  const int visits=static_cast<int>(g.study.times.size());
  const auto &r=e.get("repeated");
  const std::string kind=r.find("contrast") ? r.get("contrast").as_string() : "final";
  Vec w(visits);
  if(kind=="final") w[integer(r,"visit",0,visits-1,visits-1)]=1;
  else if(kind=="average") std::fill(w.begin(),w.end(),1.0/visits);
  else if(kind=="change") { if(visits<2) fail("Change needs at least two visits."); w.front()=-1; w.back()=1; }
  else if(kind=="slope") {
    if(visits<2) fail("A trend needs at least two visits.");
    const auto m=moments(g.study.times); double sum=0;
    for(int v=0;v<visits;++v) { w[v]=g.study.times[v]-m.mean; sum+=w[v]*w[v]; }
    for(double &v:w) v/=sum;
  } else if(kind=="custom") {
    if(!r.get("weights").is_array() || r.get("weights").size()!=static_cast<std::size_t>(visits)) fail("Custom contrast weights must match the visit schedule.");
    for(int v=0;v<visits;++v) { const auto &x=r.get("weights").at(v); if(!x.is_number() || !std::isfinite(x.as_number())) fail("Contrast weights must be finite."); w[v]=x.as_number(); }
  } else fail("Choose final visit, average, change, linear trend or custom contrast.");
  if(std::all_of(w.begin(),w.end(),[](double x){return x==0;})) fail("At least one contrast weight must be nonzero.");
  return w;
}
J study_plan(const Graph &g, const J &e, const J &spec) {
  const int a=g.index(e.get("exposure").as_string()), y=g.index(e.get("outcome").as_string());
  const auto &exposure=g.nodes[a], &outcome=g.nodes[y];
  if(e.get("kind").as_string()!="marginal" || (e.get("measure").as_string()!="mean_difference" && e.get("measure").as_string()!="risk_difference"))
    fail("Structured studies currently estimate marginal mean/risk contrasts across visits.");
  if(exposure.kind!="binary" || level(exposure)>1)
    fail("Structured power requires a binary baseline subject- or cluster-level exposure. Time-varying treatment needs a separate longitudinal estimator.");
  if(outcome.kind=="survival" || outcome.scope=="time") fail("Select a nonsurvival measured outcome for a repeated contrast.");
  if(e.get("active").as_number()==e.get("reference").as_number()) fail("Comparison levels must differ.");
  const auto weights=study_weights(g,e);
  std::string method=spec.get("method").as_string();
  const bool grouped=!g.study.clusterSizes.empty(), continuous=outcome.kind=="continuous";
  const int visits=static_cast<int>(g.study.times.size());
  if(method.empty() || method=="auto") method=!grouped && continuous && outcome.scope=="visit" && visits>1 && g.study.missing!="intermittent" ? "mmrm" : "trajectory";
  if(method!="tmle" && !exposure.terms.empty()) fail("For baseline confounding, select TMLE with complete outcomes. Other structured analyses require randomized exposure.");
  if(method=="tmle" && g.study.missing!="none") fail("Trajectory TMLE requires complete outcomes. Dropout and intermittent missingness need a longitudinal missing-data estimator; this TMLE does not provide one.");
  if(method!="mmrm" && method!="lmm" && method!="trajectory" && method!="tmle") fail("Choose MMRM, linear mixed model or trajectory means for a structured study.");
  if((method=="mmrm" || method=="lmm") && (grouped || !continuous || outcome.scope!="visit")) fail("Mixed-model trajectory analyses currently require continuous outcomes and independent subjects. Use trajectory means for clustered or discrete outcomes.");
  if(method=="mmrm" && g.study.missing=="intermittent") fail("This MMRM fitter supports monotone dropout. Use LMM or trajectory means for intermittent missing visits.");
  const std::string timeBasis=spec.find("timeBasis") ? spec.get("timeBasis").as_string() : "categorical";
  if(timeBasis!="categorical" && timeBasis!="linear") fail("Choose categorical or linear time in the planned LMM.");
  if(spec.find("randomSlope") && !spec.get("randomSlope").is_bool()) fail("The random-slope analysis setting must be true or false.");
  if(method=="lmm" && (visits<2 || (timeBasis=="categorical" && 2*visits>sim::lmm::PMAX))) fail("The LMM supports 2–6 categorical visits or 2–12 visits with linear time.");
  if(method!="tmle" && integer(spec,"folds",1,10,1)!=1) fail("These structured estimators do not use cross-fitting.");
  J notes=J::array();
  notes.push("Power is conditional on the generating distributions, response functions, dependence and observation process. Monte Carlo intervals measure simulation error, not uncertainty in those assumptions.");
  notes.push(grouped ? "Sample size counts independent clusters. Effects average over subjects; unequal cluster sizes follow the saved repeating size pattern." : "Sample size counts independent subjects; visits are repeated observations of those subjects.");
  if(method=="tmle") notes.push("TMLE targets the subject-weighted mean of the complete outcome summary under baseline treatment. It adjusts only for baseline covariates. Folds and influence-function inference preserve independent subjects or entire clusters. It is not longitudinal TMLE for time-varying treatment or dropout.");
  if(method=="mmrm") notes.push("MMRM fits separate arm-by-visit means and a common unstructured covariance. The t reference uses last-visit sample size minus two, not Kenward–Roger degrees of freedom.");
  if(method=="lmm") notes.push("LMM estimates its random-intercept/slope covariance and residual variance by REML, with Satterthwaite inference. Its covariance assumptions are independent of the generating model.");
  if(method=="trajectory") notes.push("Trajectory means use arm-by-visit averages and a CR1 sandwich clustered at the highest independent unit, with units minus two degrees of freedom. This is an approximation for few independent units.");
  if(g.study.missing!="none") notes.push(method=="trajectory" ? "Available-visit means require noninformative missingness within arm and visit. History-dependent dropout can bias this analysis; inspect simulated bias and null rejection." : "Likelihood analyses rely on the planned mean/covariance model and ignorable missingness. Selection on an unobserved current outcome can bias them.");
  if(method!="tmle" && (spec.get("basis").as_string()!="match" || spec.find("propensity") || spec.find("models"))) fail("Structured analyses use their trajectory settings, not cross-sectional nuisance-model terms.");
  return object({{"method",method},{"label",method=="mmrm" ? "MMRM contrast of arm-by-visit means" : method=="lmm" ? "Linear mixed-model trajectory contrast" : "Trajectory means with cluster-robust inference"},
    {"formula",outcome.id+" ~ arm × "+(method=="lmm" && timeBasis=="linear" ? "time" : "visit")},
    {"notes",notes},{"weights",array(weights)},{"timeBasis",timeBasis}});
}
Fit analyze_study(const Graph &g, const J &e, const Data &d, const J &spec) {
  const auto plan=study_plan(g,e,spec); const auto method=plan.get("method").as_string();
  const int a=g.index(e.get("exposure").as_string()), y=g.index(e.get("outcome").as_string());
  const int visits=static_cast<int>(g.study.times.size()); const Vec w=study_weights(g,e);
  auto arm=[&](int r) { const double x=d.values[a][r]; if(x!=e.get("active").as_number() && x!=e.get("reference").as_number()) throw Error("invalid_data","Unexpected treatment level."); return x==e.get("active").as_number()?1:0; };
  Fit out;
  if(method=="mmrm") {
    procs::mmrm::Stats s; s.visits=visits;
    for(int i=0;i<d.subjects;++i) {
      const int start=i*visits; int k=0;
      while(k<visits && d.observed[y][start+k]) ++k;
      for(int v=k;v<visits;++v) if(d.observed[y][start+v]) throw Error("invalid_data","MMRM requires monotone observation patterns.");
      if(k) s.add(arm(start),k,d.values[y].data()+start);
    }
    const auto f=procs::mmrm::fit(s,true,w.data());
    if(!f.ok) throw Error("nonconvergence","The repeated-measures covariance fit failed.");
    // The existing fitter returns group 0 minus group 1.
    out.coefficient=-f.delta; out.se=f.se; out.df=f.observed_last-2;
  } else if(method=="lmm") {
    sim::lmm::Data s; const bool linear=plan.get("timeBasis").as_string()=="linear";
    s.p=linear ? 4 : 2*visits; s.q=spec.get("randomSlope").as_bool(false) ? 2 : 1;
    for(int i=0;i<d.subjects;++i) {
      Vec x,z,yv;
      for(int v=0;v<visits;++v) {
        const int r=i*visits+v; if(!d.observed[y][r]) continue;
        const int group=arm(r); const double time=g.study.times[v]-g.study.times.front();
        if(linear) { x.insert(x.end(),{1.0,static_cast<double>(group),time,group*time}); }
        else { const auto start=x.size(); x.resize(start+s.p); x[start+v]=1; x[start+visits+v]=group; }
        z.push_back(1); if(s.q==2) z.push_back(time); yv.push_back(d.values[y][r]);
      }
      if(yv.empty()) continue;
      sim::lmm::Pattern p; p.set_design(static_cast<long>(yv.size()),s.p,s.q,x.data(),z.data()); p.add_cluster(s.p,s.q,x.data(),z.data(),yv.data()); s.patterns.push_back(p);
    }
    const auto f=sim::lmm::fit(s); if(!f.ok) throw Error("nonconvergence","The linear mixed model did not converge.");
    Vec c(s.p);
    for(int v=0;v<visits;++v) if(linear) { c[1]+=w[v]; c[3]+=w[v]*(g.study.times[v]-g.study.times.front()); } else c[visits+v]=w[v];
    const auto test=sim::lmm::wald(s,f,c.data(),0,s.clusters()-2);
    out.coefficient=test.estimate; out.se=test.se; out.df=test.df;
  } else {
    Vec sum(2*visits), count(2*visits), mean(2*visits), influence(d.units);
    for(int r=0;r<d.n;++r) if(d.observed[y][r]) { const int j=arm(r)*visits+d.visit[r]; sum[j]+=d.values[y][r]; ++count[j]; }
    for(int j=0;j<2*visits;++j) {
      if(count[j]<2 && w[j%visits]!=0) throw Error("empty_arm","An arm has fewer than two observed outcomes at a requested visit.");
      mean[j]=count[j] ? sum[j]/count[j] : 0;
    }
    for(int v=0;v<visits;++v) out.coefficient+=w[v]*(mean[visits+v]-mean[v]);
    for(int r=0;r<d.n;++r) if(d.observed[y][r]) {
      const int group=arm(r),v=d.visit[r],j=group*visits+v;
      influence[d.cluster[r]]+=(group?1:-1)*w[v]*(d.values[y][r]-mean[j])/count[j];
    }
    if(d.units<=2) throw Error("few_units","Cluster-robust inference needs at least three independent units.");
    double variance=0; for(double v:influence) variance+=v*v;
    out.se=std::sqrt(variance*d.units/(d.units-1.0)); out.df=d.units-2;
  }
  if(!std::isfinite(out.coefficient) || !std::isfinite(out.se) || out.se<=0 || !(out.df>0)) throw Error("zero_variance","The trajectory contrast has no finite positive standard error or degrees of freedom.");
  return out;
}
J study_population(const Graph &g, const J &e, int units, std::uint64_t seed, const Token *token) {
  const Vec w=study_weights(g,e); const int visits=static_cast<int>(w.size());
  const int a=g.index(e.get("exposure").as_string()), y=g.index(e.get("outcome").as_string());
  if(e.get("kind").as_string()!="marginal" || (e.get("measure").as_string()!="mean_difference" && e.get("measure").as_string()!="risk_difference")) fail("Structured population targets require marginal mean/risk contrasts.");
  Vec totals, sizes;
  const int maxSize=g.study.clusterSizes.empty()?1:*std::max_element(g.study.clusterSizes.begin(),g.study.clusterSizes.end());
  const int batch=std::max(1,50000/(visits*maxSize));
  for(int begin=0;begin<units;begin+=batch) {
    Graph local=g;
    if(!local.study.clusterSizes.empty()) std::rotate(local.study.clusterSizes.begin(),local.study.clusterSizes.begin()+begin%local.study.clusterSizes.size(),local.study.clusterSizes.end());
    const int n=std::min(batch,units-begin); const auto stream=sim::stream_seed(seed,begin,0x706f7075ULL);
    const auto lo=simulate_study(local,n,stream,{{a,e.get("reference").as_number()}},false,token), hi=simulate_study(local,n,stream,{{a,e.get("active").as_number()}},false,token);
    Vec sums(n), counts(n);
    for(int r=0;r<lo.n;++r) { sums[lo.cluster[r]]+=w[lo.visit[r]]*(hi.values[y][r]-lo.values[y][r]); if(lo.visit[r]==0) ++counts[lo.cluster[r]]; }
    totals.insert(totals.end(),sums.begin(),sums.end()); sizes.insert(sizes.end(),counts.begin(),counts.end());
  }
  const double size=moments(sizes).mean, estimate=moments(totals).mean/size;
  Vec influence(units); for(int k=0;k<units;++k) influence[k]=(totals[k]-estimate*sizes[k])/size;
  const double se=std::sqrt(moments(influence).variance/units);
  if(!std::isfinite(estimate) || !std::isfinite(se)) fail("The structured population effect overflowed.");
  return object({{"estimate",estimate},{"mcSE",se},{"interval",array({estimate-1.96*se,estimate+1.96*se})},{"n",units},{"exact",false},{"weights",array(w)}});
}
}
