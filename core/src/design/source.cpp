// Compiler for the documented, closed subset of Lean study declarations.
// It does not run Lean tactics or claim that parsing is a kernel proof check.
#include "plan.hpp"
#include "util/numtext.hpp"
#include <cctype>
#include <cmath>
#include <cstdint>
#include <map>
#include <numeric>
#include <set>
#include <sstream>

namespace ssc::design {
namespace {
using J=json::Value;
J obj(std::initializer_list<std::pair<std::string,J>> xs) {J j=J::object();for(const auto& x:xs)j.set(x.first,x.second);return j;}
struct TokenText {std::string text;std::size_t position=0;bool quoted=false;};
class Parser {
  const std::string& source;
  std::size_t position=0;
 public:
  TokenText token;
  explicit Parser(const std::string& s):source(s) {next();}
  [[noreturn]] void error(const std::string& message) const {
    int line=1,column=1;
    for(std::size_t i=0;i<token.position && i<source.size();++i) {if(source[i]=='\n'){++line;column=1;}else ++column;}
    throw std::invalid_argument("Lean study, line "+std::to_string(line)+", column "+std::to_string(column)+": "+message);
  }
  void next() {
    while(position<source.size()) {
      if(std::isspace(static_cast<unsigned char>(source[position]))) {++position;continue;}
      if(source.compare(position,2,"--")==0) {while(position<source.size()&&source[position]!='\n')++position;continue;}
      if(source.compare(position,2,"/-")==0) {
        position+=2;int depth=1;
        while(position<source.size()&&depth) {
          if(source.compare(position,2,"/-")==0) {++depth;position+=2;}
          else if(source.compare(position,2,"-/")==0) {--depth;position+=2;} else ++position;
        }
        if(depth) {token.position=position;error("unterminated block comment.");} continue;
      }
      break;
    }
    token={"",position,false};if(position==source.size())return;
    const char c=source[position];
    if(c=='"') {
      const auto begin=position++;
      bool closed=false;
      while(position<source.size()) {
        if(source[position]=='"'){++position;closed=true;break;}
        if(source[position]=='\\') {
          ++position;
          if(position==source.size()||std::string("\"\\nrt").find(source[position])==std::string::npos)
            error("use only quote, backslash, n, r, or t string escapes.");
        }
        ++position;
      }
      if(!closed)error("unterminated string.");
      token.text=json::parse(source.substr(begin,position-begin)).as_string();token.quoted=true;return;
    }
    if(std::isalpha(static_cast<unsigned char>(c))||c=='_') {
      const auto begin=position++;
      while(position<source.size()&&(std::isalnum(static_cast<unsigned char>(source[position]))||source[position]=='_'||source[position]=='.'))++position;
      token.text=source.substr(begin,position-begin);return;
    }
    if(std::isdigit(static_cast<unsigned char>(c))) {
      const auto begin=position++;
      while(position<source.size()&&(std::isdigit(static_cast<unsigned char>(source[position]))||source[position]=='.'))++position;
      token.text=source.substr(begin,position-begin);return;
    }
    if(source.compare(position,2,":=")==0){token.text=":=";position+=2;return;}
    if(std::string("{}:,.()/+-").find(c)!=std::string::npos) {token.text=source.substr(position++,1);return;}
    error("unsupported token. The browser accepts a single NormalInput or InfluenceInput declaration.");
  }
  bool take(const std::string& x) {if(!token.quoted&&token.text==x){next();return true;}return false;}
  void expect(const std::string& x) {if(!take(x))error("expected '"+x+"'.");}
  std::pair<std::int64_t,std::int64_t> decimal() {
    const auto t=token.text;
    if(token.quoted||t.empty()||!std::isdigit(static_cast<unsigned char>(t[0])))error("expected a numeric literal.");
    const auto dot=t.find('.');
    if(t.find('.',dot==std::string::npos?t.size():dot+1)!=std::string::npos||t.back()=='.')error("invalid decimal literal.");
    std::string digits=t; if(dot!=std::string::npos)digits.erase(dot,1);
    if(digits.size()>15)error("numeric literals may have at most 15 digits.");
    std::int64_t numerator=std::stoll(digits),denominator=1;
    if(dot!=std::string::npos)for(std::size_t k=dot+1;k<t.size();++k)denominator*=10;
    next();return {numerator,denominator};
  }
  J number() {
    const bool paren=take("(");const bool negative=take("-");
    auto value=decimal();
    if(take("/")) {
      const auto divisor=decimal();
      if(divisor.second!=1||divisor.first<=0||value.second!=1)error("fractions must have an integer numerator and positive integer denominator.");
      value.second=divisor.first;
    }
    if(paren)expect(")");
    if(negative)value.first=-value.first;
    const auto common=std::gcd(value.first,value.second);value.first/=common;value.second/=common;
    return std::to_string(value.first)+"/"+std::to_string(value.second);
  }
  J parse() {
    expect("import");expect("DesignCompiler.Input");expect("open");expect("DesignCompiler");
    expect("def");expect("study");expect(":");
    const auto type=token.text; if(token.quoted||(type!="NormalInput"&&type!="InfluenceInput"))error("use NormalInput or InfluenceInput.");next();expect(":=");expect("{");
    const bool normal=type=="NormalInput";
    J x=obj({{"alternative","1/2"},{"nullValue","0/1"},{"alpha","1/20"},{"targetPower","4/5"}});
    if(normal)x.set("sd","1/1");
    else {x.set("variance","4/1");x.set("independentUnit","participant");x.set("estimand","ATE");x.set("estimator","AIPW");x.set("varianceSource","Supplied per-unit planning variance");}
    J search=obj({{"minN",1},{"maxN",1000000}});
    std::string name=normal?"Normal study":"Influence-function study";
    std::set<std::string> seen;
    while(!take("}")) {
      if(token.text.empty())error("expected a field or '}'.");
      const auto key=token.text;
      if(token.quoted||!seen.insert(key).second)error("invalid or duplicate field '"+key+"'.");
      next();expect(":=");
      if(key=="name"||(!normal&&(key=="estimand"||key=="estimator"||key=="varianceSource"))) {
        if(!token.quoted)error("expected a quoted string.");
        if(key=="name")name=token.text;else x.set(key,token.text);next();
      } else if(key=="unit"&&!normal) {
        expect(".");const auto unit=token.text;
        if(token.quoted||(unit!="participant"&&unit!="cluster"))error("unit must be .participant or .cluster.");
        x.set("independentUnit",unit);next();
      } else if(key=="minN"||key=="maxN") {
        const auto n=decimal();if(n.second!=1||n.first>10000000)error("sample-size limits must be natural numbers at most 10000000.");
        search.set(key,static_cast<long long>(n.first));
      } else if(x.find(key)) x.set(key,number());
      else error("unsupported field '"+key+"' for "+type+".");
      // Lean structure fields are separated by a comma or a new line.
      if(!take(",")&&token.text!="}") {
        auto previous=token.position;
        while(previous>0&&std::isspace(static_cast<unsigned char>(source[previous-1])))--previous;
        if(source.substr(previous,token.position-previous).find('\n')==std::string::npos)
          error("separate fields with a comma or new line.");
      }
    }
    if(!token.text.empty())error("unexpected code after the study. Only the documented declaration subset is accepted in the browser.");
    return obj({{"schema","ssc.design-plan"},{"version",1},{"method",normal?"normal_known":"influence_wald"},
      {"name",name},{"parameters",x},{"search",search}});
  }
};

std::string cpp_source(const J& plan,const J& checked) {
  const auto number=[](double x){return text::format_g(x,17);};
  const auto& solve=plan.get("solve");const auto mode=solve.get("for").as_string();
  const bool perArm=plan.get("method").as_string()=="normal_known";
  std::ostringstream code;
  code << "// Generated from the supported Lean study declaration. C++17, standard library only.\n"
          "// Build: c++ -std=c++17 -O2 study.cpp -o study\n"
          "// Ordinary floating-point evaluation; no certified numerical error bound.\n"
       << (perArm ? "// Normal outcomes, independent participants, equal arms, common known SD. n is PER ARM.\n"
                  : "// Asymptotic influence-function approximation. Variance must match the planned estimator.\n")
       << "// Count unit: " << checked.get("sizeUnit").as_string() << ". Variance is held fixed when solving for effect.\n"
       <<
          "#include <cmath>\n#include <iomanip>\n#include <iostream>\n\n"
       << "constexpr double variance = " << number(checked.get("varianceCoefficient").as_number()) << ";\n"
       << "constexpr double alternative = " << number(checked.get("alternative").as_number()) << ";\n"
       << "constexpr double null_value = " << number(checked.get("nullValue").as_number()) << ";\n"
       << "constexpr double alpha = " << number(checked.get("alpha").as_number()) << ";\n"
       << "constexpr double target = " << number(checked.get("targetPower").as_number()) << ";\n"
       << "double upper_tail(double x) { return 0.5 * std::erfc(x / std::sqrt(2.0)); }\n"
          "double critical_value() {\n  double lo = 0, hi = 40;\n  for (int i=0;i<100;++i) { double m=(lo+hi)/2; if(upper_tail(m)>alpha/2) lo=m; else hi=m; }\n  return (lo+hi)/2;\n}\n"
          "double power(double effect, int n) {\n  static const double z=critical_value();\n  const double lambda=std::abs(effect)/std::sqrt(variance)*std::sqrt(double(n));\n  return upper_tail(z+lambda)+upper_tail(z-lambda);\n}\n"
          "int main() {\n  std::cout << std::setprecision(17);\n";
  if(mode=="sample_size") {
    code << "  for(int n="<<static_cast<int>(plan.get("search").get("minN").as_number())<<"; n<="<<static_cast<int>(plan.get("search").get("maxN").as_number())<<"; ++n) {\n"
      "    const double p=power(alternative-null_value,n);\n"
      "    if(std::abs(p-target)<=1e-12) { std::cout << \"{\\\"status\\\":\\\"numerically_unresolved\\\"}\\n\"; return 0; }\n"
      "    if(p>target) { std::cout << \"{\\\"selectedN\\\":\" << n << \"";
    if(perArm)code<<",\\\"totalParticipants\\\":"<<"\" << 2*n << \"";
    code << ",\\\"power\\\":\" << p << \"}\\n\"; return 0; }\n  }\n"
      "  std::cout << \"{\\\"status\\\":\\\"not_found_in_domain\\\"}\\n\";\n";
  } else if(mode=="power") {
    code << "  const int n="<<static_cast<int>(solve.get("n").as_number())<<";\n"
      "  std::cout << \"{\\\"power\\\":\" << power(alternative-null_value,n) << \"}\\n\";\n";
  } else {
    code << "  const int n="<<static_cast<int>(solve.get("n").as_number())<<";\n"
      "  double lo=0,hi=64;\n  const double scale=std::sqrt(variance)/std::sqrt(double(n));\n"
      "  for(int i=0;i<100;++i) {double m=(lo+hi)/2; if(m==lo||m==hi)break; if(power(m*scale,n)>=target)hi=m;else lo=m;}\n"
      "  std::cout << \"{\\\"detectableMagnitude\\\":\" << hi*scale << \"}\\n\";\n";
  }
  code<<"}\n";return code.str();
}
} // namespace

J compile_source(const J& request) {
  if(!request.is_object())throw std::invalid_argument("Expected a compilation request.");
  for(const auto& kv:request.items())if(kv.first!="op"&&kv.first!="source"&&kv.first!="solveFor"&&kv.first!="n")
    throw std::invalid_argument("Unsupported compilation field: "+kv.first);
  const auto& source=request.get("source");
  if(!source.is_string()||source.as_string().empty()||source.as_string().size()>65536)
    throw std::invalid_argument("Supply a Lean study declaration of 1–65536 bytes.");
  J plan=Parser(source.as_string()).parse();
  const auto mode=request.find("solveFor")?request.get("solveFor").as_string():"sample_size";
  if(mode!="sample_size"&&mode!="power"&&mode!="effect")throw std::invalid_argument("Choose sample_size, power, or effect.");
  J solve=obj({{"for",mode}});
  if(mode!="sample_size") {
    const auto& n=request.get("n");double value=n.as_number(-1);
    if(!n.is_number()||!std::isfinite(value)||value<1||value>10000000||value!=std::floor(value))
      throw std::invalid_argument("n must be a whole number from 1 to 10000000.");
    solve.set("n",n);
  } else if(request.find("n"))throw std::invalid_argument("Do not supply n when solving for sample size.");
  plan.set("solve",solve);
  // Validate without performing a sample-size search. The probe also supplies the
  // normalized variance coefficient for independent, downloadable C++ source.
  J probe=plan;probe.set("solve",obj({{"for","power"},{"n",1}}));
  const auto checked=execute(probe);
  return obj({{"plan",plan},{"cpp",cpp_source(plan,checked)},
    {"compiler","ssc-lean-study-subset-v1"},{"leanKernelChecked",false},
    {"language","Lean 4 structure declaration subset; native Lean checking is available through scripts/compile_lean.py"},
    {"execution","The browser executes this plan in the compiled C++/WebAssembly engine. Downloaded C++ can be compiled separately."}});
}
} // namespace ssc::design
