// ssc_cli — command-line access to the computational core.
//
//   ssc_cli list                      list procedures
//   ssc_cli schema <procedure-id>     print a procedure's JSON schema
//   ssc_cli run <request.json | ->    run a request (file or stdin), print JSON
//   ssc_cli eval <fn> <args...>       evaluate a distribution function
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "engine/procedure.hpp"
#include "ssc/ssc.h"
#include "util/numtext.hpp"

namespace {
int usage() {
  std::cerr << "usage:\n"
               "  ssc_cli list\n"
               "  ssc_cli schema <procedure-id>\n"
               "  ssc_cli search <words...>\n"
               "  ssc_cli run <request.json | ->\n"
               "  ssc_cli eval <function> <args...>\n";
  return 2;
}

// NULL: not even memory for the answer (see ssc.h)
int print_and_free(char* s) {
  if (!s) {
    std::cerr << "there is not enough memory for the answer\n";
    return 1;
  }
  std::cout << s << "\n";
  ssc_free(s);
  return 0;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) return usage();
  const std::string cmd = argv[1];
  if (cmd == "search") {
    std::string query;
    for (int i = 2; i < argc; ++i) query += (i > 2 ? " " : "") + std::string(argv[i]);
    return print_and_free(ssc_search_json(query.c_str(), 20));
  }
  if (cmd == "list") return print_and_free(ssc_catalog_json());
  if (cmd == "schema" && argc >= 3) {
    char* s = ssc_procedure_json(argv[2]);
    if (!s && !ssc::find_procedure(argv[2])) {
      std::cerr << "unknown procedure: " << argv[2] << "\n";
      return 1;
    }
    return print_and_free(s);
  }
  if (cmd == "run" && argc >= 3) {
    std::string text;
    if (std::string(argv[2]) == "-") {
      std::stringstream ss;
      ss << std::cin.rdbuf();
      text = ss.str();
    } else {
      std::ifstream in(argv[2]);
      if (!in) {
        std::cerr << "cannot open " << argv[2] << "\n";
        return 1;
      }
      std::stringstream ss;
      ss << in.rdbuf();
      text = ss.str();
    }
    return print_and_free(ssc_run_json(text.c_str(), nullptr));
  }
  if (cmd == "eval" && argc >= 3) {
    std::vector<double> args;
    for (int i = 3; i < argc; ++i) {
      double v;
      if (!ssc::text::parse_double(argv[i], v)) {
        std::cerr << "bad number: " << argv[i] << "\n";
        return 1;
      }
      args.push_back(v);
    }
    double out = 0;
    const int rc = ssc_eval(argv[2], args.data(), static_cast<int>(args.size()), &out);
    if (rc != 0) {
      std::cerr << "unknown function or wrong number of arguments\n";
      return 1;
    }
    std::cout << ssc::text::format_g(out, 17) << "\n";
    return 0;
  }
  return usage();
}
