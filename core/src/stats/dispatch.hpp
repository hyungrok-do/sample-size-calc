#pragma once

#include <string>
#include <vector>

namespace ssc::stats {

/// Evaluate a distribution function by name (see dispatch.cpp for the list).
bool eval(const std::string& name, const double* args, int nargs, double& out);

/// Names of all functions available through eval().
std::vector<std::string> function_names();

}  // namespace ssc::stats
