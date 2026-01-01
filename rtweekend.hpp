#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

// C++ Std Usings

using std::make_shared;
using std::shared_ptr;

// Constants

constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr double pi = 3.1415926535897932385;

// Utility Functions

inline double degrees_to_radians(const double degrees)
{
  return degrees * pi / 180.0;
}

// Common Headers

#include "color.hpp"
#include "ray.hpp"
#include "vec3.hpp"