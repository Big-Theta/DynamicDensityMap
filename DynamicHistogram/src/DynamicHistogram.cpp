#include "DynamicHistogram.h"

namespace dhist {

bool in_range(double val, double a, double b) {
  return ((val <= a) ^ (val <= b)) || val == a || val == b;
}

} // namespace dhist
