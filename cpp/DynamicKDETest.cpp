#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "DensityMap.pb.h"
#include "cpp/DynamicKDE.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using dhist::DynamicKDE;
using dynamic_histogram::DensityMap;
using testing::ContainerEq;
using testing::Eq;
using testing::Types;

// From https://en.wikipedia.org/wiki/Geometric_series#Formula
// Args:
//   a: First term in series.
//   r: Rate. Decay rate is 1 - rate in this formula.
//   n: Number of terms in series.
double exponential_sum(double a, double r, int n) {
  return a * (1.0 - pow(r, n)) / (1 - r);
}

TEST(ReferenceTest, addNoDecay) {
  static constexpr int kNumValues = 100000;
  static constexpr double kDecayRate = 0.0;
  DynamicKDE uut(/*max_num_buckets=*/31, /*decay_rate=*/kDecayRate);
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(norm(gen));
  }

  EXPECT_EQ(uut.computeTotalCount(), kNumValues);
  EXPECT_NEAR(uut.getQuantileEstimate(0.5), 0.0, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.05), -1.644854, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.95), 1.644854, 1e-1);

  EXPECT_NEAR(uut.getMean(), 0.0, 1e-1);
}