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

using dhist::Kernel;
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

TEST(KernelTest, addValue) {
  Kernel k;
  k.addValue(1.0, 1.0);
  k.addValue(3.0, 1.0);
  k.addValue(5.0, 1.0);
  k.addValue(7.0, 1.0);

  EXPECT_EQ(k.mean(), 4.0);
  EXPECT_EQ(k.variance(), 5.0);
  EXPECT_EQ(k.count(), 4.0);
}

TEST(KernelTest, decay) {
  Kernel k;
  k.addValue(1.0, 1.0);
  k.addValue(3.0, 1.0);
  k.addValue(5.0, 1.0);
  k.addValue(7.0, 1.0);
  k.decay(0.9, 5);

  EXPECT_EQ(k.mean(), 4.0);
  EXPECT_EQ(k.variance(), 5.0);
  EXPECT_EQ(k.count(), 3.6);
}

TEST(KernelTest, weightedAddValue) {
  Kernel k;
  k.addValue(1.0, 0.5);
  k.addValue(3.0, 0.5);
  k.addValue(5.0, 0.5);
  k.addValue(7.0, 0.5);

  EXPECT_EQ(k.mean(), 4.0);
  EXPECT_EQ(k.variance(), 5.0);
  EXPECT_EQ(k.count(), 2.0);
}

TEST(DynamicKDETest, addNoDecay) {
  static constexpr int kNumValues = 100000;
  static constexpr double kDecayRate = 0.0;
  DynamicKDE uut(/*max_num_buckets=*/31, /*decay_rate=*/kDecayRate);
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  setbuf(stdout, 0);
  printf("here\n");

  for (int i = 0; i < kNumValues; i++) {
    printf("here %d\n", i);
    uut.addValue(norm(gen));
  }
  printf("there\n%s\n");

  EXPECT_EQ(uut.computeTotalCount(), kNumValues);
  printf("there 2\n");
  EXPECT_NEAR(uut.getQuantileEstimate(0.5), 0.0, 1e-1);
  printf("there 3\n");
  EXPECT_NEAR(uut.getQuantileEstimate(0.05), -1.644854, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.95), 1.644854, 1e-1);

  EXPECT_NEAR(uut.getMean(), 0.0, 1e-1) << uut.debugString();
}
