#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "DynamicDensity.pb.h"
#include "cpp/DynamicKDE.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using dhist::DynamicKDE;
using dhist::Kernel;
using dynamic_density::DensityMap;
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

TEST(KernelTest, cdf) {
  Kernel k;
  k.addValue(1.0, 0.5);
  k.addValue(3.0, 0.5);
  k.addValue(5.0, 0.5);
  k.addValue(7.0, 0.5);

  EXPECT_EQ(k.cdf(4.0), 0.5);

  // In R:
  // > pnorm(5, 4, sqrt(5))
  // [1] 0.6726396
  EXPECT_NEAR(k.cdf(5.0), 0.6726396, 1e-6);

  EXPECT_NEAR(k.cdf(-100.0), 0.0, 1e-10);
  EXPECT_NEAR(k.cdf(100.0), 1.0, 1e-10);
}

TEST(KernelTest, populateProto) {
  Kernel k;
  k.addValue(1.0, 0.5);
  k.addValue(3.0, 0.5);
  k.addValue(5.0, 0.5);
  k.addValue(7.0, 0.5);

  ::dynamic_density::DynamicKDE::Kernel proto;
  k.populateProto(&proto);
  EXPECT_EQ(proto.coord()[0], 4.0);
  EXPECT_EQ(proto.variance()[0], 5.0);
  EXPECT_EQ(proto.count(), 2.0);
}

TEST(DynamicKDETest, addNoDecay) {
  static constexpr int kNumValues = 100000;
  static constexpr double kDecayRate = 0.0;
  DynamicKDE uut(/*num_kernels=*/31, /*decay_rate=*/kDecayRate);
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(norm(gen));
  }

  EXPECT_EQ(uut.computeTotalCount(), kNumValues);
  EXPECT_NEAR(uut.getQuantileEstimate(0.5), 0.0, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.05), -1.644854, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.95), 1.644854, 1e-1);

  EXPECT_NEAR(uut.getMean(), 0.0, 1e-1) << uut.debugString();
}

TEST(DynamicKDETest, addWithDecay) {
  static constexpr int kNumValues = 100000;
  static constexpr double kDecayRate = 0.00001;
  DynamicKDE uut(/*num_kernels=*/31, /*decay_rate=*/kDecayRate);
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(norm(gen));
  }

  EXPECT_NEAR(uut.computeTotalCount(),
              exponential_sum(1, 1 - kDecayRate, kNumValues), 1e-6);
  EXPECT_NEAR(uut.getQuantileEstimate(0.5), 0.0, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.05), -1.644854, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.95), 1.644854, 1e-1);

  EXPECT_NEAR(uut.getMean(), 0.0, 1e-1) << uut.debugString();
}

TEST(DynamicKDETest, toProto) {
  DynamicKDE uut(/*num_kernels=*/61);
  std::normal_distribution<double> norm(10000.0, 1.0);
  std::default_random_engine gen;

  static constexpr int kNumValues = 1000000;
  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(norm(gen));
  }

  DensityMap dm = uut.toProto("test", "x-value");

  EXPECT_EQ(dm.dynamic_kde().title(), "test");
  EXPECT_EQ(dm.dynamic_kde().label()[0], "x-value");

  size_t size = dm.dynamic_kde().kernels().size();
  EXPECT_EQ(size, 61);
  EXPECT_LT(dm.dynamic_kde().kernels()[0].coord()[0], 10000.0);
  EXPECT_GT(dm.dynamic_kde().kernels()[size - 1].coord()[0], 10000.0);

  // std::ofstream myfile("/tmp/DynamicKDE.pbuf");
  // ASSERT_TRUE(myfile.is_open());
  // ASSERT_TRUE(dm.SerializeToOstream(&myfile));
  // myfile.close();
}

TEST(DynamicKDETest, count) {
  DynamicKDE uut(/*num_kernels=*/61);
  std::normal_distribution<double> norm(10000.0, 1.0);
  std::default_random_engine gen;

  static constexpr int kNumValues = 1000000;
  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(norm(gen));
    ASSERT_EQ(uut.computeTotalCount(), i + 1) << uut.debugString();
  }
}
