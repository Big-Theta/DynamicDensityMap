#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "DynamicDensity.pb.h"
#include "DynamicKDE2D.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using dyden::DynamicKDE2D;
using dyden::DynamicKDE2DOpts;
using dyden::Kernel2D;
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
  Kernel2D k;
  k.addValue(1.0, 10.0, 1.0);
  k.addValue(3.0, 30.0, 1.0);
  k.addValue(5.0, 50.0, 1.0);
  k.addValue(7.0, 70.0, 1.0);

  EXPECT_EQ(k.mean_x(), 4.0);
  EXPECT_EQ(k.mean_y(), 40.0);
  EXPECT_EQ(k.variance_xx(), 5.0);
  EXPECT_EQ(k.variance_yy(), 500.0);
  EXPECT_EQ(k.covariance(), 50.0);
  EXPECT_EQ(k.count(), 4.0);
}

TEST(KernelTest, decay) {
  Kernel2D k;
  k.addValue(1.0, 10.0, 1.0);
  k.addValue(3.0, 30.0, 1.0);
  k.addValue(5.0, 50.0, 1.0);
  k.addValue(7.0, 70.0, 1.0);
  k.decay(0.9, 5);

  EXPECT_EQ(k.mean_x(), 4.0);
  EXPECT_EQ(k.mean_y(), 40.0);
  EXPECT_EQ(k.variance_xx(), 5.0);
  EXPECT_EQ(k.variance_yy(), 500.0);
  EXPECT_EQ(k.covariance(), 50.0);
  EXPECT_EQ(k.count(), 3.6);
}

TEST(KernelTest, weightedAddValue) {
  Kernel2D k;
  k.addValue(1.0, 10.0, 0.5);
  k.addValue(3.0, 30.0, 0.5);
  k.addValue(5.0, 50.0, 0.5);
  k.addValue(7.0, 70.0, 0.5);

  EXPECT_EQ(k.mean_x(), 4.0);
  EXPECT_EQ(k.mean_y(), 40.0);
  EXPECT_EQ(k.variance_xx(), 5.0);
  EXPECT_EQ(k.variance_yy(), 500.0);
  EXPECT_EQ(k.covariance(), 50.0);
  EXPECT_EQ(k.count(), 2.0);
}

TEST(KernelTest, populateProto) {
  Kernel2D k;
  k.addValue(1.0, 10.0, 0.5);
  k.addValue(3.0, 30.0, 0.5);
  k.addValue(5.0, 50.0, 0.5);
  k.addValue(7.0, 70.0, 0.5);

  ::dynamic_density::DynamicKDE::Kernel proto;
  k.populateProto(&proto);
  EXPECT_EQ(proto.coord()[0], 4.0);
  EXPECT_EQ(proto.coord()[1], 40.0);
  EXPECT_EQ(proto.variance()[0], 5.0);
  EXPECT_EQ(proto.variance()[1], 500.0);
  EXPECT_EQ(proto.covariance(), 50.0);
  EXPECT_EQ(proto.count(), 2.0);
}

TEST(DynamicKDE2DTest, addNoDecay) {
  static constexpr int kNumValues = 1000000;
  static constexpr double kDecayRate = 0.0;
  DynamicKDE2D uut(
      DynamicKDE2DOpts().set_num_kernels(200).set_decay_rate(kDecayRate));
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);
  std::uniform_real_distribution unif(0.0, 1.0);
  std::discrete_distribution<int> discrete{1, 1, 1, 1};

  setbuf(stdout, 0);

  for (int i = 0; i < kNumValues; i++) {
    switch (discrete(gen)) {
    case 0:
      uut.addValue(norm(gen), norm(gen));
      break;
    case 1:
      uut.addValue(norm(gen) + 1, norm(gen) + 1);
      break;
    case 2:
      uut.addValue(norm(gen) - 10, norm(gen) - 7);
      break;
    default:
      uut.addValue(2 * unif(gen) + 10, 3 * unif(gen) + 10);
      break;
    }
    ASSERT_EQ(uut.computeTotalCount(), i + 1);
  }

  printf("??? %s\n", uut.asProto().DebugString().c_str());

  EXPECT_EQ(uut.computeTotalCount(), kNumValues);

  std::ofstream myfile("/tmp/DynamicKDE2D.pbuf");
  ASSERT_TRUE(myfile.is_open());
  ASSERT_TRUE(uut.asProto().SerializeToOstream(&myfile));
  myfile.close();
}
