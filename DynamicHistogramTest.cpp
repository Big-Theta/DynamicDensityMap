#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "DynamicDensity.pb.h"
#include "DynamicHistogram.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using dyden::DynamicHistogram;
using dyden::DynamicHistogramOpts;
using dyden::in_range;
using dynamic_density::DensityMap;
using testing::ContainerEq;
using testing::Eq;
using testing::Types;

TEST(InRange, in_range) {
  struct tcase {
    double val;
    double a;
    double b;
    bool in_range_expected;
  };

  std::vector<struct tcase> cases{
      {0.0, -1.0, 1.0, true},  {0.0, 1.0, -1.0, true},
      {1.0, 1.0, -1.0, true},  {1.0, -1.0, 1.0, true},
      {-1.0, 1.0, -1.0, true}, {-1.0, -1.0, 1.0, true},
      {2.0, 1.0, -1.0, false}, {-2.0, 1.0, -1.0, false},
      {-1.0, 0.0, 0.0, false}, {1.0, 0.0, 0.0, false},
      {0.0, 0.0, 0.0, true}};
  for (auto tcase : cases) {
    EXPECT_THAT(in_range(tcase.val, tcase.a, tcase.b),
                Eq(tcase.in_range_expected))
        << "in_range(" << tcase.val << ", " << tcase.a << ", " << tcase.b
        << ") == " << in_range(tcase.val, tcase.a, tcase.b)
        << "; expected = " << tcase.in_range_expected;
  }
}

// From https://en.wikipedia.org/wiki/Geometric_series#Formula
// Args:
//   a: First term in series.
//   r: Rate. Decay rate is 1 - rate in this formula.
//   n: Number of terms in series.
double exponential_sum(double a, double r, int n) {
  return a * (1.0 - pow(r, n)) / (1 - r);
}

TEST(DynamicHistogramTest, addNoDecay) {
  static constexpr int kNumValues = 1000;
  static constexpr double kDecayRate = 0.0;
  DynamicHistogram uut(
      DynamicHistogramOpts().set_num_buckets(31).set_decay_rate(kDecayRate));
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(norm(gen));
  }

  EXPECT_EQ(uut.computeTotalCount(), kNumValues);
  EXPECT_NEAR(uut.getQuantileEstimate(0.5), 0.0, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.05), -1.644854, 1e-1)
      << uut.debugString();
  EXPECT_NEAR(uut.getQuantileEstimate(0.95), 1.644854, 1e-1)
      << uut.debugString();
}

TEST(DynamicHistogramTest, addWithDecay) {
  static constexpr int kNumValues = 100000;
  static constexpr double kDecayRate = 0.00001;
  DynamicHistogram uut(
      DynamicHistogramOpts().set_num_buckets(31).set_decay_rate(kDecayRate));
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(norm(gen));
  }

  EXPECT_NEAR(uut.computeTotalCount(),
              exponential_sum(1, 1 - kDecayRate, kNumValues), 1e-6);

  // Using R:
  // $ R -q
  // > qnorm(0.5, 0, 1)
  // [1] 0
  // > qnorm(0.05, 0, 1)
  // [1] -1.644854
  // > qnorm(0.95, 0, 1)
  // [1] 1.644854
  EXPECT_NEAR(uut.getQuantileEstimate(0.5), 0.0, 1e-1) << uut.debugString();
  EXPECT_NEAR(uut.getQuantileEstimate(0.05), -1.644854, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.95), 1.644854, 1e-1);
}

TEST(DynamicHistogramTest, multithreadedStress) {
  static const int kInsertsPerThread = 100000;
  static const int kNumThreads = 8;
  std::thread threads[kNumThreads];

  DynamicHistogram uut(
      DynamicHistogramOpts().set_num_buckets(31).set_decay_rate(0.00001));

  for (int tx = 0; tx < kNumThreads; tx++) {
    threads[tx] = std::thread(
        [&](int tx) {
          std::default_random_engine gen;
          gen.seed(tx);
          std::normal_distribution<double> norm(0.0, 1.0);
          std::bernoulli_distribution bern(0.01);

          for (int i = 0; i < kInsertsPerThread; i++) {
            uut.addValue(norm(gen));

            if (bern(gen)) {
              EXPECT_GT(uut.json().size(), 0);
            }

            if (i > 100000 && bern(gen)) {
              EXPECT_NEAR(uut.getQuantileEstimate(0.5), 0, 1.0);
            }
          }
        },
        tx);
  }

  for (int tx = 0; tx < kNumThreads; tx++) {
    threads[tx].join();
  }
}

TEST(DynamicHistogramTest, getMean) {
  DynamicHistogram uut(DynamicHistogramOpts().set_num_buckets(61));
  std::normal_distribution<double> norm(10000.0, 1.0);
  std::default_random_engine gen;

  static constexpr int kNumValues = 100000;
  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(norm(gen));
  }

  EXPECT_NEAR(uut.getMean(), 10000.0, 1e-1);
}

TEST(DynamicHistogramTest, asProto) {
  DynamicHistogram uut(DynamicHistogramOpts().set_num_buckets(61));
  uut.mutable_description()->set_title("test");
  uut.mutable_description()->set_labels({"x-value"});
  std::normal_distribution<double> norm(10000.0, 1.0);
  std::default_random_engine gen;

  static constexpr int kNumValues = 100000;
  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(norm(gen));
  }

  DensityMap dm = uut.asProto();

  EXPECT_EQ(dm.dynamic_histogram().description().title(), "test");
  EXPECT_EQ(dm.dynamic_histogram().description().labels()[0], "x-value");
  EXPECT_EQ(dm.dynamic_histogram().bounds()[0], uut.getMin());
  EXPECT_EQ(
      dm.dynamic_histogram().bounds()[dm.dynamic_histogram().bounds_size() - 1],
      uut.getMax());

  std::ofstream myfile("/tmp/DynamicHistogram.pbuf");
  ASSERT_TRUE(myfile.is_open());
  ASSERT_TRUE(dm.SerializeToOstream(&myfile));
  myfile.close();
}
