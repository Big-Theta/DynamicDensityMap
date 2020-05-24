#include <algorithm>
#include <cstddef>
#include <random>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "DynamicHistogram.h"
#include "DynamicHistogramReference.h"

using dhist::DynamicHistogramReference;
using dhist::in_range;
using testing::ContainerEq;
using testing::Eq;

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
  static constexpr int kNumValues = 100000;
  static constexpr double kDecayRate = 0.0;
  static constexpr double kMean = 0.0;
  static constexpr double kStdDev = 1.0;
  DynamicHistogramReference uut(/*decay_rate=*/kDecayRate,
                                /*max_num_buckets=*/31);
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(norm(gen));
  }

  EXPECT_EQ(uut.computeTotalCount(), kNumValues);
  EXPECT_NEAR(uut.getQuantileEstimate(0.5), 0.0, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.05), -1.644854, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.95), 1.644854, 1e-1);
}

TEST(DynamicHistogramTest, addWithDecay) {
  static constexpr int kNumValues = 100000;
  static constexpr double kDecayRate = 0.00001;
  static constexpr double kMean = 0.0;
  static constexpr double kStdDev = 1.0;
  DynamicHistogramReference uut(/*decay_rate=*/kDecayRate,
                                /*max_num_buckets=*/31);
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(norm(gen));
  }

  EXPECT_NEAR(uut.computeTotalCount(),
              exponential_sum(1, 1 - kDecayRate, kNumValues), 1e-10);

  // Using R:
  // $ R -q
  // > qnorm(0.5, 0, 1)
  // [1] 0
  // > qnorm(0.05, 0, 1)
  // [1] -1.644854
  // > qnorm(0.95, 0, 1)
  // [1] 1.644854
  EXPECT_NEAR(uut.getQuantileEstimate(0.5), 0.0, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.05), -1.644854, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.95), 1.644854, 1e-1);
}

TEST(DynamicHistogramTest, addRandomNoDecay) {
  static constexpr int kNumValues = 100000;
  DynamicHistogramReference uut(/*decay_rate=*/0.0,
                                /*max_num_buckets=*/31);
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(norm(gen));
  }

  EXPECT_NEAR(uut.getQuantileEstimate(0.5), 0.0, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.05), -1.644854, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.95), 1.644854, 1e-1);
}

TEST(DynamicHistogramTest, trackQuantiles) {
  static constexpr int kNumValues = 1000000;
  DynamicHistogramReference uut(
      /*decay_rate=*/0.0, /*max_num_buckets=*/4,
      /*ubounds=*/std::vector<double>({0.25, 0.5, 0.75}),
      /*counts=*/std::vector<double>({100.0, 100.0, 100.0, 100.0}));
  uut.trackQuantiles(
      /*quantiles=*/std::vector<double>({0.01, 0.1, 0.25, 0.5, 0.9, 0.99}));

  std::map<double, double> quantiles = uut.getTrackedQuantiles();
  EXPECT_THAT(uut.getTrackedQuantiles(),
              ContainerEq(std::map<double, double>{{0.01, 0.01},
                                                   {0.1, 0.1},
                                                   {0.25, 0.25},
                                                   {0.5, 0.5},
                                                   {0.9, 0.9},
                                                   {0.99, 0.99}}));
}

TEST(DynamicHistogramTest, quantilesTracking) {
  static constexpr int kNumValues = 100000;
  DynamicHistogramReference uut(
      /*decay_rate=*/0.0, /*max_num_buckets=*/31);
  uut.trackQuantiles(
      /*quantiles=*/std::vector<double>({0.01, 0.1, 0.25, 0.5, 0.9, 0.99}));

  std::map<double, double> quantiles = uut.getTrackedQuantiles();
  EXPECT_THAT(uut.getTrackedQuantiles(),
              ContainerEq(std::map<double, double>{{0.01, 0.0},
                                                   {0.1, 0.0},
                                                   {0.25, 0.0},
                                                   {0.5, 0.0},
                                                   {0.9, 0.0},
                                                   {0.99, 0.0}}));

  std::default_random_engine gen;
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(unif(gen));
  }

  for (const auto &kv : uut.getTrackedQuantiles()) {
    EXPECT_NEAR(kv.first, kv.second, 0.01);
  }
}
