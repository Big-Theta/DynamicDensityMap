#include <algorithm>
#include <cstddef>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "cpp/DynamicHistogram.h"
#include "cpp/DynamicHistogramReference.h"

using dhist::DynamicHistogram;
using dhist::DynamicHistogramReference;
using dhist::in_range;
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

TEST(ReferenceTest, addNoDecay) {
  static constexpr int kNumValues = 100000;
  static constexpr double kDecayRate = 0.0;
  DynamicHistogramReference uut(
      /*max_num_buckets=*/31,
      /*decay_rate=*/kDecayRate);
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

TEST(ReferenceTest, addWithDecay) {
  static constexpr int kNumValues = 100000;
  static constexpr double kDecayRate = 0.00001;
  DynamicHistogramReference uut(
      /*max_num_buckets=*/31,
      /*decay_rate=*/kDecayRate);
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
  EXPECT_NEAR(uut.getQuantileEstimate(0.5), 0.0, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.05), -1.644854, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.95), 1.644854, 1e-1);
}

TEST(ReferenceTest, addRandomNoDecay) {
  static constexpr int kNumValues = 100000;
  DynamicHistogramReference uut(
      /*max_num_buckets=*/31,
      /*decay_rate=*/0.0);
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(norm(gen));
  }

  EXPECT_NEAR(uut.getQuantileEstimate(0.5), 0.0, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.05), -1.644854, 1e-1);
  EXPECT_NEAR(uut.getQuantileEstimate(0.95), 1.644854, 1e-1);
}

TEST(ReferenceTest, trackQuantiles) {
  DynamicHistogramReference uut(
      /*max_num_buckets=*/4,
      /*decay_rate=*/0.0,
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

TEST(ReferenceTest, quantilesTracking) {
  static constexpr int kNumValues = 100000;
  DynamicHistogramReference uut(
      /*max_num_buckets=*/31,
      /*decay_rate=*/0.0);
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

TEST(DynamicHistogramTest, addNoDecay) {
  static constexpr int kNumValues = 100000;
  static constexpr double kDecayRate = 0.0;
  DynamicHistogram</*kUseDecay=*/false, /*kThreadsafe=*/false> uut(
      /*max_num_buckets=*/31,
      /*decay_rate=*/kDecayRate);
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
  DynamicHistogram</*kUseDecay=*/true, /*kThreadsafe=*/false> uut(
      /*max_num_buckets=*/31,
      /*decay_rate=*/kDecayRate);
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

  DynamicHistogram</*kUseDecay=*/true, /*kThreadsafe=*/true> uut(
      /*max_num_buckets=*/31, /*decay_rate=*/0.00001);

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

struct TypedTestTrue {
  static constexpr bool val = true;
};

struct TypedTestFalse {
  static constexpr bool val = false;
};

template <typename T> class DynamicHistogramTypedTest : public testing::Test {};

using test_types = testing::Types<TypedTestTrue, TypedTestFalse>;
TYPED_TEST_SUITE(DynamicHistogramTypedTest, test_types);

TYPED_TEST(DynamicHistogramTypedTest, referenceEquivalence) {
  static constexpr double kDecayRate = 0.0001;
  static constexpr int kMaxNumBuckets = 31;
  static constexpr int kNumValues = 1;

  TypeParam kThreadsafe;

  DynamicHistogramReference ref(
      /*max_num_buckets=*/kMaxNumBuckets,
      /*decay_rate=*/kDecayRate);
  DynamicHistogram</*kUseDecay=*/true, /*kThreadsafe=*/kThreadsafe.val> dyn(
      /*max_num_buckets=*/kMaxNumBuckets,
      /*decay_rate=*/kDecayRate);

  std::default_random_engine gen;
  gen.seed(1);
  std::normal_distribution<double> norm(0.0, 1.0);
  for (int i = 0; i < kNumValues; i++) {
    double val = norm(gen);
    ref.addValue(val);
    dyn.addValue(val);
  }

  EXPECT_NEAR(ref.computeTotalCount(), dyn.computeTotalCount(), 1e-6);
  for (int i = 1; i < 100; i++) {
    ASSERT_NEAR(ref.getQuantileEstimate(i / 100.0),
                dyn.getQuantileEstimate(i / 100.0), 1e-6)
        << "(tested percentile: " << i / 100.0 << ")";
  }

  EXPECT_EQ(ref.getNumBuckets(), dyn.getNumBuckets());
  for (int bx = 0; bx < ref.getNumBuckets(); bx++) {
    auto ref_bucket = ref.getBucketByIndex(bx);
    auto dyn_bucket = dyn.getBucketByIndex(bx);
    EXPECT_NEAR(ref_bucket.min(), dyn_bucket.min(), 1e-6);
    EXPECT_NEAR(ref_bucket.max(), dyn_bucket.max(), 1e-6);
    EXPECT_NEAR(ref_bucket.count(), dyn_bucket.count(), 1e-6);
  }
}
