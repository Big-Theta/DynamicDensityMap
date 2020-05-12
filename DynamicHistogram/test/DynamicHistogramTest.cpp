#include <algorithm>
#include <cstddef>
#include <random>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "DynamicHistogram.h"

using testing::ContainerEq;
using testing::Eq;

bool in_range(double val, double a, double b) {
  return (val <= a) ^ (val <= b) || val == a || val == b;
}

TEST(InRange, inRange) {
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
  };
  for (auto tcase : cases) {
    EXPECT_THAT(in_range(tcase.val, tcase.a, tcase.b),
                Eq(tcase.in_range_expected))
        << "in_range(" << tcase.val << ", " << tcase.a << ", " << tcase.b
        << ") == " << in_range(tcase.val, tcase.a, tcase.b)
        << "; expected = " << tcase.in_range_expected;
    ;
  }
}

class DynamicHistogramReference {
public:
  DynamicHistogramReference(double decay_rate, size_t max_num_buckets)
      : decay_rate_(decay_rate), max_num_buckets_(max_num_buckets),
        generation_(0) {
    ubounds_.resize(max_num_buckets_ - 1);
    counts_.resize(max_num_buckets_);
  }

  DynamicHistogramReference(double decay_rate, size_t max_num_buckets,
                            const std::vector<double> &quantiles)
      : DynamicHistogramReference(decay_rate, max_num_buckets) {
    std::copy(quantiles.begin(), quantiles.end(), back_inserter(quantiles_));
  }

  virtual ~DynamicHistogramReference() {}

  void addValue(double val) {
    generation_++;
    decay();
    // shift_quantiles will treat the new value as a point with mass 1 which
    // only makes sense after the decay and before the point has been
    // integrated into the histogram.
    shift_quantiles(val);
    int bucket_idx = insert_value(val);

    if (counts_[bucket_idx] < splitThreshold()) {
      return;
    }

    split(bucket_idx);
    merge();
  }

  void addRepeatedValue(double val) {}
  void clear() {}
  void merge(const DynamicHistogramReference &other) {}
  void copy(const DynamicHistogramReference &other) {}

  size_t getNumBuckets() { return counts_.size(); }

  Bucket getBucketByIndex(size_t idx) {
    double min;
    if (idx > 0) {
      min = ubounds_[idx - 1];
    } else {
      min = getMin();
    }

    double max;
    if (idx < ubounds_.size()) {
      max = ubounds_[idx];
    } else {
      max = getMax();
    }

    return Bucket(/*min=*/min, /*max=*/max, /*count=*/counts_[idx]);
  }

  double getMin() const {
    if (counts_[1] == 0) {
      return ubounds_[0];
    }
    return ubounds_[0] -
           (counts_[0] / counts_[1]) * (ubounds_[1] - ubounds_[0]);
  }

  double getMax() const {
    const int idx = ubounds_.size() - 1;
    if (counts_[idx] == 0.0) {
      return ubounds_.back();
    }
    return ubounds_[idx] + (counts_[idx + 1] / counts_[idx]) *
                               (ubounds_[idx] - ubounds_[idx - 1]);
  }

  double computeTotalCount() const {
    double total_count = 0.0;
    for (double count : counts_) {
      total_count += count;
    }
    return total_count;
  }

  // pct is in [0, 1]
  double getPercentileEstimate(double pct) {
    double total_count = computeTotalCount();
    if (total_count == 0.0) {
      return 0.0;
    }

    double cdf = 0.0;
    double next_cdf = 0.0;
    int bucket_idx = 0;
    for (int i = 0; i < counts_.size(); i++) {
      next_cdf = cdf + counts_[i] / total_count;
      if (next_cdf > pct) {
        bucket_idx = i;
        break;
      }
      cdf = next_cdf;
    }

    double frac = (pct - cdf) / (next_cdf - cdf);

    double lower_bound;
    if (bucket_idx == 0) {
      lower_bound = getMin();
    } else {
      lower_bound = ubounds_[bucket_idx - 1];
    }

    double upper_bound;
    if (bucket_idx + 1 == ubounds_.size()) {
      upper_bound = getMax();
    } else {
      upper_bound = ubounds_[bucket_idx];
    }

    return frac * (upper_bound - lower_bound) + lower_bound;
  }

  std::map<double, double> getTrackedQuantiles() {
    std::map<double, double> qs;
    for (auto quantile : quantiles_) {
      qs[quantile] = getPercentileEstimate(quantile);
    }
    return qs;
  }

  bool isValid() const {
    if (ubounds_.size() + 1 != counts_.size()) {
      return false;
    }

    if (getMin() > ubounds_[0]) {
      return false;
    }

    for (int i = 0; i < ubounds_.size() - 1; i++) {
      if (ubounds_[i] > ubounds_[i + 1]) {
        return false;
      }
    }

    if (ubounds_.back() > getMax()) {
      return false;
    }

    for (auto count : counts_) {
      if (count < 0.0) {
        return false;
      }
    }

    return true;
  }

  std::string debugString() const {
    std::string s;
    s.resize(50 * (counts_.size() + 1));
    int cursor = 0;

    cursor += snprintf(&s[cursor], s.size() - cursor,
                       "generation: %lu\ntotal_count: %lf\n", generation_,
                       computeTotalCount());

    cursor += snprintf(&s[cursor], s.size() - cursor, "  %d [%lf, %lf): %lf\n",
                       0, getMin(), ubounds_[1], counts_[0]);
    int i;
    for (i = 0; i < ubounds_.size() - 1; i++) {
      cursor +=
          snprintf(&s[cursor], s.size() - cursor, "  %d [%lf, %lf): %lf\n",
                   i + 1, ubounds_[i], ubounds_[i + 1], counts_[i + 1]);
    }
    cursor += snprintf(&s[cursor], s.size() - cursor, "  %d [%lf, %lf): %lf\n",
                       i + 1, ubounds_.back(), getMax(), counts_[i + 1]);
    s.resize(cursor);

    return s;
  }

  std::string json() const {
    std::string s;
    s.resize((counts_.size() + ubounds_.size() + 2) * 16);
    int cursor = 0;
    cursor += snprintf(&s[cursor], s.size() - cursor, "{\n  \"bounds\": [%lf, ",
                       getMin());
    for (auto bound : ubounds_) {
      cursor += snprintf(&s[cursor], s.size() - cursor, "%lf, ", bound);
    }
    cursor += snprintf(&s[cursor], s.size() - cursor, "%lf],\n  \"counts\": [",
                       getMax());

    for (int i = 0; i < counts_.size() - 1; i++) {
      cursor += snprintf(&s[cursor], s.size() - cursor, "%lf, ", counts_[i]);
    }
    cursor +=
        snprintf(&s[cursor], s.size() - cursor, "%lf]\n}\n", counts_.back());
    s.resize(cursor);
    return s;
  }

protected:
  const double decay_rate_;
  const size_t max_num_buckets_;

  uint64_t generation_;

  // ubounds_ records the upper bound of a bucket. There isn't an upper bound
  // for the last bucket, so the length of counts_ will generally be one greater
  // than the length of counts_.
  std::vector<double> ubounds_;
  std::vector<double> counts_;
  std::vector<double> quantiles_;

  double splitThreshold() {
    double total_count = computeTotalCount();
    return 2 * total_count / getNumBuckets();
  }

  void decay() {
    for (int i = 0; i < counts_.size(); i++) {
      counts_[i] = counts_[i] * decay_rate_;
    }
  }

  // Adjust bounds and add.
  // Returns:
  //   The index of the bucket that the new value landed in.
  int insert_value(double val) {
    int bucket_idx;
    for (bucket_idx = 0; bucket_idx < ubounds_.size(); bucket_idx++) {
      if (val < ubounds_[bucket_idx]) {
        break;
      }
    }

    if (bucket_idx > 0) {
      double count_with_below = counts_[bucket_idx - 1] + counts_[bucket_idx];
      ubounds_[bucket_idx - 1] =
          (count_with_below * ubounds_[bucket_idx - 1] + val) /
          (count_with_below + 1.0);
    }

    if (bucket_idx < ubounds_.size()) {
      double count_with_above = counts_[bucket_idx] + counts_[bucket_idx + 1];
      ubounds_[bucket_idx] = (count_with_above * ubounds_[bucket_idx] + val) /
                             (count_with_above + 1.0);
    }

    // Add value
    counts_[bucket_idx] += 1.0;

    return bucket_idx;
  }

  void split(int bucket_idx) {
    double lower_bound;
    if (bucket_idx == 0) {
      lower_bound = getMin();
    } else {
      lower_bound = ubounds_[bucket_idx - 1];
    }

    if (bucket_idx == ubounds_.size()) {
      double upper_bound = getMax();
      ubounds_.push_back((lower_bound + upper_bound) / 2.0);
      counts_[bucket_idx] /= 2.0;
      counts_.push_back(counts_.back());
    } else {
      double upper_bound = ubounds_[bucket_idx];

      ubounds_.push_back(0.0);
      for (int i = ubounds_.size() - 1; i > bucket_idx; i--) {
        ubounds_[i] = ubounds_[i - 1];
      }

      ubounds_[bucket_idx] = (lower_bound + upper_bound) / 2.0;

      counts_.push_back(0.0);
      for (int i = counts_.size() - 1; i > bucket_idx; i--) {
        counts_[i] = counts_[i - 1];
      }
      counts_[bucket_idx] /= 2.0;
      counts_[bucket_idx + 1] = counts_[bucket_idx];
    }
  }

  void merge() {
    int merge_idx = 0;
    double merged_count = counts_[0] + counts_[1];
    for (int i = 1; i < counts_.size() - 1; i++) {
      double pos_count = counts_[i] + counts_[i + 1];
      if (pos_count < merged_count) {
        merge_idx = i;
        merged_count = pos_count;
      }
    }

    for (int i = merge_idx; i < ubounds_.size() - 1; i++) {
      ubounds_[i] = ubounds_[i + 1];
    }
    ubounds_.pop_back();

    counts_[merge_idx] = merged_count;
    for (int i = merge_idx + 1; i < counts_.size() - 1; i++) {
      counts_[i] = counts_[i + 1];
    }
    counts_.pop_back();
  }

  void shift_quantiles(double val) {
    for (auto quantile : quantiles_) {

    }
  }

  // Shifts the quantile inside the bucket. Returns the amount that remains to
  // be shifted. If 0, the shift is done.
  double shift_quantile_in_bucket(int quantile_idx, int bucket_idx, double val,
                                  double shift_remaining) {
    if (shift_remaining == 0.0) {
      return 0.0;
    }

    Bucket bucket = getBucketByIndex(bucket_idx);

    double target_location = quantiles_[quantile_idx] + shift_remaining *
                                                            bucket.diameter() /
                                                            bucket.count();

    assert(false);
    return 0.0;
  }
};

// From https://en.wikipedia.org/wiki/Geometric_series#Formula
// Args:
//   a: First term in series.
//   r: Decay rate.
//   n: Number of terms in series.
double exponential_sum(double a, double r, int n) {
  return a * (1.0 - pow(r, n)) / (1 - r);
}

TEST(DynamicHistogramTest, addNoDecay) {
  static constexpr int kNumValues = 100;
  static constexpr double kVal = 100.0;
  DynamicHistogramReference uut(/*decay_rate=*/1.0, /*max_num_buckets=*/10);

  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(kVal);
    ASSERT_TRUE(uut.isValid()) << uut.debugString();
    ASSERT_EQ(uut.computeTotalCount(), i + 1);
  }

  EXPECT_EQ(uut.computeTotalCount(), kNumValues) << uut.debugString();
  EXPECT_EQ(uut.getPercentileEstimate(0.5), kVal) << uut.debugString();
}

TEST(DynamicHistogramTest, addWithDecay) {
  static constexpr int kNumValues = 100000;
  static constexpr double kDecayRate = 0.9999;
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
              exponential_sum(1, kDecayRate, kNumValues), 1e-10);

  // Using R:
  // $ R -q
  // > qnorm(0.5, 0, 1)
  // [1] 0
  // > qnorm(0.05, 0, 1)
  // [1] -1.644854
  // > qnorm(0.95, 0, 1)
  // [1] 1.644854
  EXPECT_NEAR(uut.getPercentileEstimate(0.5), 0.0, 1e-1);
  EXPECT_NEAR(uut.getPercentileEstimate(0.05), -1.644854, 1e-1);
  EXPECT_NEAR(uut.getPercentileEstimate(0.95), 1.644854, 1e-1);
}

TEST(DynamicHistogramTest, addRandomNoDecay) {
  static constexpr int kNumValues = 100000;
  DynamicHistogramReference uut(/*decay_rate=*/1.0,
                                /*max_num_buckets=*/31);
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(norm(gen));
  }

  EXPECT_NEAR(uut.getPercentileEstimate(0.5), 0.0, 1e-1);
  EXPECT_NEAR(uut.getPercentileEstimate(0.05), -1.644854, 1e-1);
  EXPECT_NEAR(uut.getPercentileEstimate(0.95), 1.644854, 1e-1);
}

TEST(DynamicHistogramTest, quantilesConstructor) {
  DynamicHistogramReference uut(
      /*decay_rate=*/0.9999, /*max_num_buckets=*/31,
      /*quantiles=*/std::vector<double>({0.01, 0.1, 0.25, 0.5, 0.9, 0.99}));

  std::map<double, double> quantiles = uut.getTrackedQuantiles();
  EXPECT_THAT(uut.getTrackedQuantiles(),
              ContainerEq(std::map<double, double>{{0.01, 0.0},
                                                   {0.1, 0.0},
                                                   {0.25, 0.0},
                                                   {0.5, 0.0},
                                                   {0.9, 0.0},
                                                   {0.99, 0.0}}));
}
