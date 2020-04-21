#include <cstddef>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "DynamicHistogram.h"


using testing::Eq;

class DynamicHistogramReference {
public:
  DynamicHistogramReference(double decay_rate, size_t max_num_buckets)
      : decay_rate_(decay_rate), max_num_buckets_(max_num_buckets),
        generation_(0) {
    ubounds_.resize(max_num_buckets_ - 1);
    counts_.resize(max_num_buckets_);
  }

  virtual ~DynamicHistogramReference() {}

  void addValue(double val) {
    // Decay
    for (int i = 0; i < counts_.size(); i++) {
      counts_[i] = counts_[i] * decay_rate_;
    }

    // Adjust bounds
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

    // Decide whether to split
    if (counts_[bucket_idx] < splitThreshold()) {
      return;
    }

    // Split
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
      for (int i = ubounds_.size() - 1; i >= bucket_idx; i--) {
        ubounds_[i] = ubounds_[i - 1];
      }

      ubounds_[bucket_idx - 1] = (lower_bound + upper_bound) / 2.0;

      counts_.push_back(0.0);
      for (int i = counts_.size() - 1; i > bucket_idx; i--) {
        counts_[i] = counts_[i - 1];
      }
      counts_[bucket_idx] /= 2.0;
      counts_[bucket_idx + 1] = counts_[bucket_idx];
    }

    // Merge
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

  void addRepeatedValue(double val) {}
  void clear() {}
  void merge(const DynamicHistogramReference &other) {}
  void copy(const DynamicHistogramReference &other) {}

  size_t getNumBuckets() { return counts_.size(); }
  Bucket getBucketByIndex(size_t idx) { return Bucket(0, 0, 0, 0); }

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

  double computeTotalCount() {
    double total_count = 0.0;
    for (double count : counts_) {
      total_count += count;
    }
    return total_count;
  }

  double getPercentileEstimate(double pct) {
    double total_count = computeTotalCount();
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
    s.resize(50 * counts_.size());
    int cursor = 0;

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
    cursor += snprintf(&s[cursor], s.size() - cursor, "{\n  \"bounds\": [%lf, ", getMin());
    for (auto bound : ubounds_) {
      cursor += snprintf(&s[cursor], s.size() - cursor, "%lf, ", bound);
    }
    cursor += snprintf(&s[cursor], s.size() - cursor, "%lf],\n  \"counts\": [", getMax());

    for (int i = 0; i < counts_.size() - 1; i++) {
      cursor += snprintf(&s[cursor], s.size() - cursor, "%lf, ", counts_[i]);
    }
    cursor += snprintf(&s[cursor], s.size() - cursor, "%lf]\n}\n", counts_.back());
    s.resize(cursor);
    return s;
  }

protected:
  const double decay_rate_;
  const size_t max_num_buckets_;

  uint64_t generation_;

  // ubounds_ records the upper bound of a bucket. There isn't an upper bound
  // for the last bucket, so the length of counts_ will generally be one greater than the
  // length of counts_.
  std::vector<double> ubounds_;
  std::vector<double> counts_;

  double splitThreshold() {
    double total_count = computeTotalCount();
    return 2 * total_count / getNumBuckets();
  }
};

TEST(DynamicHistogramTest, addNoDecay) {
  static constexpr int kNumValues = 100;
  static constexpr double kVal = 100.0;
  DynamicHistogramReference uut(/*decay_rate=*/1.0, 10);

  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(kVal);
    ASSERT_TRUE(uut.isValid()) << uut.debugString();
    ASSERT_EQ(uut.computeTotalCount(), i + 1);
  }

  EXPECT_EQ(uut.computeTotalCount(), kNumValues) << uut.debugString();
  EXPECT_EQ(uut.getPercentileEstimate(0.5), kVal) << uut.json();
}
