#include "DynamicHistogram.h"

namespace dhist {

class DynamicHistogramReference {
public:
  DynamicHistogramReference(size_t max_num_buckets, double decay_rate = 0.0)
      : max_num_buckets_(max_num_buckets), decay_rate_(decay_rate),
        generation_(0) {
    ubounds_.resize(max_num_buckets_ - 1);
    counts_.resize(max_num_buckets_);
  }

  DynamicHistogramReference(size_t max_num_buckets, double decay_rate,
                            const std::vector<double> &ubounds,
                            const std::vector<double> &counts)
      : DynamicHistogramReference(max_num_buckets, decay_rate) {
    assert(ubounds.size() + 1 == counts.size());

    ubounds_.clear();
    std::copy(ubounds.begin(), ubounds.end(), back_inserter(ubounds_));

    counts_.clear();
    std::copy(counts.begin(), counts.end(), back_inserter(counts_));
  }

  virtual ~DynamicHistogramReference() {}

  void addValue(double val) {
    generation_++;
    decay();
    // shift_quantiles will treat the new value as a point with mass 1 which
    // only makes sense after the decay and before the point has been
    // integrated into the histogram.
    shift_quantiles(val);
    int bucket_idx = insertValue(val);

    if (counts_[bucket_idx] < splitThreshold()) {
      return;
    }

    split(bucket_idx);
    merge();
  }

  size_t getNumBuckets() { return counts_.size(); }

  Bucket getBucketByIndex(size_t idx) const {
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

  int getBucketIndexByValue(double val) const {
    int bucket_idx = 0;
    for (; bucket_idx < ubounds_.size(); bucket_idx++) {
      if (val < ubounds_[bucket_idx]) {
        return bucket_idx;
      }
    }
    return ubounds_.size();
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

  double getMean() {
    double acc = counts_[0] * (ubounds_[0] + getMin()) / 2.0;
    int i = 1;
    for (; i < counts_.size() - 1; i++) {
      acc += counts_[i] * (ubounds_[i] + ubounds_[i - 1]) / 2.0;
    }
    acc += counts_[i] * (getMax() + ubounds_[i]) / 2.0;
    return acc / computeTotalCount();
  }

  // quantile is in [0, 1]
  double getQuantileEstimate(double quantile) const {
    double total_count = computeTotalCount();
    if (total_count == 0.0) {
      return 0.0;
    }

    double cdf = 0.0;
    double next_cdf = 0.0;
    int bucket_idx = 0;
    for (int i = 0; i < counts_.size(); i++) {
      next_cdf = cdf + counts_[i] / total_count;
      if (next_cdf > quantile) {
        bucket_idx = i;
        break;
      }
      cdf = next_cdf;
    }

    double frac = (quantile - cdf) / (next_cdf - cdf);

    double lower_bound;
    if (bucket_idx == 0) {
      lower_bound = getMin();
    } else {
      lower_bound = ubounds_[bucket_idx - 1];
    }

    double upper_bound;
    if (bucket_idx < ubounds_.size()) {
      upper_bound = ubounds_[bucket_idx];
    } else {
      upper_bound = getMax();
    }

    return frac * (upper_bound - lower_bound) + lower_bound;
  }

  void trackQuantiles(const std::vector<double> &quantiles) {
    std::copy(quantiles.begin(), quantiles.end(), back_inserter(quantiles_));
    quantile_locations_.clear();
    for (int i = 0; i < quantiles_.size(); i++) {
      quantile_locations_.push_back(getQuantileEstimate(quantiles_[i]));
    }
  }

  std::map<double, double> getTrackedQuantiles() const {
    std::map<double, double> qs;
    for (auto quantile : quantiles_) {
      qs[quantile] = getQuantileEstimate(quantile);
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

    if (!std::is_sorted(quantiles_.begin(), quantiles_.end())) {
      return false;
    }

    if (!std::is_sorted(quantile_locations_.begin(),
                        quantile_locations_.end())) {
      return false;
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
  const size_t max_num_buckets_;
  const double decay_rate_;

  uint64_t generation_;

  // ubounds_ records the upper bound of a bucket. There isn't an upper bound
  // for the last bucket, so the length of counts_ will generally be one greater
  // than the length of counts_.
  std::vector<double> ubounds_;
  std::vector<double> counts_;
  std::vector<double> quantiles_;
  std::vector<double> quantile_locations_;

  double splitThreshold() {
    double total_count = computeTotalCount();
    return 2 * total_count / getNumBuckets();
  }

  void decay() {
    for (int i = 0; i < counts_.size(); i++) {
      counts_[i] = counts_[i] * (1.0 - decay_rate_);
    }
  }

  // Adjust bounds and add.
  // Returns:
  //   The index of the bucket that the new value landed in.
  int insertValue(double val) {
    int bucket_idx = getBucketIndexByValue(val);
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
    for (int idx = 0; idx < quantiles_.size(); idx++) {
      double shift_remaining = quantiles_[idx];
      double location = quantile_locations_[idx];
      if (val < location) {
        shift_remaining = -shift_remaining;
      }

      int bucket_idx = getBucketIndexByValue(location);
      Bucket bucket = getBucketByIndex(bucket_idx);
      do {
        double target_location = val;
        if (bucket.count() != 0.0) {
          target_location = quantiles_[idx] + shift_remaining *
                                                  bucket.diameter() /
                                                  bucket.count();
        }

        if (in_range(val, bucket.min(), bucket.max()) &&
            in_range(val, location, target_location)) {
          // The inserted value is in between the current location and the
          // target location, which is in this bucket. Since the inserted value
          // has a point mass of 1, that's the farthest we can shift the value.
          quantile_locations_[idx] = val;
          break;
        }

        if (bucket_idx > 0 &&
            in_range(bucket.min(), location, target_location)) {
          // Shift to the lower bucket.
          shift_remaining = shift_remaining * (location - bucket.min()) /
                            (location - target_location);
          location = bucket.min();
          --bucket_idx;
          bucket = getBucketByIndex(bucket_idx);
        } else if (bucket_idx + 1 < ubounds_.size() &&
                   in_range(bucket.max(), location, target_location)) {
          // Shift to the upper bucket.
          shift_remaining = shift_remaining * (bucket.max() - location) /
                            (target_location - location);
          location = bucket.max();
          ++bucket_idx;
          bucket = getBucketByIndex(bucket_idx);
        } else {
          quantile_locations_[idx] = target_location;
          shift_remaining = 0.0;
        }
      } while (shift_remaining);
    }
  }

  // Shifts the quantile inside the bucket. Returns the amount that remains to
  // be shifted. If 0, the shift is done.
  double shift_quantile_in_bucket(int quantile_idx, double val,
                                  double shift_remaining) {
    return 0.0;
  }
};

} // namespace dhist
