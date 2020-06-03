#ifndef DYNAMIC_HISTOGRAM_H
#define DYNAMIC_HISTOGRAM_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <vector>

namespace dhist {

bool in_range(double val, double a, double b);

class Bucket {
public:
  // Represents data in the half-open interval [min, max).
  Bucket(double min, double max, double count)
      : min_(min), max_(max), count_(count) {}

  double diameter() const { return max() - min(); }

  double min() const { return min_; }

  double max() const { return max_; }

  double count() const { return count_; }

private:
  double min_;
  double max_;
  double count_;
};


class DynamicHistogram {
public:
  DynamicHistogram(double decay_rate, size_t max_num_buckets)
      : use_decay_(decay_rate > 0.0), decay_rate_(decay_rate),
        max_num_buckets_(max_num_buckets), generation_(0) {
    ubounds_.resize(max_num_buckets_ - 1);
    counts_.resize(max_num_buckets_);
    bucket_generation_.resize(max_num_buckets_);
  }

  virtual ~DynamicHistogram() {}

  void addValue(double val) {
    generation_++;
    //decay();
    // shift_quantiles will treat the new value as a point with mass 1 which
    // only makes sense after the decay and before the point has been
    // integrated into the histogram.
    //shift_quantiles(val);
    int bx = insertValue(val);

    if (counts_[bx] < splitThreshold()) {
      return;
    }

    decayAll();
    split(bx);
    merge();
  }

  void addRepeatedValue(double val, uint64_t nValues) {
    for (int i = 0; i < nValues; i++) {
      addValue(val);
    }
  }

  void clear() {
  }

  // TODO(lpe)
  void merge(const DynamicHistogram &other) { assert(false); }

  // TODO(lpe)
  void copy(const DynamicHistogram &other) { assert(false); }

  size_t getNumBuckets() { return counts_.size(); }

  Bucket getBucketByIndex(size_t idx) const {
    return Bucket(0, 0, 0);
  }

  int getBucketIndexByValue(double val) const {
    return -1;
  }

  double getMin() const {
    assert(bucket_generation_[0] == bucket_generation_[1]);
    if (counts_[1] == 0) {
      return ubounds_[0];
    }
    return ubounds_[0] -
           (counts_[0] / counts_[1]) * (ubounds_[1] - ubounds_[0]);
  }

  double getMax() const {
    const int bx = ubounds_.size() - 1;
    assert(bucket_generation_[bx] == bucket_generation_[bx - 1]);
    if (counts_[bx] == 0.0) {
      return ubounds_.back();
    }
    return ubounds_[bx] + (counts_[bx + 1] / counts_[bx]) *
                               (ubounds_[bx] - ubounds_[bx - 1]);
  }

  double computeTotalCount() const {
    if (!use_decay_) {
      return generation_;
    }
    printf("> computeTotalCount: 1.0 - pow(%lf, %ld)) / (1 - %lf) == %lf\n",
           decay_rate_, generation_, decay_rate_,
           (1.0 - pow(1 - decay_rate_, generation_)) / (1 - decay_rate_));
    return (1.0 - pow(1 - decay_rate_, generation_)) / (1 - decay_rate_);
  }

  // quantile is in [0, 1]
  double getQuantileEstimate(double quantile) {
    decayAll();

    double total_count = computeTotalCount();
    if (total_count == 0.0) {
      return 0.0;
    }

    double cdf = 0.0;
    double next_cdf = 0.0;
    int bx = 0;
    for (int i = 0; i < counts_.size(); i++) {
      next_cdf = cdf + counts_[i] / total_count;
      if (next_cdf > quantile) {
        bx = i;
        break;
      }
      cdf = next_cdf;
    }

    double frac = (quantile - cdf) / (next_cdf - cdf);

    double lower_bound;
    if (bx == 0) {
      lower_bound = getMin();
    } else {
      lower_bound = ubounds_[bx - 1];
    }

    double upper_bound;
    if (bx < ubounds_.size()) {
      upper_bound = ubounds_[bx];
    } else {
      upper_bound = getMax();
    }

    return frac * (upper_bound - lower_bound) + lower_bound;
  }

  void trackQuantiles(const std::vector<double> &quantiles) {
  }

  std::map<double, double> getTrackedQuantiles() const {
    return std::map<double, double>();
  }

  std::string debugString() {
    decayAll();
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

  std::string json() {
    decayAll();
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
  const bool use_decay_;
  const double decay_rate_;
  const size_t max_num_buckets_;

  uint64_t generation_;

  // ubounds_ records the upper bound of a bucket. There isn't an upper bound
  // for the last bucket, so the length of counts_ will generally be one greater
  // than the length of counts_.
  std::vector<double> ubounds_;
  std::vector<double> counts_;
  std::vector<uint64_t> bucket_generation_;
  std::vector<double> quantiles_;
  std::vector<double> quantile_locations_;

  double splitThreshold() {
    double total_count = computeTotalCount();
    return 2 * total_count / getNumBuckets();
  }

  double computeDecay(double original, uint64_t generations) {
    if (!use_decay_) {
      return original;
    }
    printf("> computeDecay(%lf, %ld) == %lf\n", original, generations,
           original * pow(decay_rate_, generations));
    return original * pow(decay_rate_, generations);
  }

  void decay(int bx) {
    counts_[bx] =
        computeDecay(counts_[bx], generation_ - bucket_generation_[bx]);
    bucket_generation_[bx] = generation_;
  }

  void decayAll() {
    for (int bx = 0; bx < counts_.size(); bx++) {
      decay(bx);
    }
  }

  // Adjust bounds and add.
  // Returns:
  //   The index of the bucket that the new value landed in.
  int insertValue(double val) {
    int bx =
        std::distance(ubounds_.begin(),
                      std::upper_bound(ubounds_.begin(), ubounds_.end(), val));

    decay(bx);

    if (bx > 0) {
      decay(bx - 1);
      double count_with_below = counts_[bx - 1] + counts_[bx];
      ubounds_[bx - 1] = (count_with_below * ubounds_[bx - 1] + val) /
                         (count_with_below + 1.0);
    }

    if (bx < ubounds_.size()) {
      decay(bx + 1);
      double count_with_above = counts_[bx] + counts_[bx + 1];
      ubounds_[bx] =
          (count_with_above * ubounds_[bx] + val) / (count_with_above + 1.0);
    }

    counts_[bx] += 1;

    return bx;
  }

  void split(int bx) {
    double lower_bound;
    if (bx == 0) {
      lower_bound = getMin();
    } else {
      lower_bound = ubounds_[bx - 1];
    }

    if (bx == ubounds_.size()) {
      double upper_bound = getMax();
      ubounds_.push_back((lower_bound + upper_bound) / 2.0);
      counts_[bx] /= 2.0;
      counts_.push_back(counts_.back());
    } else {
      double upper_bound = ubounds_[bx];

      ubounds_.push_back(0.0);
      for (int i = ubounds_.size() - 1; i > bx; i--) {
        ubounds_[i] = ubounds_[i - 1];
      }

      ubounds_[bx] = (lower_bound + upper_bound) / 2.0;

      counts_.push_back(0.0);
      for (int i = counts_.size() - 1; i > bx; i--) {
        counts_[i] = counts_[i - 1];
      }
      counts_[bx] /= 2.0;
      counts_[bx + 1] = counts_[bx];
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

      int bx = getBucketIndexByValue(location);
      Bucket bucket = getBucketByIndex(bx);
      do {
        double target_location = val;
        if (bucket.count() != 0.0) {
          target_location = quantiles_[idx] + shift_remaining *
                                                  bucket.diameter() /
                                                  bucket.count();
        }

        //printf("idx: %d, quantile: %lf, location: %lf, target_location: %lf, "
        //       "bucket.min(): %lf, bucket.max(): %lf\n",
        //       idx, quantiles_[idx], location, target_location, bucket.min(),
        //       bucket.max());

        if (in_range(val, bucket.min(), bucket.max()) &&
            in_range(val, location, target_location)) {
          // The inserted value is in between the current location and the
          // target location, which is in this bucket. Since the inserted value
          // has a point mass of 1, that's the farthest we can shift the value.
          quantile_locations_[idx] = val;
          break;
        }

        if (bx > 0 &&
            in_range(bucket.min(), location, target_location)) {
          // Shift to the lower bucket.
          shift_remaining = shift_remaining * (location - bucket.min()) /
                            (location - target_location);
          location = bucket.min();
          --bx;
          bucket = getBucketByIndex(bx);
        } else if (bx + 1 < ubounds_.size() &&
                   in_range(bucket.max(), location, target_location)) {
          // Shift to the upper bucket.
          shift_remaining = shift_remaining * (bucket.max() - location) /
                            (target_location - location);
          location = bucket.max();
          ++bx;
          bucket = getBucketByIndex(bx);
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

}

#endif  // DYNAMIC_HISTOGRAM_H
