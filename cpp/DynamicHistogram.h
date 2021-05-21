// MIT License
//
// Copyright (c) 2020 Logan Evans
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <google/protobuf/util/time_util.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <map>
#include <vector>

#include "DynamicDensity.pb.h"
#include "cpp/InsertionBuffer.h"

namespace dhist {

using ::dynamic_density::DensityMap;
using ::google::protobuf::Timestamp;

bool in_range(double val, double a, double b) {
  return ((val <= a) ^ (val <= b)) || val == a || val == b;
}

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
  DynamicHistogram(size_t num_buckets, double decay_rate = 0.0,
                   size_t refresh_interval = 512)
      : decay_rate_(decay_rate),
        refresh_interval_(refresh_interval),
        generation_(0),
        refresh_generation_(0),
        total_count_(0.0),
        insertion_buffer_(/*buffer_size=*/2 * refresh_interval) {
    ubounds_.resize(num_buckets);
    ubounds_.back() = std::numeric_limits<double>::max();
    counts_.resize(num_buckets);
    bucket_generation_.resize(num_buckets);

    if (decay_rate_ != 0.0) {
      decay_factors_.resize(refresh_interval_);
      double decay = 1.0;
      for (size_t i = 0; i < refresh_interval_; i++) {
        decay_factors_[i] = decay;
        decay *= 1.0 - decay_rate_;
      }
    }
  }

  virtual ~DynamicHistogram() {}

  void addValue(double val) {
    size_t unflushed = insertion_buffer_.addValue(val);
    if (unflushed >= refresh_interval_) {
      auto flush_it = insertion_buffer_.lockedIterator();
      flush(&flush_it);
    }
  }

  size_t getNumBuckets() const { return counts_.size(); }

  Bucket getBucketByIndex(size_t bx) {
    double min;
    if (bx > 0) {
      min = ubounds_[bx - 1];
    } else {
      min = getMin();
    }

    double max;
    if (bx + 1 < ubounds_.size()) {
      max = ubounds_[bx];
    } else {
      max = getMax();
    }

    decay(bx);

    return Bucket(/*min=*/min, /*max=*/max, /*count=*/counts_[bx]);
  }

  int getBucketIndexByValue(double val) const { return -1; }

  double getMin() const {
    assert(bucket_generation_[0] == bucket_generation_[1]);
    if (counts_[1] == 0) {
      return ubounds_[0];
    }
    return ubounds_[0] -
           (counts_[0] / counts_[1]) * (ubounds_[1] - ubounds_[0]);
  }

  double getMax() const {
    const size_t bx = ubounds_.size() - 2;
    assert(bucket_generation_[bx] == bucket_generation_[bx - 1]);
    if (counts_[bx] == 0.0) {
      // The last ubounds_ value is a fake value that allows std::upper_bound
      // to work.
      return ubounds_[ubounds_.size() - 2];
    }
    return ubounds_[bx] +
           (counts_[bx + 1] / counts_[bx]) * (ubounds_[bx] - ubounds_[bx - 1]);
  }

  double getUpperBound(int i) const {
    if (i == -1) {
      return getMin();
    }
    if (static_cast<size_t>(i) == ubounds_.size()) {
      return getMax();
    }
    return ubounds_[i];
  }

  double computeTotalCount() {
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);
    return total_count_;
  }

  double getMean() {
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);

    double mean = 0.0;
    mean += counts_[0] * (getMin() + ubounds_[0]) / 2 / total_count_;

    size_t i = 1;
    for (; i < counts_.size() - 1; i++) {
      const double new_val = (ubounds_[i - 1] + ubounds_[i]) / 2;
      mean += counts_[i] * new_val / total_count_;
    }

    mean += counts_[i] * (ubounds_[i - 1] + getMax()) / 2 / total_count_;
    return mean;
  }

  // quantile is in [0, 1]
  double getQuantileEstimate(double quantile) {
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);

    if (total_count_ == 0.0) {
      return 0.0;
    }

    double cdf = 0.0;
    double next_cdf = 0.0;
    size_t bx = 0;
    for (size_t i = 0; i < counts_.size(); i++) {
      next_cdf = cdf + counts_[i] / total_count_;
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
    if (bx + 1 < ubounds_.size()) {
      upper_bound = ubounds_[bx];
    } else {
      upper_bound = getMax();
    }

    return frac * (upper_bound - lower_bound) + lower_bound;
  }

  double getQuantileOfValue(double value) {
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);

    if (value <= getMin()) {
      return getMin();
    }

    if (value >= getMax()) {
      return getMax();
    }

    double count = 0.0;
    for (size_t i = 0; i < counts_.size(); i++) {
      if (ubounds_[i] >= value) {
        count += counts_[i];
      } else {
        count += (getUpperBound(i) - value) /
                 (getUpperBound(i) - getUpperBound(static_cast<int>(i) - 1));
        break;
      }
    }

    return count / computeTotalCount();
  }

  void trackQuantiles(const std::vector<double> &quantiles) {}

  std::map<double, double> getTrackedQuantiles() const {
    return std::map<double, double>();
  }

  std::string debugString() {
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);

    std::string s;
    s += "generation: " + std::to_string(total_count_) +
         "\n"
         "total_count: " +
         std::to_string(total_count_) + "\n";
    s += "  " + std::to_string(0) + " [" + std::to_string(getMin()) + ", " +
         std::to_string(ubounds_[0]) + "): " + std::to_string(counts_[0]) +
         "\n";
    size_t i = 1;
    for (; i < counts_.size() - 1; i++) {
      s += "  " + std::to_string(i) + " [" + std::to_string(ubounds_[i - 1]) +
           ", " + std::to_string(ubounds_[i]) +
           "): " + std::to_string(counts_[i]) + "\n";
    }
    s += "  " + std::to_string(i) + " [" + std::to_string(ubounds_[i - 1]) +
         ", " + std::to_string(getMax()) + "): " + std::to_string(counts_[i]);
    return s;
  }

  std::string json(std::string title = "", std::string label = "") {
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);

    std::string s("{\n");
    if (!title.empty()) {
      s += "  \"title\": \"" + title + "\",\n";
    }
    if (!label.empty()) {
      s += "  \"label\": \"" + label + "\",\n";
    }
    s += "  \"bounds\": [" + std::to_string(getMin()) + ", ";
    for (size_t i = 0; i + 1 < ubounds_.size(); i++) {
      s += std::to_string(ubounds_[i]) + ", ";
    }
    s += std::to_string(getMax()) + "],\n  \"counts\": [";

    size_t i = 0;
    for (; i < counts_.size() - 1; i++) {
      s += std::to_string(counts_[i]) + ", ";
    }
    s += std::to_string(counts_[i]) + "]\n}";
    return s;
  }

  DensityMap toProto(std::string title = "", std::string label = "") {
    DensityMap dm;
    auto *dhist = dm.mutable_dynamic_histogram();

    if (!title.empty()) {
      dhist->set_title(title);
    }

    if (!label.empty()) {
      dhist->set_label(label);
    }

    struct timeval tv;
    gettimeofday(&tv, NULL);
    dhist->mutable_timestamp()->set_seconds(tv.tv_sec);
    dhist->mutable_timestamp()->set_nanos(tv.tv_usec * 1000);

    dhist->set_decay_rate(decay_rate_);

    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);

    dhist->add_bounds(getMin());
    for (size_t i = 0; i + 1 < ubounds_.size(); i++) {
      dhist->add_bounds(ubounds_[i]);
    }
    dhist->add_bounds(getMax());

    for (size_t i = 0; i < counts_.size(); i++) {
      dhist->add_counts(counts_[i]);
    }

    return dm;
  }

 protected:
  const double decay_rate_;
  const size_t refresh_interval_;

  uint64_t generation_;
  uint64_t refresh_generation_;
  double total_count_;

  InsertionBuffer<double> insertion_buffer_;

  // ubounds_ records the upper bound of a bucket. There isn't an upper bound
  // for the last bucket, so the length of counts_ will generally be one greater
  // than the length of counts_.
  std::vector<double> ubounds_;
  std::vector<double> counts_;
  std::vector<uint64_t> bucket_generation_;
  std::vector<double> quantiles_;
  std::vector<double> quantile_locations_;
  std::vector<double> decay_factors_;

  void flush(FlushIterator<double>* flush_it) {
    for (; *flush_it; ++(*flush_it)) {
      flushValue(**flush_it);
    }
    refresh();
  }

  void flushValue(double val) {
    generation_++;

    // shift_quantiles will treat the new value as a point with mass 1 which
    // only makes sense after the decay and before the point has been
    // integrated into the histogram.
    // shift_quantiles(val);
    size_t bx = insertValue(val);

    if (decay_rate_ != 0.0) {
      total_count_ = total_count_ * (1.0 - decay_rate_) + 1.0;
    } else {
      total_count_ = generation_;
    }

    if (generation_ - refresh_generation_ + 1 == refresh_interval_) {
      refresh();
    }

    if (counts_[bx] < splitThreshold()) {
      return;
    }

    refresh();
    if (decay_rate_ != 0.0) {
      for (size_t i = 0; i < bucket_generation_.size(); i++) {
        assert(bucket_generation_[i] == generation_);
      }
    }

    split(bx);
    merge();
  }

  double splitThreshold() const { return 2 * total_count_ / getNumBuckets(); }

  double countAfterDecay(double original, uint64_t generations) {
    if (decay_rate_ == 0.0) {
      return original;
    }
    // This invariant should be maintained by addValue calling refresh()
    // periodically.
    assert(generations < decay_factors_.size());
    return original * decay_factors_[generations];
  }

  void decay(size_t bx) {
    if (decay_rate_ == 0.0) {
      return;
    }
    counts_[bx] =
        countAfterDecay(counts_[bx], generation_ - bucket_generation_[bx]);
    bucket_generation_[bx] = generation_;
  }

  void refresh() {
    if (refresh_generation_ == generation_) {
      return;
    }

    total_count_ = 0.0;
    for (size_t bx = 0; bx < counts_.size(); bx++) {
      decay(bx);
      total_count_ += counts_[bx];
    }
    refresh_generation_ = generation_;
  }

  // Adjust bounds and add.
  // Returns:
  //   The index of the bucket that the new value landed in.
  size_t insertValue(double val) {
    size_t bx =
        std::distance(ubounds_.begin(),
                      std::upper_bound(ubounds_.begin(), ubounds_.end(), val));

    decay(bx);

    double count_bx = counts_[bx];
    if (bx > 0) {
      decay(bx - 1);
      double count_with_below = counts_[bx - 1] + count_bx;
      ubounds_[bx - 1] = (count_with_below * ubounds_[bx - 1] + val) /
                         (count_with_below + 1.0);
    }

    if (bx + 1 < counts_.size()) {
      decay(bx + 1);
      double count_with_above = count_bx + counts_[bx + 1];
      ubounds_[bx] =
          (count_with_above * ubounds_[bx] + val) / (count_with_above + 1.0);
    }

    counts_[bx] = count_bx + 1;

    return bx;
  }

  void split(size_t bx) {
    double lower_bound;
    if (bx == 0) {
      lower_bound = getMin();
    } else {
      lower_bound = ubounds_[bx - 1];
    }

    if (bx + 1 == ubounds_.size()) {
      double upper_bound = getMax();
      ubounds_.back() = (lower_bound + upper_bound) / 2.0;
      ubounds_.push_back(std::numeric_limits<double>::max());
      counts_[bx] /= 2.0;
      counts_.push_back(counts_.back());
    } else {
      double upper_bound = ubounds_[bx];

      ubounds_.push_back(0.0);
      for (size_t i = ubounds_.size() - 1; i > bx; i--) {
        ubounds_[i] = ubounds_[i - 1];
      }

      ubounds_[bx] = (lower_bound + upper_bound) / 2.0;

      counts_.push_back(0.0);
      for (size_t i = counts_.size() - 1; i > bx; i--) {
        counts_[i] = counts_[i - 1];
      }
      counts_[bx] /= 2.0;
      counts_[bx + 1] = counts_[bx];
    }
  }

  void merge() {
    int merge_idx = 0;
    double merged_count = counts_[0] + counts_[1];
    for (size_t i = 1; i < counts_.size() - 1; i++) {
      double pos_count = counts_[i] + counts_[i + 1];
      if (pos_count < merged_count) {
        merge_idx = i;
        merged_count = pos_count;
      }
    }

    for (size_t i = merge_idx; i < ubounds_.size() - 1; i++) {
      ubounds_[i] = ubounds_[i + 1];
    }
    ubounds_.pop_back();

    counts_[merge_idx] = merged_count;
    for (size_t i = merge_idx + 1; i < counts_.size() - 1; i++) {
      counts_[i] = counts_[i + 1];
    }
    counts_.pop_back();
  }

  void shift_quantiles(double val) {
    for (size_t idx = 0; idx < quantiles_.size(); idx++) {
      double shift_remaining = quantiles_[idx];
      double location = quantile_locations_[idx];
      if (val < location) {
        shift_remaining = -shift_remaining;
      }

      size_t bx = getBucketIndexByValue(location);
      Bucket bucket = getBucketByIndex(bx);
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

        if (bx > 0 && in_range(bucket.min(), location, target_location)) {
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

}  // namespace dhist
