#ifndef DYNAMIC_HISTOGRAM_H
#define DYNAMIC_HISTOGRAM_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <vector>

namespace dhist {

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

template <bool kUseDecay = true, bool kThreadsafe = true>
class DynamicHistogram {
public:
  DynamicHistogram(size_t max_num_buckets, double decay_rate = 0.0,
                   int refresh_interval = 512)
      : max_num_buckets_(max_num_buckets), decay_rate_(decay_rate),
        refresh_interval_(refresh_interval), generation_(0),
        refresh_generation_(0), total_count_(0.0), insert_queue_begin_(0),
        insert_queue_end_(0), insert_queue_to_flush_(0) {
    assert(kUseDecay == (decay_rate != 0.0));
    ubounds_.resize(max_num_buckets_);
    ubounds_.back() = std::numeric_limits<double>::max();
    counts_.resize(max_num_buckets_);
    bucket_generation_.resize(max_num_buckets_);
    insert_queue_.resize(refresh_interval_);

    if (kUseDecay) {
      decay_factors_.resize(refresh_interval_);
      double decay = 1.0;
      for (int i = 0; i < refresh_interval_; i++) {
        decay_factors_[i] = decay;
        decay *= 1.0 - decay_rate_;
      }
    }
  }

  virtual ~DynamicHistogram() {}

  void addValue(double val) { addValueImpl<kThreadsafe>(val); }

  void addRepeatedValue(double val, uint64_t nValues) {
    for (int i = 0; i < nValues; i++) {
      addValue(val);
    }
  }

  // TODO(lpe)
  void clear() { assert(false); }

  // TODO(lpe)
  void merge(const DynamicHistogram &other) { assert(false); }

  // TODO(lpe)
  void copy(const DynamicHistogram &other) { assert(false); }

  size_t getNumBuckets() { return counts_.size(); }

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
    const int bx = ubounds_.size() - 2;
    assert(bucket_generation_[bx] == bucket_generation_[bx - 1]);
    if (counts_[bx] == 0.0) {
      // The last ubounds_ value is a fake value that allows std::upper_bound
      // to work.
      return ubounds_[ubounds_.size() - 2];
    }
    return ubounds_[bx] +
           (counts_[bx + 1] / counts_[bx]) * (ubounds_[bx] - ubounds_[bx - 1]);
  }

  double computeTotalCount() {
    int to_flush = reserve_flush_items();
    std::unique_ptr<std::scoped_lock<std::mutex>> lp;
    if (kThreadsafe) {
      lp.reset(new std::scoped_lock(flush_mu_));
    }
    flush<kThreadsafe>(to_flush);
    return total_count_;
  }

  double getMean() {
    // TODO(lpe): This snippet repeats a few times. Why?
    int to_flush = reserve_flush_items();
    std::unique_ptr<std::scoped_lock<std::mutex>> lp;
    if (kThreadsafe) {
      lp.reset(new std::scoped_lock(flush_mu_));
    }
    flush<kThreadsafe>(to_flush);

    double acc = counts_[0] * (ubounds_[0] + getMin());
    int i = 1;
    for (; i < counts_.size() - 1; i++) {
      acc += counts_[i] * (ubounds_[i] + ubounds_[i - 1]);
    }
    acc += counts_[i] * (getMax() + ubounds_[i]);
    return acc / total_count_ / 2.0;
  }

  // quantile is in [0, 1]
  double getQuantileEstimate(double quantile) {
    int to_flush = reserve_flush_items();
    std::unique_ptr<std::scoped_lock<std::mutex>> lp;
    if (kThreadsafe) {
      lp.reset(new std::scoped_lock(flush_mu_));
    }
    flush<kThreadsafe>(to_flush);

    if (total_count_ == 0.0) {
      return 0.0;
    }

    double cdf = 0.0;
    double next_cdf = 0.0;
    int bx = 0;
    for (int i = 0; i < counts_.size(); i++) {
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

  void trackQuantiles(const std::vector<double> &quantiles) {}

  std::map<double, double> getTrackedQuantiles() const {
    return std::map<double, double>();
  }

  std::string debugString() {
    int to_flush = reserve_flush_items();
    std::unique_ptr<std::scoped_lock<std::mutex>> lp;
    if (kThreadsafe) {
      lp.reset(new std::scoped_lock(flush_mu_));
    }
    flush<kThreadsafe>(to_flush);

    std::string s;
    s += "generation: " + std::to_string(total_count_) + "\n"
         "total_count: " + std::to_string(computeTotalCount()) + "\n";
    s += "  " + std::to_string(0) + " [" + std::to_string(getMin()) + ", " +
         std::to_string(ubounds_[0]) + "): " + std::to_string(counts_[0]) +
         "\n";
    int i = 1;
    for (; i < ubounds_.size() - 1; i++) {
      s += "  " + std::to_string(i) + " [" + std::to_string(ubounds_[i - 1]) +
           ", " + std::to_string(ubounds_[i]) +
           "): " + std::to_string(counts_[i]) + "\n";
    }
    s += "  " + std::to_string(i) + " [" + std::to_string(ubounds_[i - 1]) +
         ", " + std::to_string(getMax()) +
         "): " + std::to_string(counts_[i]);
    return s;
  }

  std::string json(std::string title = "", std::string label = "") {
    int to_flush = reserve_flush_items();
    std::unique_ptr<std::scoped_lock<std::mutex>> lp;
    if (kThreadsafe) {
      lp.reset(new std::scoped_lock(flush_mu_));
    }
    flush<kThreadsafe>(to_flush);

    std::string s("{\n");
    if (!title.empty()) {
      s += "  \"title\": \"" + title + "\",\n";
    }
    if (!label.empty()) {
      s += "  \"label\": \"" + label + "\",\n";
    }
    s += "  \"bounds\": [" + std::to_string(getMin()) + ", ";
    for (int i = 0; i + 1 < ubounds_.size(); i++) {
      s += std::to_string(ubounds_[i]) + ", ";
    }
    s += std::to_string(getMax()) + "],\n  \"counts\": [";

    int i = 0;
    for (; i < counts_.size() - 1; i++) {
      s += std::to_string(counts_[i]) + ", ";
    }
    s += std::to_string(counts_[i]) + "]\n}";
    return s;
  }

protected:
  const size_t max_num_buckets_;
  const double decay_rate_;
  const int refresh_interval_;

  uint64_t generation_;
  uint64_t refresh_generation_;
  double total_count_;
  int insert_queue_begin_;
  int insert_queue_end_;
  int insert_queue_to_flush_;

  std::mutex insert_mu_;
  std::mutex flush_mu_;

  // ubounds_ records the upper bound of a bucket. There isn't an upper bound
  // for the last bucket, so the length of counts_ will generally be one greater
  // than the length of counts_.
  std::vector<double> ubounds_;
  std::vector<double> counts_;
  std::vector<uint64_t> bucket_generation_;
  std::vector<double> insert_queue_;
  std::vector<double> quantiles_;
  std::vector<double> quantile_locations_;
  std::vector<double> decay_factors_;

  template <bool threadsafe> int insert_queue_count() const;

  template <> int insert_queue_count</*threadsafe=*/true>() const {
    int size = insert_queue_end_ - insert_queue_begin_;
    if (size < 0) {
      size += insert_queue_.size();
    }
    return size;
  }

  template <> int insert_queue_count</*threadsafe=*/false>() const { return 0; }

  template <bool threadsafe> void flush(int items);

  template <> void flush</*threadsafe=*/true>(int items) {
    const int size = insert_queue_.size();
    int qx = insert_queue_begin_;
    for (int i = 0; i < items; i++) {
      // Use the unlocked version of addValueImpl now that we're under a lock.
      addValueImpl</*threadsafe=*/false>(insert_queue_[qx]);
      qx++;
      if (qx == size) {
        qx = 0;
      }
    }
    insert_queue_begin_ = qx;
    refresh();
  }

  template <> void flush</*threadsafe=*/false>(int items) { refresh(); }

  template <bool threadsafe> void addValueImpl(double val);

  template <> void addValueImpl</*threadsafe=*/true>(double val) {
    std::scoped_lock insert_lock(insert_mu_);

    while (kThreadsafe &&
           (2 * insert_queue_to_flush_ >= insert_queue_.size())) {
      int items = insert_queue_to_flush_;
      insert_queue_to_flush_ = 0;

      insert_mu_.unlock();
      {
        std::scoped_lock flush_lock(flush_mu_);
        flush</*kThreadsafe=*/true>(items);
      }
      insert_mu_.lock();
    }

    insert_queue_[insert_queue_end_] = val;
    insert_queue_end_++;
    insert_queue_to_flush_++;

    if (insert_queue_end_ == insert_queue_.size()) {
      insert_queue_end_ = 0;
    }
  }

  template <> void addValueImpl</*threadsafe=*/false>(double val) {
    generation_++;

    // decay();
    // shift_quantiles will treat the new value as a point with mass 1 which
    // only makes sense after the decay and before the point has been
    // integrated into the histogram.
    // shift_quantiles(val);
    int bx = insertValue(val);

    if (kUseDecay) {
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
    if (kUseDecay) {
      for (int i = 0; i < bucket_generation_.size(); i++) {
        assert(bucket_generation_[i] == generation_);
      }
    }

    split(bx);
    merge();
  }

  int reserve_flush_items() {
    // TODO: Could this just use an atomic for insert_queue_to_flush_?
    std::scoped_lock l(insert_mu_);
    int to_flush = insert_queue_to_flush_;
    insert_queue_to_flush_ = 0;
    return to_flush;
  }

  double splitThreshold() { return 2 * total_count_ / getNumBuckets(); }

  double countAfterDecay(double original, uint64_t generations) {
    if (!kUseDecay) {
      return original;
    }
    // This invariant should be maintained by addValue calling refresh()
    // periodically.
    assert(generations < decay_factors_.size());
    return original * decay_factors_[generations];
  }

  void decay(int bx) {
    if (!kUseDecay) {
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
    for (int bx = 0; bx < counts_.size(); bx++) {
      decay(bx);
      total_count_ += counts_[bx];
    }
    refresh_generation_ = generation_;
  }

  // Adjust bounds and add.
  // Returns:
  //   The index of the bucket that the new value landed in.
  int insertValue(double val) {
    int bx =
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

  void split(int bx) {
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

} // namespace dhist

#endif // DYNAMIC_HISTOGRAM_H
