// MIT License
//
// Copyright (c) 2021 Logan Evans
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

#include "DynamicHistogram.h"

namespace dyden {

DynamicHistogram::DynamicHistogram(const DynamicHistogramOpts& opts)
    : DensityMapBase(DescriptionOpts()
                         .set_type(MapType::HISTOGRAM)
                         .set_decay_rate(opts.decay_rate())
                         .set_refresh_interval(opts.refresh_interval())
                         .set_title(opts.title())
                         .set_labels({opts.label()})
                         .set_num_containers(opts.num_buckets())),
      generation_(0),
      refresh_generation_(0),
      total_count_(0.0),
      split_threshold_(0.0),
      insertion_buffer_(/*buffer_size=*/2 * opts.refresh_interval()) {
  ubounds_.resize(opts.num_buckets());
  ubounds_.back() = std::numeric_limits<double>::max();
  counts_.resize(opts.num_buckets());
  bucket_generation_.resize(opts.num_buckets());
  if (opts.register_with_server()) {
    registerWithServer();
  }
}

void DynamicHistogram::addValue(double val) {
  size_t unflushed = insertion_buffer_.addValue(val);
  if (unflushed >= description().refresh_interval()) {
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);
  }
}

Bucket DynamicHistogram::getBucketByIndex(size_t bx) {
  double min;
  if (bx > 0) {
    min = ubounds_[bx - 1];
  } else {
    min = getMinNoLock();
  }

  double max;
  if (bx + 1 < ubounds_.size()) {
    max = ubounds_[bx];
  } else {
    max = getMaxNoLock();
  }

  decay(bx);

  return Bucket(/*min=*/min, /*max=*/max, /*count=*/counts_[bx]);
}

double DynamicHistogram::computeTotalCount() {
  auto flush_it = insertion_buffer_.lockedIterator();
  flush(&flush_it);
  return total_count_;
}

double DynamicHistogram::getMin() {
  auto flush_it = insertion_buffer_.lockedIterator();
  flush(&flush_it);
  return getMinNoLock();
}

double DynamicHistogram::getMax() {
  auto flush_it = insertion_buffer_.lockedIterator();
  flush(&flush_it);
  return getMaxNoLock();
}

double DynamicHistogram::getMean() {
  auto flush_it = insertion_buffer_.lockedIterator();
  flush(&flush_it);

  double mean = 0.0;
  mean += counts_[0] * (getMinNoLock() + ubounds_[0]) / 2 / total_count_;

  size_t i = 1;
  for (; i < counts_.size() - 1; i++) {
    const double new_val = (ubounds_[i - 1] + ubounds_[i]) / 2;
    mean += counts_[i] * new_val / total_count_;
  }

  mean += counts_[i] * (ubounds_[i - 1] + getMaxNoLock()) / 2 / total_count_;
  return mean;
}

double DynamicHistogram::getQuantileEstimate(double quantile) {
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
    lower_bound = getMinNoLock();
  } else {
    lower_bound = ubounds_[bx - 1];
  }

  double upper_bound;
  if (bx + 1 < ubounds_.size()) {
    upper_bound = ubounds_[bx];
  } else {
    upper_bound = getMaxNoLock();
  }

  return frac * (upper_bound - lower_bound) + lower_bound;
}

double DynamicHistogram::getQuantileOfValue(double value) {
  auto flush_it = insertion_buffer_.lockedIterator();
  flush(&flush_it);

  if (value <= getMinNoLock()) {
    return getMinNoLock();
  }

  if (value >= getMaxNoLock()) {
    return getMaxNoLock();
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

std::string DynamicHistogram::debugString() {
  auto flush_it = insertion_buffer_.lockedIterator();
  flush(&flush_it);

  std::string s;
  s += "generation: " + std::to_string(total_count_) +
       "\n"
       "total_count: " +
       std::to_string(total_count_) + "\n";
  s += "  " + std::to_string(0) + " [" + std::to_string(getMinNoLock()) + ", " +
       std::to_string(ubounds_[0]) + "): " + std::to_string(counts_[0]) + "\n";
  size_t i = 1;
  for (; i < counts_.size() - 1; i++) {
    s += "  " + std::to_string(i) + " [" + std::to_string(ubounds_[i - 1]) +
         ", " + std::to_string(ubounds_[i]) +
         "): " + std::to_string(counts_[i]) + "\n";
  }
  s += "  " + std::to_string(i) + " [" + std::to_string(ubounds_[i - 1]) +
       ", " + std::to_string(getMaxNoLock()) + "): " + std::to_string(counts_[i]);
  return s;
}

std::string DynamicHistogram::json() {
  auto flush_it = insertion_buffer_.lockedIterator();
  flush(&flush_it);

  std::string s("{\n");
  if (!description().title().empty()) {
    s += "  \"title\": \"" + description().title() + "\",\n";
  }
  if (!description().labels().empty()) {
    s += "  \"label\": \"" + description().labels()[0] + "\",\n";
  }
  s += "  \"bounds\": [" + std::to_string(getMinNoLock()) + ", ";
  for (size_t i = 0; i + 1 < ubounds_.size(); i++) {
    s += std::to_string(ubounds_[i]) + ", ";
  }
  s += std::to_string(getMaxNoLock()) + "],\n  \"counts\": [";

  size_t i = 0;
  for (; i < counts_.size() - 1; i++) {
    s += std::to_string(counts_[i]) + ", ";
  }
  s += std::to_string(counts_[i]) + "]\n}";
  return s;
}

DensityMap DynamicHistogram::asProto() {
  DensityMap dm;
  toProto(&dm);
  return dm;
}

void DynamicHistogram::toProto(DensityMap* proto) {
  auto* dhist = proto->mutable_dynamic_histogram();

  auto* desc = dhist->mutable_description();
  description().toProto(desc);

  auto flush_it = insertion_buffer_.lockedIterator();
  flush(&flush_it);

  dhist->add_bounds(getMinNoLock());
  for (size_t i = 0; i + 1 < ubounds_.size(); i++) {
    dhist->add_bounds(ubounds_[i]);
  }
  dhist->add_bounds(getMaxNoLock());

  for (size_t i = 0; i < counts_.size(); i++) {
    dhist->add_counts(counts_[i]);
  }
}

void DynamicHistogram::registerWithServer() {
  DensityMapRegistry::getInstance().registerDensityMap(this);
}

double DynamicHistogram::splitThreshold() const { return split_threshold_; }

double DynamicHistogram::decay_rate() const {
  return description().decay_rate();
}

double DynamicHistogram::getMinNoLock() const {
  if (counts_[1] == 0) {
    return ubounds_[0];
  }
  return ubounds_[0] - (counts_[0] / counts_[1]) * (ubounds_[1] - ubounds_[0]);
}

double DynamicHistogram::getMaxNoLock() const {
  const size_t bx = ubounds_.size() - 2;
  if (counts_[bx] == 0.0) {
    // The last ubounds_ value is a fake value that allows std::upper_bound
    // to work.
    return ubounds_[ubounds_.size() - 2];
  }
  return ubounds_[bx] +
         (counts_[bx + 1] / counts_[bx]) * (ubounds_[bx] - ubounds_[bx - 1]);
}

double DynamicHistogram::getUpperBound(int i) const {
  if (i == -1) {
    return getMinNoLock();
  }
  if (static_cast<size_t>(i) == ubounds_.size()) {
    return getMaxNoLock();
  }
  return ubounds_[i];
}

void DynamicHistogram::flush(FlushIterator<double>* flush_it) {
  for (; *flush_it; ++(*flush_it)) {
    flushValue(**flush_it);
  }
  refresh();
}

void DynamicHistogram::flushValue(double val) {
  generation_++;

  size_t bx = insertValue(val);

  if (decay_rate() != 0.0) {
    total_count_ = total_count_ * (1.0 - decay_rate()) + 1.0;
  } else {
    total_count_ = generation_;
  }

  if (counts_[bx] < splitThreshold()) {
    return;
  }

  refresh();

  split(bx);
  merge();
}

void DynamicHistogram::decay(size_t bx) {
  uint64_t old_generation = bucket_generation_[bx];
  bucket_generation_[bx] = generation_;
  if (decay_rate() == 0.0) {
    return;
  }
  counts_[bx] *= description().decay_factor(generation_ - old_generation);
  bucket_generation_[bx] = generation_;
}

void DynamicHistogram::refresh() {
  if (refresh_generation_ == generation_) {
    return;
  }

  refresh_generation_ = generation_;

  total_count_ = 0.0;
  double min_count = std::numeric_limits<double>::max();
  double max_count = 0.0;
  for (size_t bx = 0; bx < counts_.size(); bx++) {
    decay(bx);
    double count = counts_[bx];
    total_count_ += count;
    if (count < min_count) {
      min_count = count;
    }
    if (count > max_count) {
      max_count = count;
    }
  }

  if (min_count * 4 < max_count) {
    split_threshold_ = max_count;
  } else {
    split_threshold_ = 2 * total_count_ / getNumBuckets();
  }

  if (static_cast<int32_t>(getNumBuckets()) != description().num_containers()) {
    while (static_cast<int32_t>(getNumBuckets()) <
           description().num_containers()) {
      size_t bx = 0;
      double highest_count = counts_[0];
      for (size_t i = 1; i < counts_.size(); i++) {
        if (counts_[i] > highest_count) {
          bx = i;
          highest_count = counts_[i];
        }
      }
      split(bx);
    }

    while (static_cast<int32_t>(getNumBuckets()) >
           description().num_containers()) {
      merge();
    }

    split_threshold_ = 2 * total_count_ / getNumBuckets();
  }
}

size_t DynamicHistogram::insertValue(double val) {
  size_t bx =
      std::distance(ubounds_.begin(),
                    std::upper_bound(ubounds_.begin(), ubounds_.end(), val));

  decay(bx);

  double count_bx = counts_[bx];
  if (bx > 0) {
    decay(bx - 1);
    double count_with_below = counts_[bx - 1] + count_bx;
    // Roundoff issues can happen, so make sure that the bound does *not* move
    // left.
    ubounds_[bx - 1] = std::max(
        (count_with_below * ubounds_[bx - 1] + val) / (count_with_below + 1.0),
        ubounds_[bx - 1]);
  }

  if (bx + 1 < counts_.size()) {
    decay(bx + 1);
    double count_with_above = count_bx + counts_[bx + 1];
    // Roundoff issues can happen, so make sure that the bound does *not* move
    // right.
    ubounds_[bx] = std::min(
        (count_with_above * ubounds_[bx] + val) / (count_with_above + 1.0),
        ubounds_[bx]);
  }

  counts_[bx] = count_bx + 1;

  return bx;
}

void DynamicHistogram::split(size_t bx) {
  double lower_bound;
  if (bx == 0) {
    lower_bound = getMinNoLock();
  } else {
    lower_bound = ubounds_[bx - 1];
  }

  if (bx + 1 == ubounds_.size()) {
    double upper_bound = getMaxNoLock();
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
  bucket_generation_.push_back(generation_);
}

void DynamicHistogram::merge() {
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
  bucket_generation_.pop_back();
}

}  // namespace dyden
