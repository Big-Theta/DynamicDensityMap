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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <vector>

#include "DensityMap.pb.h"
#include "cpp/InsertionBuffer.h"

namespace dhist {

using ::dynamic_histogram::DensityMap;
using ::google::protobuf::Timestamp;

class Kernel {
 public:
  Kernel()
      : mean_(0.0),
        weighted_sum_(0.0),
        weighted_sum_squares_(0.0),
        S_(0.0),
        generation_(0) {}

  ~Kernel() = default;

  Kernel(const Kernel& other) {
    mean_ = other.mean_;
    weighted_sum_ = other.weighted_sum_;
    weighted_sum_squares_ = other.weighted_sum_squares_;
    S_ = other.S_;
    generation_ = other.generation_;
  }

  Kernel& operator=(const Kernel& other) {
    mean_ = other.mean_;
    weighted_sum_ = other.weighted_sum_;
    weighted_sum_squares_ = other.weighted_sum_squares_;
    S_ = other.S_;
    generation_ = other.generation_;
    return *this;
  }

  std::string debugString() const {
    std::string s = "Kernel{mean: " + std::to_string(mean()) +
                    ", variance: " + std::to_string(variance()) +
                    ", count: " + std::to_string(count()) + "}";
    return s;
  }

  void addValue(double value, double weight) {
    weighted_sum_ += weight;
    weighted_sum_squares_ += weight * weight;
    double old_mean = mean_;
    mean_ += (weight / weighted_sum_) * (value - old_mean);
    S_ += weight * (value - old_mean) * (value - mean_);
  }

  void decay(double factor, uint64_t new_generation) {
    generation_ = new_generation;
    if (factor == 1.0) {
      return;
    }
    weighted_sum_ *= factor;
    weighted_sum_squares_ *= factor;
    S_ *= factor;
  }

  double mean() const {
    return mean_;
  }

  double standard_deviation() const {
    return sqrt(variance());
  }

  double variance() const {
    return S_ / weighted_sum_;
  }

  double count() const {
    return weighted_sum_;
  }

  uint64_t generation() const {
    return generation_;
  }

  double cdf(double value) const {
    double z = (value - mean()) / standard_deviation();
    return 0.5 * erfc(-z * M_SQRT1_2);
  }

  static bool lessThan(double value, const Kernel& b) {
    return value < b.mean_;
  }

  static Kernel merge(const Kernel& a, const Kernel& b) {
    Kernel kernel;
    kernel.weighted_sum_ = a.weighted_sum_ + b.weighted_sum_;
    kernel.weighted_sum_squares_ =
        a.weighted_sum_squares_ + b.weighted_sum_squares_;
    assert(a.generation_ == b.generation_);
    kernel.generation_ = a.generation_;
    if (a.weighted_sum_ <= 1.0 || b.weighted_sum_ <= 1.0) {
      kernel.mean_ = (a.mean_ + b.mean_) / 2.0;
      kernel.S_ = 0.0;
    } else {
      kernel.mean_ = (a.mean_ * a.weighted_sum_ + b.mean_ * b.weighted_sum_) /
                     kernel.weighted_sum_;
      kernel.S_ = (a.weighted_sum_ - 1.0) * a.variance() +
                  (b.weighted_sum_ - 1.0) * b.variance();
    }
    return kernel;
  }

 private:
  double mean_;
  double weighted_sum_;
  double weighted_sum_squares_;
  double S_;
  uint64_t generation_;
};

class DynamicKDE {
 public:
  DynamicKDE(size_t num_kernels, double decay_rate = 0.0,
             int refresh_interval = 512)
      : decay_rate_(decay_rate),
        refresh_interval_(refresh_interval),
        generation_(0),
        refresh_generation_(0),
        total_count_(0.0),
        insertion_buffer_(/*buffer_size=*/2 * refresh_interval) {
    kernels_.reserve(num_kernels + 2);
    kernels_.resize(num_kernels);

    if (decay_rate_ != 0.0) {
      decay_factors_.resize(refresh_interval_);
      double decay = 1.0;
      for (int i = 0; i < refresh_interval_; i++) {
        decay_factors_[i] = decay;
        decay *= 1.0 - decay_rate_;
      }
    }
  }

  void addValue(double val) {
    size_t unflushed = insertion_buffer_.addValue(val);
    if (unflushed >= refresh_interval_) {
      auto flush_it = insertion_buffer_.lockedIterator();
      flush(&flush_it);
    }
  }

  size_t getNumKernels() const { return kernels_.size(); }

  double computeTotalCount() {
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);
    return total_count_;
  }

  double getMean() {
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);

    double acc = 0.0;
    for (const auto& kernel : kernels_) {
      acc += kernel.mean() * kernel.count();
    }

    return acc / total_count_;
  }

  // quantile is in [0, 1]
  double getQuantileEstimate(double quantile) {
    static constexpr double kError = 1e-9;
    double low = kernels_[0].mean();
    double high = kernels_[kernels_.size() - 1].mean();

    while (getQuantileOfValue(low) > quantile) {
      low -= high - low;
    }

    while (getQuantileOfValue(high) < quantile) {
      high += high - low;
    }

    double value;
    double estimated_quantile;
    do {
      value = (low + high) / 2.0;
      estimated_quantile = getQuantileOfValue(value);
      if (estimated_quantile < quantile) {
        low = value;
      } else {
        high = value;
      }
    } while (std::abs(quantile - estimated_quantile) > kError);

    return value;
  }

  double getQuantileOfValue(double value) {
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);

    double count = 0.0;
    for (const auto& kernel : kernels_) {
      count += kernel.cdf(value) * kernel.count();
    }

    return count / total_count_;
  }

  std::string debugString() {
    std::string s = "DynamicKDE{mean: " + std::to_string(getMean()) +
                    ", count: " + std::to_string(computeTotalCount()) + ",\n";
    s += "splitThreshold: " + std::to_string(splitThreshold()) + "\n";
    s += "generation: " + std::to_string(generation_) + "\n";
    for (auto kernel : kernels_) {
      s += "  " + kernel.debugString() + ",\n";
    }
    s += "}";
    return s;
  }

  std::string json(std::string title = "", std::string label = "") {
    return "";
  }

  DensityMap to_proto(std::string title = "", std::string label = "") {
    DensityMap dm;
    return dm;
  }

 private:
  const double decay_rate_;
  const size_t refresh_interval_;

  uint64_t generation_;
  uint64_t refresh_generation_;
  double total_count_;

  InsertionBuffer<double> insertion_buffer_;
  std::vector<Kernel> kernels_;
  std::vector<double> decay_factors_;

  double splitThreshold() const { return 2 * total_count_ / getNumKernels(); }

  void flush(FlushIterator<double>* flush_it) {
    for (; *flush_it != insertion_buffer_.end(); ++(*flush_it)) {
      flushValue(**flush_it);
    }
    refresh();
  }

  void flushValue(double val) {
    size_t kx = std::distance(kernels_.begin(),
                              std::upper_bound(kernels_.begin(), kernels_.end(),
                                               val, Kernel::lessThan));
    generation_++;
    int num_splits = 0;
    if (kx == 0) {
      kernels_[kx].decay(decay_factor(kernels_[kx].generation()), generation_);
      kernels_[kx].addValue(val, 1.0);
      if (kernels_[kx].count() > splitThreshold()) {
        split(kx);
        num_splits++;
      }
    } else if (kx == getNumKernels()) {
      kernels_[kx - 1].decay(decay_factor(kernels_[kx - 1].generation()),
                             generation_);
      kernels_[kx - 1].addValue(val, 1.0);
      if (kernels_[kx - 1].count() > splitThreshold()) {
        split(kx - 1);
        num_splits++;
      }
    } else {
      kernels_[kx].decay(decay_factor(kernels_[kx].generation()), generation_);
      kernels_[kx].addValue(val, 0.5);
      if (kernels_[kx].count() > splitThreshold()) {
        split(kx);
        num_splits++;
      }

      kernels_[kx - 1].decay(decay_factor(kernels_[kx - 1].generation()),
                             generation_);
      kernels_[kx - 1].addValue(val, 0.5);
      if (kernels_[kx - 1].count() > splitThreshold()) {
        split(kx - 1);
        num_splits++;
      }
    }

    while (num_splits--) {
      merge();
    }

    if (decay_rate_ != 0.0) {
      total_count_ = total_count_ * (1.0 - decay_rate_) + 1.0;
    } else {
      total_count_ = generation_;
    }
  }

  double decay_factor(uint64_t kernel_generation) const {
    if (decay_rate_ == 0.0) {
      return 1.0;
    }
    uint64_t generations = generation_ - kernel_generation;
    assert(generations < decay_factors_.size());
    return decay_factors_[generations];
  }

  void refresh() {
    if (refresh_generation_ == generation_) {
      return;
    }

    total_count_ = 0.0;
    for (size_t kx = 0; kx < kernels_.size(); kx++) {
      kernels_[kx].decay(decay_factor(kernels_[kx].generation()), generation_);
      total_count_ += kernels_[kx].count();
    }
    refresh_generation_ = generation_;
  }

  void split(size_t kx) {
    kernels_.resize(kernels_.size() + 1);

    for (size_t x = kernels_.size() - 1; x > kx; --x) {
      kernels_[x] = kernels_[x - 1];
    }

    kernels_[kx].decay(0.5, generation_);
    kernels_[kx + 1].decay(0.5, generation_);
  }

  void merge() {
    refresh();
    size_t best_kx = 0;
    double smallest_sum = kernels_[0].count() + kernels_[1].count();
    for (size_t kx = 1; kx + 1 < kernels_.size(); kx++) {
      double sum = kernels_[kx].count() + kernels_[kx + 1].count();
      if (sum < smallest_sum) {
        smallest_sum = sum;
        best_kx = kx;
      }
    }

    kernels_[best_kx] = Kernel::merge(kernels_[best_kx], kernels_[best_kx + 1]);

    for (size_t kx = best_kx + 1; kx + 1 < kernels_.size(); kx++) {
      kernels_[kx] = kernels_[kx + 1];
    }

    kernels_.pop_back();
  }
};

}  // namespace dhist
