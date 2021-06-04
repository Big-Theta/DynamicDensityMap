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

#pragma once

#include <google/protobuf/util/time_util.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <vector>

#include "DynamicDensity.pb.h"
#include "DensityMapBase.h"
#include "DensityMapDescription.h"
#include "DensityMapServer.h"
#include "InsertionBuffer.h"

namespace dyden {

using ::dynamic_density::DensityMap;
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
                    ", count: " + std::to_string(count()) +
                    ", ws: " + std::to_string(weighted_sum_) +
                    ", ws2: " + std::to_string(weighted_sum_squares_) + "}";
    return s;
  }

  void addValue(double value, double weight) {
    if (weighted_sum_ == 0.0) {
      mean_ = value;
      weighted_sum_ = weight;
      weighted_sum_squares_ = weight * weight;
      S_ = 0.0;
      return;
    }

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

  double mean() const { return mean_; }

  double standardDeviation() const { return sqrt(variance()); }

  double variance() const {
    if (weighted_sum_ == 0.0) {
      return 0.0;
    }
    return S_ / weighted_sum_;
  }

  double count() const { return weighted_sum_; }

  uint64_t generation() const { return generation_; }

  double cdf(double value) const {
    double z = (value - mean()) / standardDeviation();
    return 0.5 * erfc(-z * M_SQRT1_2);
  }

  void populateProto(::dynamic_density::DynamicKDE_Kernel* proto) const {
    proto->add_coord(mean());
    proto->add_variance(variance());
    proto->set_count(count());
  }

  static bool lessThan(double value, const Kernel& b) {
    return value < b.mean_;
  }

  static Kernel merge(const Kernel& a, const Kernel& b) {
    assert(a.generation_ == b.generation_);
    Kernel kernel;
    kernel.weighted_sum_ = a.weighted_sum_ + b.weighted_sum_;
    kernel.weighted_sum_squares_ =
        a.weighted_sum_squares_ + b.weighted_sum_squares_;
    kernel.generation_ = a.generation_;
    if (kernel.weighted_sum_ == 0.0) {
      kernel.mean_ = (a.mean_ + b.mean_) / 2.0;
      kernel.S_ = 0.0;
    } else if (a.weighted_sum_ < 1.0) {
      kernel.mean_ = b.mean_;
      kernel.S_ = b.S_;
    } else if (b.weighted_sum_ < 1.0) {
      kernel.mean_ = a.mean_;
      kernel.S_ = a.S_;
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

class DynamicKDEOpts {
 public:
  DynamicKDEOpts()
      : num_kernels_(100),
        decay_rate_(0.0),
        refresh_interval_(512),
        title_("title"),
        label_("label"),
        register_with_server_(false) {}

  DynamicKDEOpts& set_num_kernels(size_t num_kernels) {
    num_kernels_ = num_kernels;
    return *this;
  }
  size_t num_kernels() const { return num_kernels_; }

  DynamicKDEOpts& set_decay_rate(double decay_rate) {
    decay_rate_ = decay_rate;
    return *this;
  }
  double decay_rate() const { return decay_rate_; }

  DynamicKDEOpts& set_refresh_interval(size_t refresh_interval) {
    refresh_interval_ = refresh_interval;
    return *this;
  }
  size_t refresh_interval() const { return refresh_interval_; }

  DynamicKDEOpts& set_title(std::string title) {
    title_ = title;
    return *this;
  }
  std::string title() const { return title_; }

  DynamicKDEOpts& set_label(std::string label) {
    label_ = label;
    return *this;
  }
  std::string label() const { return label_; }

  DynamicKDEOpts& set_register_with_server(bool register_with_server) {
    register_with_server_ = register_with_server;
    return *this;
  }
  bool register_with_server() const {
    return register_with_server_;
  }

 private:
  size_t num_kernels_;
  double decay_rate_;
  size_t refresh_interval_;
  std::string title_;
  std::string label_;
  bool register_with_server_;
};

class DynamicKDE : public DensityMapBase {
 public:
  DynamicKDE(const DynamicKDEOpts& opts)
      : DensityMapBase(DescriptionOpts()
                           .set_type(MapType::KDE)
                           .set_decay_rate(opts.decay_rate())
                           .set_refresh_interval(opts.refresh_interval())
                           .set_title(opts.title())
                           .set_labels({opts.label()})),
        generation_(0),
        refresh_generation_(0),
        total_count_(0.0),
        split_threshold_(0.0),
        insertion_buffer_(/*buffer_size=*/2 * opts.refresh_interval()) {
    kernels_.reserve(opts.num_kernels() + 2);
    kernels_.resize(opts.num_kernels());
    if (opts.register_with_server()) {
      registerWithServer();
    }
  }

  void addValue(double val) {
    size_t unflushed = insertion_buffer_.addValue(val);
    if (unflushed >= description().refresh_interval()) {
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

    return getMeanNoLock();
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
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);

    std::string s = "DynamicKDE{mean: " + std::to_string(getMeanNoLock()) +
                    ", count: " + std::to_string(computeTotalCount()) + ",\n";
    s += "splitThreshold: " + std::to_string(splitThreshold()) + "\n";
    s += "generation: " + std::to_string(generation_) + "\n";
    for (auto kernel : kernels_) {
      s += "  " + kernel.debugString() + ",\n";
    }
    s += "}";
    return s;
  }

  DensityMap asProto() override {
    DensityMap dm;
    toProto(&dm);
    return dm;
  }

  void toProto(DensityMap* proto) override {
    auto* dkde = proto->mutable_dynamic_kde();

    auto* desc = dkde->mutable_description();
    description().toProto(desc);

    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);

    for (const auto& kernel : kernels_) {
      ::dynamic_density::DynamicKDE::Kernel* k_proto = dkde->add_kernels();
      kernel.populateProto(k_proto);
    }
  }

  void registerWithServer() override {
    DensityMapRegistry::getInstance().registerDensityMap(this);
  }

 private:
  uint64_t generation_;
  uint64_t refresh_generation_;
  double total_count_;
  double split_threshold_;

  InsertionBuffer<double> insertion_buffer_;
  std::vector<Kernel> kernels_;

  double splitThreshold() const { return split_threshold_; }

  double decay_rate() const { return description().decay_rate(); }

  void flush(FlushIterator<double>* flush_it) {
    for (; *flush_it; ++(*flush_it)) {
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
      kernels_[kx].decay(
          description().decay_factor(generation_ - kernels_[kx].generation()),
          generation_);
      kernels_[kx].addValue(val, 1.0);
      if (kernels_[kx].count() > splitThreshold()) {
        split(kx);
        num_splits++;
      }
    } else if (kx == getNumKernels()) {
      kernels_[kx - 1].decay(description().decay_factor(
                                 generation_ - kernels_[kx - 1].generation()),
                             generation_);
      kernels_[kx - 1].addValue(val, 1.0);
      if (kernels_[kx - 1].count() > splitThreshold()) {
        split(kx - 1);
        num_splits++;
      }
    } else {
      double mean_left = kernels_[kx - 1].mean();
      double mean_right = kernels_[kx].mean();
      double count_left;
      if (mean_left < mean_right) {
        count_left = (val - mean_left) / (mean_right - mean_left);
      } else {
        count_left = 0.5;
      }

      kernels_[kx].decay(
          description().decay_factor(generation_ - kernels_[kx].generation()),
          generation_);
      kernels_[kx].addValue(val, 1.0 - count_left);
      if (kernels_[kx].count() > splitThreshold()) {
        split(kx);
        num_splits++;
      }

      kernels_[kx - 1].decay(description().decay_factor(
                                 generation_ - kernels_[kx - 1].generation()),
                             generation_);
      kernels_[kx - 1].addValue(val, count_left);
      if (kernels_[kx - 1].count() > splitThreshold()) {
        split(kx - 1);
        num_splits++;
      }
    }

    while (num_splits--) {
      merge();
    }

    if (decay_rate() != 0.0) {
      total_count_ = total_count_ * (1.0 - decay_rate()) + 1.0;
    } else {
      total_count_ = generation_;
    }
  }

  void refresh() {
    if (refresh_generation_ == generation_) {
      return;
    }

    total_count_ = 0.0;
    double min_count = std::numeric_limits<double>::max();
    double max_count = 0.0;
    for (size_t kx = 0; kx < kernels_.size(); kx++) {
      kernels_[kx].decay(
          description().decay_factor(generation_ - kernels_[kx].generation()),
          generation_);
      double count = kernels_[kx].count();
      total_count_ += count;
      if (count < min_count) {
        min_count = count;
      }
      if (count > max_count) {
        max_count = count;
      }
    }
    refresh_generation_ = generation_;
    if (min_count * 4 < max_count) {
      split_threshold_ = max_count;
    } else {
      split_threshold_ = 2 * total_count_ / getNumKernels();
    }
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

  double getMeanNoLock() const {
    double acc = 0.0;
    for (const auto& kernel : kernels_) {
      acc += kernel.mean() * kernel.count();
    }

    return acc / total_count_;
  }
};

}  // namespace dyden
