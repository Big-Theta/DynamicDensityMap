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

#include "DensityMap.pb.h"
#include "cpp/InsertionBuffer.h"

namespace dhist {

using ::dynamic_histogram::DensityMap;
using ::google::protobuf::Timestamp;

class Kernel {
 public:
  Kernel()
      : mean_(0.0), weighted_sum_(0.0), weighted_sum_squares_(0.0), S_(0.0) {}

  ~Kernel() = default;

  Kernel(const Kernel& other) {
    mean_ = other.mean_;
    weighted_sum_ = other.weighted_sum_;
    weighted_sum_squares_ = other.weighted_sum_squares_;
    S_ = other.S_;
  }

  Kernel& operator=(const Kernel& other) {
    mean_ = other.mean_;
    weighted_sum_ = other.weighted_sum_;
    weighted_sum_squares_ = other.weighted_sum_squares_;
    S_ = other.S_;
    return *this;
  }

  void addValue(double value, double weight) {
    weighted_sum_ += weight;
    weighted_sum_squares_ += weight * weight;
    double old_mean = mean_;
    mean_ += (weight / weighted_sum_) * (value - old_mean);
    S_ += weight * (value - old_mean) * (value - mean_);
  }

  void decay(double factor) {
    weighted_sum_ *= factor;
    weighted_sum_squares_ *= factor;
    S_ *= factor;
  }

  double mean() const {
    return mean_;
  }

  double variance() const {
    return S_ / weighted_sum_;
  }

  double count() const {
    return weighted_sum_;
  }

 private:
  double mean_;
  double weighted_sum_;
  double weighted_sum_squares_;
  double S_;
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

  void addValue(double val) {}

  double computeTotalCount() { return -1.0; }

  double getMean() { return -1.0; }

  // quantile is in [0, 1]
  double getQuantileEstimate(double quantile) { return -1.0; }

  double getQuantileOfValue(double value) { return -1.0; }

  std::string debugString() { return ""; }

  std::string json(std::string title = "", std::string label = "") {
    return "";
  }

  DensityMap to_proto(std::string title = "", std::string label = "") {
    DensityMap dm;
    return dm;
  }

 private:
  const double decay_rate_;
  const int refresh_interval_;

  uint64_t generation_;
  uint64_t refresh_generation_;
  double total_count_;

  InsertionBuffer<double> insertion_buffer_;
  std::vector<Kernel> kernels_;
  std::vector<double> decay_factors_;
};

}  // namespace dhist
