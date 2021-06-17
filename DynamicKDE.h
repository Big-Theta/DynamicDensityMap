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

#include "DensityMapBase.h"
#include "DensityMapDescription.h"
#include "DensityMapServer.h"
#include "DynamicDensity.pb.h"
#include "LocklessInsertionBuffer.h"

namespace dyden {

using ::dynamic_density::DensityMap;
using ::google::protobuf::Timestamp;

class Kernel {
 public:
  Kernel();
  ~Kernel() = default;

  Kernel(const Kernel& other);

  Kernel& operator=(const Kernel& other);

  std::string debugString() const;

  void addValue(double value, double weight);

  void decay(double factor, uint64_t new_generation);

  double mean() const;

  double standardDeviation() const;

  double variance() const;

  double count() const;

  uint64_t generation() const;

  double cdf(double value) const;

  void populateProto(::dynamic_density::DynamicKDE_Kernel* proto) const;

  static bool lessThan(double value, const Kernel& b);

  static Kernel merge(const Kernel& a, const Kernel& b);

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

  DynamicKDEOpts& setNumKernels(size_t num_kernels) {
    num_kernels_ = num_kernels;
    return *this;
  }
  size_t numKernels() const { return num_kernels_; }

  DynamicKDEOpts& setDecayRate(double decay_rate) {
    decay_rate_ = decay_rate;
    return *this;
  }
  double decayRate() const { return decay_rate_; }

  DynamicKDEOpts& setRefreshInterval(size_t refresh_interval) {
    refresh_interval_ = refresh_interval;
    return *this;
  }
  size_t refreshInterval() const { return refresh_interval_; }

  DynamicKDEOpts& setTitle(std::string title) {
    title_ = title;
    return *this;
  }
  std::string title() const { return title_; }

  DynamicKDEOpts& setLabel(std::string label) {
    label_ = label;
    return *this;
  }
  std::string label() const { return label_; }

  DynamicKDEOpts& setRegisterWithServer(bool register_with_server) {
    register_with_server_ = register_with_server;
    return *this;
  }
  bool registerWithServer() const { return register_with_server_; }

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
  DynamicKDE(const DynamicKDEOpts& opts);

  void addValue(double val);

  size_t getNumKernels() const;

  double computeTotalCount();

  double getMean();

  // quantile is in [0, 1]
  double getQuantileEstimate(double quantile);

  double getQuantileOfValue(double value);

  std::string debugString();

  DensityMap asProto() override;

  void toProto(DensityMap* proto) override;

  void registerWithServer() override;

 private:
  uint64_t generation_;
  uint64_t refresh_generation_;
  double total_count_;
  double split_threshold_;

  LocklessInsertionBuffer<double> insertion_buffer_;
  std::vector<Kernel> kernels_;

  double splitThreshold() const;

  double decayRate() const;

  void flush(LockedFlushIterator<double>* flush_it);

  void flushValue(double val);

  void refresh();

  void decayAll();

  void split(size_t kx);

  void merge();

  double getMeanNoLock() const;
};

}  // namespace dyden
