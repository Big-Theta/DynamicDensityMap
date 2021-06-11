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
#include <limits>
#include <map>
#include <utility>
#include <vector>

#include "DensityMapBase.h"
#include "DensityMapDescription.h"
#include "DensityMapServer.h"
#include "DynamicDensity.pb.h"
#include "InsertionBuffer.h"

namespace dyden {

using ::dynamic_density::DensityMap;
using ::dynamic_density::DensityMapDescription;
using ::dynamic_density::DensityMapIdentifier;
using ::google::protobuf::Timestamp;

class Kernel2D {
 public:
  Kernel2D();

  ~Kernel2D() = default;

  Kernel2D(const Kernel2D& other);

  Kernel2D& operator=(const Kernel2D& other);

  std::string debugString() const;

  void addValue(double x, double y, double weight);

  void decay(double factor, uint64_t new_generation);

  double meanX() const;

  double meanY() const;

  double varianceXX() const;

  double varianceYY() const;

  double covariance() const;

  double count() const;

  uint64_t generation() const;

  void populateProto(::dynamic_density::DynamicKDE_Kernel* proto) const;

  static Kernel2D merge(const Kernel2D& a, const Kernel2D& b);

 private:
  double mean_x_;
  double mean_y_;
  double weighted_sum_;
  double weighted_sum_squares_;
  double Sxx_;
  double Syy_;
  double Sxy_;
  uint64_t generation_;
};

class DynamicKDE2DOpts {
 public:
  DynamicKDE2DOpts()
      : num_kernels_(100),
        decay_rate_(0.0),
        refresh_interval_(512),
        title_("title"),
        labels_({"x-label", "y-label"}),
        register_with_server_(false) {}

  DynamicKDE2DOpts& setNumKernels(size_t num_kernels) {
    num_kernels_ = num_kernels;
    return *this;
  }
  size_t numKernels() const { return num_kernels_; }

  DynamicKDE2DOpts& setDecayRate(double decay_rate) {
    decay_rate_ = decay_rate;
    return *this;
  }
  double decayRate() const { return decay_rate_; }

  DynamicKDE2DOpts& setRefreshInterval(size_t refresh_interval) {
    refresh_interval_ = refresh_interval;
    return *this;
  }
  size_t refreshInterval() const { return refresh_interval_; }

  DynamicKDE2DOpts& setTitle(std::string title) {
    title_ = title;
    return *this;
  }
  std::string title() const { return title_; }

  DynamicKDE2DOpts& setLabels(std::vector<std::string> labels) {
    labels_ = labels;
    return *this;
  }
  std::vector<std::string> labels() const { return labels_; }

  DynamicKDE2DOpts& setRegisterWithServer(bool register_with_server) {
    register_with_server_ = register_with_server;
    return *this;
  }
  bool registerWithServer() const { return register_with_server_; }

 private:
  size_t num_kernels_;
  double decay_rate_;
  size_t refresh_interval_;
  std::string title_;
  std::vector<std::string> labels_;
  bool register_with_server_;
};

class DynamicKDE2D : public DensityMapBase {
 public:
  DynamicKDE2D(const DynamicKDE2DOpts& opts);

  void addValue(double x, double y);

  size_t getNumKernels() const;

  double computeTotalCount();

  std::string debugString();

  DensityMap asProto() override;

  void toProto(DensityMap* proto) override;

  void registerWithServer() override;

 private:
  uint64_t generation_;
  uint64_t refresh_generation_;
  double total_count_;
  double split_threshold_;

  InsertionBuffer<std::pair<double, double>> insertion_buffer_;
  std::vector<Kernel2D> kernels_;

  double splitThreshold() const;

  double decayRate() const;

  void flush(FlushIterator<std::pair<double, double>>* flush_it);

  void flushValue(const std::pair<double, double>& val);

  void refresh();

  void decayAll();

  void split(size_t kx);

  void merge();

  size_t closest(size_t kx);
};

}  // namespace dyden
