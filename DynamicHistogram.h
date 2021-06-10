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
#include <cstddef>
#include <cstdint>
#include <map>
#include <vector>

#include "DensityMapBase.h"
#include "DensityMapDescription.h"
#include "DensityMapServer.h"
#include "DynamicDensity.pb.h"
#include "InsertionBuffer.h"

namespace dyden {

using ::dynamic_density::DensityMap;

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

class DynamicHistogramOpts {
 public:
  DynamicHistogramOpts()
      : num_buckets_(100),
        decay_rate_(0.0),
        refresh_interval_(512),
        title_("title"),
        label_("label"),
        register_with_server_(false) {}

  DynamicHistogramOpts& set_num_buckets(size_t num_buckets) {
    num_buckets_ = num_buckets;
    return *this;
  }
  size_t num_buckets() const { return num_buckets_; }

  DynamicHistogramOpts& set_decay_rate(double decay_rate) {
    decay_rate_ = decay_rate;
    return *this;
  }
  double decay_rate() const { return decay_rate_; }

  DynamicHistogramOpts& set_refresh_interval(size_t refresh_interval) {
    refresh_interval_ = refresh_interval;
    return *this;
  }
  size_t refresh_interval() const { return refresh_interval_; }

  DynamicHistogramOpts& set_title(std::string title) {
    title_ = title;
    return *this;
  }
  std::string title() const { return title_; }

  DynamicHistogramOpts& set_label(std::string label) {
    label_ = label;
    return *this;
  }
  std::string label() const { return label_; }

  DynamicHistogramOpts& set_register_with_server(bool register_with_server) {
    register_with_server_ = register_with_server;
    return *this;
  }
  bool register_with_server() const { return register_with_server_; }

 private:
  size_t num_buckets_;
  double decay_rate_;
  size_t refresh_interval_;
  std::string title_;
  std::string label_;
  bool register_with_server_;
};

class DynamicHistogram : public DensityMapBase {
 public:
  DynamicHistogram(const DynamicHistogramOpts& opts);

  virtual ~DynamicHistogram() {}

  void addValue(double val);

  size_t getNumBuckets() const { return counts_.size(); }

  Bucket getBucketByIndex(size_t bx);

  double computeTotalCount();

  double getMin();

  double getMax();

  double getMean();

  // quantile is in [0, 1]
  double getQuantileEstimate(double quantile);

  double getQuantileOfValue(double value);

  std::string debugString();

  std::string json();

  DensityMap asProto() override;

  void toProto(DensityMap* proto) override;

  void registerWithServer() override;

 protected:
  uint64_t generation_;
  uint64_t refresh_generation_;
  double total_count_;
  double split_threshold_;

  InsertionBuffer<double> insertion_buffer_;

  // ubounds_ records the upper bound of a bucket. There isn't an upper bound
  // for the last bucket, so the length of counts_ will generally be one greater
  // than the length of counts_.
  std::vector<double> ubounds_;
  std::vector<double> counts_;
  std::vector<uint64_t> bucket_generation_;

  double splitThreshold() const;

  double decay_rate() const;

  double getMinRaw() const;

  double getMaxRaw() const;

  double getUpperBound(int i) const;

  void flush(FlushIterator<double>* flush_it);

  void flushValue(double val);

  void decay(size_t bx);

  void refresh();

  // Adjust bounds and add.
  // Returns:
  //   The index of the bucket that the new value landed in.
  size_t insertValue(double val);

  void split(size_t bx);

  void merge();
};

}  // namespace dyden
