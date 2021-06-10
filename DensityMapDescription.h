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

#include <vector>

#include "DynamicDensity.pb.h"

namespace dyden {

using ::dynamic_density::DensityMapIdentifier;
using ::dynamic_density::DensityMapDescription;

class Identifier {
 public:
  Identifier(): identity_(0) {}

  int32_t identity() const { return identity_; }

 protected:
  friend class Description;

  void set_identity(int32_t identity) { identity_ = identity; }

 private:
  int32_t identity_;
};

enum class MapType { UNKNOWN, HISTOGRAM, KDE, KDE2D };

class Description;

class DescriptionOpts {
 public:
  DescriptionOpts()
      : type_(MapType::UNKNOWN),
        decay_rate_(0.0),
        refresh_interval_(512),
        title_("title"),
        labels_({""}) {}

  DescriptionOpts& set_type(MapType type) {
    type_ = type;
    return *this;
  }
  MapType type() const { return type_; }

  DescriptionOpts& set_decay_rate(double decay_rate) {
    decay_rate_ = decay_rate;
    return *this;
  }
  double decay_rate() const { return decay_rate_; }

  DescriptionOpts& set_refresh_interval(size_t refresh_interval) {
    refresh_interval_ = refresh_interval;
    return *this;
  }
  size_t refresh_interval() const { return refresh_interval_; }

  DescriptionOpts& set_title(std::string title) {
    title_ = title;
    return *this;
  }
  const std::string& title() const { return title_; }

  DescriptionOpts& set_labels(std::vector<std::string> labels) {
    labels_ = labels;
    return *this;
  }
  std::vector<std::string> labels() const { return labels_; }

 private:
  MapType type_;
  double decay_rate_;
  size_t refresh_interval_;
  std::string title_;
  std::vector<std::string> labels_;
};

class Description {
 public:
  Description(const DescriptionOpts& opts)
      : title_(opts.title()),
        labels_(opts.labels()),
        refresh_interval_(opts.refresh_interval()),
        type_(opts.type()) {
    set_decay_rate(opts.decay_rate());
  }

  const Identifier& identifier() const { return identifier_; }

  void setFromProto(const DensityMapDescription& proto);

  const std::string& title() const { return title_; }
  void set_title(std::string title) { title_ = title; }

  const std::vector<std::string> labels() const { return labels_; }
  void set_labels(std::vector<std::string> labels) { labels_ = labels; }

  double decay_rate() const { return decay_rate_; }
  void set_decay_rate(double rate);

  size_t refresh_interval() const { return refresh_interval_; }

  DensityMapDescription asProto() const;
  void toProto(DensityMapDescription* proto) const;

  static void copyToProto(const DensityMapDescription& from_proto,
                          DensityMapDescription* to_proto);

 protected:
  friend class DensityMapBase;
  friend class DensityMapRegistry;
  friend class DynamicHistogram;
  friend class DynamicKDE;
  friend class DynamicKDE2D;

  void set_identity(int32_t identity) {
    identifier_.set_identity(identity);
  }

  double decay_factor(uint64_t generations) const {
    assert(generations < decay_factors_.size());
    return decay_factors_[generations];
  }

 private:
  Identifier identifier_;
  std::string title_;
  std::vector<std::string> labels_;
  double decay_rate_;
  size_t refresh_interval_;
  MapType type_;

  std::vector<double> decay_factors_;
};

}  // namespace dyden
