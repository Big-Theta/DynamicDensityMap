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
using ::google::protobuf::Timestamp;

class Identifier {
 public:
  Identifier(): identity_(0) {}

  int32_t identity() const { return identity_; }

 protected:
  friend class Description;

  void setIdentity(int32_t identity) {
    identity_ = identity;
  }

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
  std::string title_;
  std::vector<std::string> labels_;
};

class Description {
 public:
  Description(const DescriptionOpts& opts)
      : title_(opts.title()),
        decay_rate_(opts.decay_rate()),
        type_(opts.type()) {}

  const Identifier& identifier() const { return identifier_; }

  void setFromProto(const DensityMapDescription& proto) {
    printf("> setFromProto(%s)\n", proto.DebugString().c_str());

    title_ = proto.title();

    labels_.clear();
    for (auto label : proto.labels()) {
      labels_.push_back(label);
    }

    decay_rate_ = proto.decay_rate();

    printf("< %s\n", asProto().DebugString().c_str());
  }

  const std::string& title() const { return title_; }
  void set_title(std::string title) { title_ = title; }

  const std::vector<std::string> labels() const { return labels_; }
  void set_labels(std::vector<std::string> labels) { labels_ = labels; }

  double decay_rate() const { return decay_rate_; }
  void set_decay_rate(double decay_rate) { decay_rate_ = decay_rate; }

  DensityMapDescription asProto() const {
    DensityMapDescription desc;
    toProto(&desc);
    return desc;
  }

  void toProto(DensityMapDescription* proto) const {
    DensityMapIdentifier* identifier = proto->mutable_identifier();
    identifier->set_identity(identifier_.identity());

    proto->set_title(title_);

    for (const auto& label : labels_) {
      proto->add_labels(label);
    }

    struct timeval tv;
    gettimeofday(&tv, NULL);
    proto->mutable_timestamp()->set_seconds(tv.tv_sec);
    proto->mutable_timestamp()->set_nanos(tv.tv_usec * 1000);

    proto->set_decay_rate(decay_rate());

    if (type_ == MapType::HISTOGRAM) {
      proto->set_type(DensityMapDescription::HISTOGRAM);
    } else if (type_ == MapType::KDE) {
      proto->set_type(DensityMapDescription::KDE);
    } else {
      proto->set_type(DensityMapDescription::KDE2D);
    }
  }

  static void copyToProto(const DensityMapDescription& from_proto,
                          DensityMapDescription* to_proto) {
    DensityMapIdentifier* identifier = to_proto->mutable_identifier();
    identifier->set_identity(from_proto.identifier().identity());
    to_proto->set_title(from_proto.title());
    for (const auto& label : from_proto.labels()) {
      to_proto->add_labels(label);
    }
    to_proto->mutable_timestamp()->set_seconds(
        from_proto.timestamp().seconds());
    to_proto->mutable_timestamp()->set_nanos(from_proto.timestamp().nanos());
    to_proto->set_decay_rate(from_proto.decay_rate());
    to_proto->set_type(from_proto.type());
  }

 protected:
  friend class DynamicHistogram;
  friend class DynamicKDE;
  friend class DynamicKDE2d;

  void setIdentity(int32_t identity) {
    identifier_.setIdentity(identity);
  }

 private:
  Identifier identifier_;
  std::string title_;
  std::vector<std::string> labels_;
  double decay_rate_;
  MapType type_;
};

}  // namespace dyden
