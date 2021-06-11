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

#include "DensityMapDescription.h"

#include <sys/time.h>

namespace dyden {

void Description::setFromProto(const DensityMapDescription& proto) {
  setTitle(proto.title());
  labels_.clear();
  for (auto label : proto.labels()) {
    labels_.push_back(label);
  }
  setDecayRate(proto.decay_rate());
  setNumContainers(proto.num_containers());
}

void Description::setDecayRate(double rate) {
  // Keep this in line with the size of the InsertionBuffer.
  decay_factors_.resize(2 * refreshInterval());
  double decay = 1.0;
  for (size_t i = 0; i < refreshInterval(); i++) {
    decay_factors_[i] = decay;
    decay *= 1.0 - rate;
  }
  decay_rate_ = rate;
}

DensityMapDescription Description::asProto() const {
  DensityMapDescription desc;
  toProto(&desc);
  return desc;
}

void Description::toProto(DensityMapDescription* proto) const {
  DensityMapIdentifier* identifier = proto->mutable_identifier();
  identifier->set_identity(identifier_.identity());

  proto->set_title(title());

  for (const auto& label : labels()) {
    proto->add_labels(label);
  }

  struct timeval tv;
  gettimeofday(&tv, NULL);
  proto->mutable_timestamp()->set_seconds(tv.tv_sec);
  proto->mutable_timestamp()->set_nanos(tv.tv_usec * 1000);

  proto->set_decay_rate(decayRate());

  if (type_ == MapType::HISTOGRAM) {
    proto->set_type(DensityMapDescription::HISTOGRAM);
  } else if (type_ == MapType::KDE) {
    proto->set_type(DensityMapDescription::KDE);
  } else {
    proto->set_type(DensityMapDescription::KDE2D);
  }

  proto->set_num_containers(numContainers());
}

/*static*/
void Description::copyToProto(const DensityMapDescription& from_proto,
                              DensityMapDescription* to_proto) {
  DensityMapIdentifier* identifier = to_proto->mutable_identifier();
  identifier->set_identity(from_proto.identifier().identity());
  to_proto->set_title(from_proto.title());
  for (const auto& label : from_proto.labels()) {
    to_proto->add_labels(label);
  }
  to_proto->mutable_timestamp()->set_seconds(from_proto.timestamp().seconds());
  to_proto->mutable_timestamp()->set_nanos(from_proto.timestamp().nanos());
  to_proto->set_decay_rate(from_proto.decay_rate());
  to_proto->set_type(from_proto.type());
  to_proto->set_num_containers(from_proto.num_containers());
}

}  // namespace dyden
