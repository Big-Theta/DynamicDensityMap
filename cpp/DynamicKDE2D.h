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

#include "DynamicDensity.pb.h"
#include "cpp/InsertionBuffer.h"

namespace dhist {

using ::dynamic_density::DensityMap;
using ::google::protobuf::Timestamp;

class Kernel2D {
 public:
  Kernel2D()
      : mean_x_(0.0),
        mean_y_(0.0),
        weighted_sum_(0.0),
        weighted_sum_squares_(0.0),
        Sxx_(0.0),
        Syy_(0.0),
        Sxy_(0.0),
        generation_(0) {}

  ~Kernel2D() = default;

  Kernel2D(const Kernel2D& other) {
    mean_x_ = other.mean_x_;
    mean_y_ = other.mean_y_;
    weighted_sum_ = other.weighted_sum_;
    weighted_sum_squares_ = other.weighted_sum_squares_;
    Sxx_ = other.Sxx_;
    Syy_ = other.Syy_;
    Sxy_ = other.Sxy_;
    generation_ = other.generation_;
  }

  Kernel2D& operator=(const Kernel2D& other) {
    mean_x_ = other.mean_x_;
    mean_y_ = other.mean_y_;
    weighted_sum_ = other.weighted_sum_;
    weighted_sum_squares_ = other.weighted_sum_squares_;
    Sxx_ = other.Sxx_;
    Syy_ = other.Syy_;
    Sxy_ = other.Sxy_;
    generation_ = other.generation_;
    return *this;
  }

  std::string debugString() const {
    std::string s = "Kernel2D{mean_x: " + std::to_string(mean_x()) +
                    ", mean_y: " + std::to_string(mean_y()) +
                    ", variance_xx: " + std::to_string(variance_xx()) +
                    ", variance_yy: " + std::to_string(variance_yy()) +
                    ", covariance: " + std::to_string(covariance()) +
                    ", weighted_sum: " + std::to_string(weighted_sum_) + "}";
    return s;
  }

  void addValue(double x, double y, double weight) {
    weighted_sum_ += weight;
    weighted_sum_squares_ += weight * weight;

    double prop = weight / weighted_sum_;
    double dx_old = x - mean_x_;
    mean_x_ += dx_old * prop;
    double dx_new = x - mean_x_;
    Sxx_ += weight * dx_old * dx_new;

    double dy_old = y - mean_y_;
    mean_y_ += dy_old * prop;
    double dy_new = y - mean_y_;
    Syy_ += weight * dy_old * dy_new;

    Sxy_ += weight * dx_old * dy_new;
  }

  void decay(double factor, uint64_t new_generation) {
    generation_ = new_generation;
    if (factor == 1.0) {
      return;
    }
    weighted_sum_ *= factor;
    weighted_sum_squares_ *= factor;
    Sxx_ *= factor;
    Syy_ *= factor;
    Sxy_ *= factor;
  }

  double mean_x() const { return mean_x_; }

  double mean_y() const { return mean_y_; }

  double variance_xx() const {
    if (weighted_sum_ == 0.0) {
      return 0.0;
    }
    return Sxx_ / weighted_sum_;
  }

  double variance_yy() const {
    if (weighted_sum_ == 0.0) {
      return 0.0;
    }
    return Syy_ / weighted_sum_;
  }

  double covariance() const {
    if (weighted_sum_ == 0.0) {
      return 0.0;
    }
    return Sxy_ / weighted_sum_;
  }

  double count() const { return weighted_sum_; }

  uint64_t generation() const { return generation_; }

  void populateProto(::dynamic_density::DynamicKDE_Kernel* proto) const {
    proto->add_coord(mean_x());
    proto->add_coord(mean_y());
    proto->add_variance(variance_xx());
    proto->add_variance(variance_yy());
    proto->set_covariance(covariance());
    proto->set_count(count());
  }

  static Kernel2D merge(const Kernel2D& a, const Kernel2D& b) {
    assert(a.generation_ == b.generation_);
    Kernel2D kernel;
    double weighted_sum = a.weighted_sum_ + b.weighted_sum_;
    double weighted_sum_squares =
        a.weighted_sum_squares_ + b.weighted_sum_squares_;

    if (weighted_sum == 0.0) {
      kernel.mean_x_ = (a.mean_x_ + b.mean_x_) / 2.0;
      kernel.mean_y_ = (a.mean_y_ + b.mean_y_) / 2.0;
    } else if (a.weighted_sum_ < 1.0) {
      kernel = b;
    } else if (b.weighted_sum_ < 1.0) {
      kernel = a;
    } else {
      kernel.mean_x_ =
          (a.mean_x_ * a.weighted_sum_ + b.mean_x_ * b.weighted_sum_) /
          weighted_sum;
      kernel.mean_y_ =
          (a.mean_y_ * a.weighted_sum_ + b.mean_y_ * b.weighted_sum_) /
          weighted_sum;

      kernel.Sxx_ = (a.weighted_sum_ - 1.0) * a.variance_xx() +
                    (b.weighted_sum_ - 1.0) * b.variance_xx();
      kernel.Syy_ = (a.weighted_sum_ - 1.0) * a.variance_yy() +
                    (b.weighted_sum_ - 1.0) * b.variance_yy();
      kernel.Sxy_ = (a.weighted_sum_ - 1.0) * a.covariance() +
                    (b.weighted_sum_ - 1.0) * b.covariance();
    }

    kernel.weighted_sum_ = weighted_sum;
    kernel.weighted_sum_squares_ = weighted_sum_squares;
    kernel.generation_ = a.generation_;

    return kernel;
  }

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

class DynamicKDE2D {
 public:
  DynamicKDE2D(size_t num_kernels, double decay_rate = 0.0,
               size_t refresh_interval = 512)
      : decay_rate_(decay_rate),
        refresh_interval_(refresh_interval),
        generation_(0),
        refresh_generation_(0),
        total_count_(0.0),
        split_threshold_(0.0),
        insertion_buffer_(/*buffer_size=*/2 * refresh_interval) {
    kernels_.reserve(num_kernels + 2);
    kernels_.resize(num_kernels);

    if (decay_rate_ != 0.0) {
      decay_factors_.resize(refresh_interval_);
      double decay = 1.0;
      for (size_t i = 0; i < refresh_interval_; i++) {
        decay_factors_[i] = decay;
        decay *= 1.0 - decay_rate_;
      }
    }
  }

  void addValue(double x, double y) {
    size_t unflushed = insertion_buffer_.addValue(std::make_pair(x, y));
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

  std::string debugString() {
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);

    std::string s =
        "DynamicKDE2D{count: " + std::to_string(total_count_) + ",\n";
    s += "splitThreshold: " + std::to_string(splitThreshold()) + "\n";
    s += "generation: " + std::to_string(generation_) + "\n";
    for (auto kernel : kernels_) {
      s += "  " + kernel.debugString() + ",\n";
    }
    s += "}";
    return s;
  }

  std::string json(std::string title = "", std::string label = "") {
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);

    return "";
  }

  DensityMap toProto(std::string title = "", std::string label = "") {
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);

    DensityMap dm;
    auto* dkde = dm.mutable_dynamic_kde();

    if (!title.empty()) {
      dkde->set_title(title);
    }

    if (!label.empty()) {
      dkde->add_label(label);
    }

    // dkde->set_timestamp(google::protobuf::util::GetCurrentTime());

    dkde->set_decay_rate(decay_rate_);

    for (const auto& kernel : kernels_) {
      ::dynamic_density::DynamicKDE::Kernel* k_proto = dkde->add_kernels();
      kernel.populateProto(k_proto);
    }

    return dm;
  }

 private:
  const double decay_rate_;
  const size_t refresh_interval_;

  uint64_t generation_;
  uint64_t refresh_generation_;
  double total_count_;
  double split_threshold_;

  InsertionBuffer<std::pair<double, double>> insertion_buffer_;
  std::vector<Kernel2D> kernels_;
  std::vector<double> decay_factors_;

  double splitThreshold() const { return split_threshold_; }

  void flush(FlushIterator<std::pair<double, double>>* flush_it) {
    for (; *flush_it; ++(*flush_it)) {
      flushValue(**flush_it);
    }
    refresh();
  }

  void flushValue(const std::pair<double, double>& val) {
    static constexpr size_t N = 4;
    std::vector<size_t> best_kxs;
    best_kxs.resize(N);
    double best_distances[N];
    for (size_t i = 0; i < N; i++) {
      best_distances[i] = std::numeric_limits<double>::max();
    }

    for (size_t kx = 0; kx < kernels_.size(); kx++) {
      const auto& kernel = kernels_[kx];
      double dx = kernel.mean_x() - val.first;
      double dy = kernel.mean_y() - val.second;
      double dist = dx * dx + dy * dy;

      if (dist >= best_distances[N - 1]) {
        continue;
      }

      size_t idx = N - 1;
      while (idx > 0 && dist < best_distances[idx - 1]) {
        best_kxs[idx] = best_kxs[idx - 1];
        best_distances[idx] = best_distances[idx - 1];
        idx--;
      }

      best_kxs[idx] = kx;
      best_distances[idx] = dist;
    }

    generation_++;
    int num_splits = 0;

    // Process the kernels from largest index first so that splitting doesn't
    // invalidate the unprocessed indeces.
    std::sort(best_kxs.begin(), best_kxs.end());
    for (auto kx_it = best_kxs.rbegin(); kx_it != best_kxs.rend(); ++kx_it) {
      auto* kernel = &kernels_[*kx_it];
      kernel->decay(decay_factor(kernel->generation()), generation_);
      kernel->addValue(val.first, val.second, 0.25);
      if (kernel->count() > splitThreshold()) {
        split(*kx_it);
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

    double min_count = std::numeric_limits<double>::max();
    double max_count = 0.0;
    total_count_ = 0.0;
    for (size_t kx = 0; kx < kernels_.size(); kx++) {
      kernels_[kx].decay(decay_factor(kernels_[kx].generation()), generation_);
      total_count_ += kernels_[kx].count();
      if (kernels_[kx].count() < min_count) {
        min_count = kernels_[kx].count();
      }
      if (kernels_[kx].count() > max_count) {
        max_count = kernels_[kx].count();
      }
    }
    if (min_count * 4 < max_count) {
      split_threshold_ = max_count;
    } else {
      split_threshold_ = 2 * total_count_ / getNumKernels();
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
    size_t best_kx[2] = {0, closest(0)};
    double smallest_sum =
        kernels_[best_kx[0]].count() + kernels_[best_kx[1]].count();

    for (size_t kx = 1; kx < kernels_.size(); kx++) {
      size_t closest_kx = closest(kx);
      double sum = kernels_[kx].count() + kernels_[closest_kx].count();

      if (sum < smallest_sum) {
        smallest_sum = sum;
        best_kx[0] = kx;
        best_kx[1] = closest_kx;
      }
    }

    double dx = kernels_[best_kx[0]].mean_x() - kernels_[best_kx[1]].mean_x();
    double dy = kernels_[best_kx[0]].mean_y() - kernels_[best_kx[1]].mean_y();

    kernels_[best_kx[0]] =
        Kernel2D::merge(kernels_[best_kx[0]], kernels_[best_kx[1]]);
    kernels_.erase(kernels_.begin() + best_kx[1]);
  }

  size_t closest(size_t kx) {
    size_t best_x;
    double best_dist = std::numeric_limits<double>::max();
    for (size_t x = 0; x < kernels_.size(); x++) {
      if (x == kx) {
        continue;
      }
      double dx = kernels_[kx].mean_x() - kernels_[x].mean_x();
      double dy = kernels_[kx].mean_y() - kernels_[x].mean_y();
      double dist = dx * dx + dy * dy;
      if (dist < best_dist) {
        best_dist = dist;
        best_x = x;
      }
    }
    return best_x;
  }
};

}  // namespace dhist
