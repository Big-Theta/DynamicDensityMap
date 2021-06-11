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

#include "DynamicKDE2D.h"

namespace dyden {

Kernel2D::Kernel2D()
    : mean_x_(0.0),
      mean_y_(0.0),
      weighted_sum_(0.0),
      weighted_sum_squares_(0.0),
      Sxx_(0.0),
      Syy_(0.0),
      Sxy_(0.0),
      generation_(0) {}

Kernel2D::Kernel2D(const Kernel2D& other) {
  mean_x_ = other.mean_x_;
  mean_y_ = other.mean_y_;
  weighted_sum_ = other.weighted_sum_;
  weighted_sum_squares_ = other.weighted_sum_squares_;
  Sxx_ = other.Sxx_;
  Syy_ = other.Syy_;
  Sxy_ = other.Sxy_;
  generation_ = other.generation_;
}

Kernel2D& Kernel2D::operator=(const Kernel2D& other) {
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

std::string Kernel2D::debugString() const {
  std::string s = "Kernel2D{mean_x: " + std::to_string(meanX()) +
                  ", mean_y: " + std::to_string(meanY()) +
                  ", variance_xx: " + std::to_string(varianceXX()) +
                  ", variance_yy: " + std::to_string(varianceYY()) +
                  ", covariance: " + std::to_string(covariance()) +
                  ", weighted_sum: " + std::to_string(weighted_sum_) + "}";
  return s;
}

void Kernel2D::addValue(double x, double y, double weight) {
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

void Kernel2D::decay(double factor, uint64_t new_generation) {
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

double Kernel2D::meanX() const { return mean_x_; }

double Kernel2D::meanY() const { return mean_y_; }

double Kernel2D::varianceXX() const {
  if (weighted_sum_ == 0.0) {
    return 0.0;
  }
  return Sxx_ / weighted_sum_;
}

double Kernel2D::varianceYY() const {
  if (weighted_sum_ == 0.0) {
    return 0.0;
  }
  return Syy_ / weighted_sum_;
}

double Kernel2D::covariance() const {
  if (weighted_sum_ == 0.0) {
    return 0.0;
  }
  return Sxy_ / weighted_sum_;
}

double Kernel2D::count() const { return weighted_sum_; }

uint64_t Kernel2D::generation() const { return generation_; }

void Kernel2D::populateProto(
    ::dynamic_density::DynamicKDE_Kernel* proto) const {
  proto->add_coord(meanX());
  proto->add_coord(meanY());
  proto->add_variance(varianceXX());
  proto->add_variance(varianceYY());
  proto->set_covariance(covariance());
  proto->set_count(count());
}

/*static*/
Kernel2D Kernel2D::merge(const Kernel2D& a, const Kernel2D& b) {
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

    kernel.Sxx_ = (a.weighted_sum_ - 1.0) * a.varianceXX() +
                  (b.weighted_sum_ - 1.0) * b.varianceXX();
    kernel.Syy_ = (a.weighted_sum_ - 1.0) * a.varianceYY() +
                  (b.weighted_sum_ - 1.0) * b.varianceYY();
    kernel.Sxy_ = (a.weighted_sum_ - 1.0) * a.covariance() +
                  (b.weighted_sum_ - 1.0) * b.covariance();
  }

  kernel.weighted_sum_ = weighted_sum;
  kernel.weighted_sum_squares_ = weighted_sum_squares;
  kernel.generation_ = a.generation_;

  return kernel;
}

DynamicKDE2D::DynamicKDE2D(const DynamicKDE2DOpts& opts)
    : DensityMapBase(DescriptionOpts()
                         .setType(MapType::KDE2D)
                         .setDecayRate(opts.decayRate())
                         .setRefreshInterval(opts.refreshInterval())
                         .setTitle(opts.title())
                         .setLabels(opts.labels())
                         .setNumContainers(opts.numKernels())),
      generation_(0),
      refresh_generation_(0),
      total_count_(0.0),
      split_threshold_(0.0),
      insertion_buffer_(/*buffer_size=*/2 * opts.refreshInterval()) {
  kernels_.reserve(opts.numKernels() + 2);
  kernels_.resize(opts.numKernels());
  if (opts.registerWithServer()) {
    registerWithServer();
  }
}

void DynamicKDE2D::addValue(double x, double y) {
  size_t unflushed = insertion_buffer_.addValue(std::make_pair(x, y));
  if (unflushed >= description().refreshInterval()) {
    auto flush_it = insertion_buffer_.lockedIterator();
    flush(&flush_it);
  }
}

size_t DynamicKDE2D::getNumKernels() const { return kernels_.size(); }

double DynamicKDE2D::computeTotalCount() {
  auto flush_it = insertion_buffer_.lockedIterator();
  flush(&flush_it);
  return total_count_;
}

std::string DynamicKDE2D::debugString() {
  auto flush_it = insertion_buffer_.lockedIterator();
  flush(&flush_it);

  std::string s = "DynamicKDE2D{count: " + std::to_string(total_count_) + ",\n";
  s += "splitThreshold: " + std::to_string(splitThreshold()) + "\n";
  s += "generation: " + std::to_string(generation_) + "\n";
  for (auto kernel : kernels_) {
    s += "  " + kernel.debugString() + ",\n";
  }
  s += "}";
  return s;
}

DensityMap DynamicKDE2D::asProto() {
  DensityMap dm;
  toProto(&dm);
  return dm;
}

void DynamicKDE2D::toProto(DensityMap* proto) {
  auto* dkde2d = proto->mutable_dynamic_kde();

  auto* desc = dkde2d->mutable_description();
  description().toProto(desc);

  auto flush_it = insertion_buffer_.lockedIterator();
  flush(&flush_it);

  for (const auto& kernel : kernels_) {
    ::dynamic_density::DynamicKDE::Kernel* k_proto = dkde2d->add_kernels();
    kernel.populateProto(k_proto);
  }
}

void DynamicKDE2D::registerWithServer() {
  DensityMapRegistry::getInstance().registerDensityMap(this);
}

double DynamicKDE2D::splitThreshold() const { return split_threshold_; }

double DynamicKDE2D::decayRate() const { return description().decayRate(); }

void DynamicKDE2D::flush(FlushIterator<std::pair<double, double>>* flush_it) {
  for (; *flush_it; ++(*flush_it)) {
    flushValue(**flush_it);
  }
  refresh();
}

void DynamicKDE2D::flushValue(const std::pair<double, double>& val) {
  static constexpr size_t N = 4;
  static constexpr double kCountPerKernel = 1.0 / N;
  std::vector<size_t> best_kxs;
  best_kxs.resize(N);
  double best_distances[N];
  for (size_t i = 0; i < N; i++) {
    best_distances[i] = std::numeric_limits<double>::max();
  }

  for (size_t kx = 0; kx < kernels_.size(); kx++) {
    const auto& kernel = kernels_[kx];
    double dx = kernel.meanX() - val.first;
    double dy = kernel.meanY() - val.second;
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
    kernel->decay(description().decayFactor(generation_ - kernel->generation()),
                  generation_);
    kernel->addValue(val.first, val.second, kCountPerKernel);
    if (kernel->count() > splitThreshold()) {
      split(*kx_it);
      num_splits++;
    }
  }

  while (num_splits--) {
    merge();
  }

  if (decayRate() != 0.0) {
    total_count_ = total_count_ * (1.0 - decayRate()) + 1.0;
  } else {
    total_count_ = generation_;
  }
}

void DynamicKDE2D::refresh() {
  if (refresh_generation_ == generation_) {
    return;
  }

  refresh_generation_ = generation_;
  decayAll();

  double min_count = std::numeric_limits<double>::max();
  double max_count = 0.0;
  total_count_ = 0.0;
  for (size_t kx = 0; kx < kernels_.size(); kx++) {
    double count = kernels_[kx].count();
    total_count_ += count;
    if (count < min_count) {
      min_count = count;
    }
    if (count > max_count) {
      max_count = count;
    }
  }
  if (min_count * 4 < max_count) {
    split_threshold_ = max_count;
  } else {
    split_threshold_ = 2 * total_count_ / getNumKernels();
  }

  if (static_cast<int32_t>(getNumKernels()) != description().numContainers()) {
    while (static_cast<int32_t>(getNumKernels()) <
           description().numContainers()) {
      size_t kx = 0;
      double highest_count = kernels_[0].count();
      for (size_t i = 1; i < kernels_.size(); i++) {
        if (kernels_[i].count() > highest_count) {
          kx = i;
          highest_count = kernels_[i].count();
        }
      }
      split(kx);
    }

    while (static_cast<int32_t>(getNumKernels()) >
           description().numContainers()) {
      merge();
    }

    split_threshold_ = 2 * total_count_ / getNumKernels();
  }
}

void DynamicKDE2D::decayAll() {
  for (size_t kx = 0; kx < kernels_.size(); kx++) {
    kernels_[kx].decay(
        description().decayFactor(generation_ - kernels_[kx].generation()),
        generation_);
  }
}

void DynamicKDE2D::split(size_t kx) {
  kernels_.push_back(kernels_[kx]);
  kernels_[kx].decay(0.5, generation_);
  kernels_.back().decay(0.5, generation_);
}

void DynamicKDE2D::merge() {
  decayAll();

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

  kernels_[best_kx[0]] =
      Kernel2D::merge(kernels_[best_kx[0]], kernels_[best_kx[1]]);
  kernels_.erase(kernels_.begin() + best_kx[1]);
}

size_t DynamicKDE2D::closest(size_t kx) {
  size_t best_x = 0;
  double best_dist = std::numeric_limits<double>::max();
  for (size_t x = 0; x < kernels_.size(); x++) {
    if (x == kx) {
      continue;
    }
    double dx = kernels_[kx].meanX() - kernels_[x].meanX();
    double dy = kernels_[kx].meanY() - kernels_[x].meanY();
    double dist = dx * dx + dy * dy;
    if (dist < best_dist) {
      best_dist = dist;
      best_x = x;
    }
  }
  return best_x;
}
}  // namespace dyden
