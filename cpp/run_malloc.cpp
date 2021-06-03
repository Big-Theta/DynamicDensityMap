#include <cmath>
#include <random>

#include "cpp/DensityMapServer.h"
#include "cpp/DynamicHistogram.h"
#include "cpp/DynamicKDE.h"
#include "cpp/DynamicKDE2D.h"

uint64_t rdtsc(){
  unsigned int lo, hi;
  __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
  return ((uint64_t)hi << 32) | lo;
}

int main() {
  dhist::StartDensityMapServer();

  static constexpr size_t kNumPtrs = 4096;
  void* ptrs[kNumPtrs];
  memset(ptrs, 0, sizeof(void*) * kNumPtrs);

  dhist::Description* description;
  auto* dynamic_histogram =
      dhist::DensityMapsRegistry::getInstance().registerDynamicHistogram(
          /*num_buckets=*/100, /*decay_rate=*/0.00001);
  description = dynamic_histogram->mutable_description();
  description->set_title("malloc");
  description->set_labels({"log(cycles)"});

  auto* dynamic_kde =
      dhist::DensityMapsRegistry::getInstance().registerDynamicKDE(
          /*num_kernels=*/100, /*decay_rate=*/0.00001);
  description = dynamic_kde->mutable_description();
  description->set_title("malloc");
  description->set_labels({"log(cycles)"});

  auto* dynamic_kde_2d =
      dhist::DensityMapsRegistry::getInstance().registerDynamicKDE2D(
          /*num_kernels=*/100, /*decay_rate=*/0.00001);
  description = dynamic_kde_2d->mutable_description();
  description->set_title("malloc");
  description->set_labels({"log(cycles)", "log(size)"});

  std::default_random_engine gen;
  std::uniform_int_distribution<size_t> ptr_idx_dist(0, kNumPtrs - 1);
  std::uniform_real_distribution log_size_dist(0.0, log(1024 * 1024));

  while (true) {
    size_t ptr_idx = ptr_idx_dist(gen);
    void* ptr = ptrs[ptr_idx];
    if (ptr != nullptr) {
      free(ptr);
      ptrs[ptr_idx] = nullptr;
    } else {
      double log_size = log_size_dist(gen);
      size_t size = static_cast<size_t>(exp(log_size));
      uint64_t begin = rdtsc();
      ptr = malloc(size);
      uint64_t end = rdtsc();
      ptrs[ptr_idx] = ptr;

      double log_cycles = log(end - begin);
      dynamic_histogram->addValue(log_cycles);
      dynamic_kde->addValue(log_cycles);
      dynamic_kde_2d->addValue(log_cycles, log_size);
    }
  }
}

