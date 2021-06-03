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
  dyden::DensityMapDaemon::startDaemon();

  static constexpr size_t kNumPtrs = 4096;
  void* ptrs[kNumPtrs];
  memset(ptrs, 0, sizeof(void*) * kNumPtrs);

  dyden::DynamicHistogram dynamic_histogram(dyden::DynamicHistogramOpts()
                                                .set_num_buckets(100)
                                                .set_decay_rate(0.00001)
                                                .set_title("malloc")
                                                .set_label("log(cycles)"));
  dyden::DensityMapRegistry::getInstance().registerDynamicHistogram(
      &dynamic_histogram);

  dyden::DynamicKDE dynamic_kde(dyden::DynamicKDEOpts()
                                                .set_num_kernels(100)
                                                .set_decay_rate(0.00001)
                                                .set_title("malloc")
                                                .set_label("log(cycles)"));
  dyden::DensityMapRegistry::getInstance().registerDynamicKDE(&dynamic_kde);

  dyden::DynamicKDE2D dynamic_kde_2d(
      dyden::DynamicKDE2DOpts()
          .set_num_kernels(100)
          .set_decay_rate(0.00001)
          .set_title("malloc")
          .set_labels({"log(cycles)", "log(size)"}));
  dyden::DensityMapRegistry::getInstance().registerDynamicKDE2D(
      &dynamic_kde_2d);

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
      dynamic_histogram.addValue(log_cycles);
      dynamic_kde.addValue(log_cycles);
      dynamic_kde_2d.addValue(log_cycles, log_size);
    }
  }
}

