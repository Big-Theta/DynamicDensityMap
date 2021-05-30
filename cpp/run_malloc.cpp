#include <cmath>
#include <random>

#include "cpp/DensityMapServer.h"
#include "cpp/DynamicKDE2D.h"

using dhist::DynamicKDE2D;

uint64_t rdtsc(){
  unsigned int lo, hi;
  __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
  return ((uint64_t)hi << 32) | lo;
}

int main() {
  static constexpr size_t kNumPtrs = 4096;
  void* ptrs[kNumPtrs];
  memset(ptrs, 0, sizeof(void*) * kNumPtrs);

  DynamicKDE2D dynamic_kde_2d(/*num_kernels=*/200, /*decay_rate=*/0.00001);

  std::default_random_engine gen;
  std::uniform_int_distribution<size_t> ptr_idx_dist(0, kNumPtrs - 1);
  // After taking e^(rand), the value will be between 1B and 1MB.
  std::uniform_real_distribution log_size_dist(0.0, 14.555);

  for (size_t i = 1; i < 1000000000; i++) {
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

      dynamic_kde_2d.addValue(log_size, end - begin);
    }

    if (i % 100000 == 0) {
      printf("%s\n", dynamic_kde_2d.asProto().DebugString().c_str());
    }
  }
}

