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

#include <atomic>
#include <cmath>
#include <random>
#include <thread>

#include "DensityMapServer.h"
#include "DynamicHistogram.h"
#include "DynamicKDE.h"
#include "DynamicKDE2D.h"

uint64_t rdtsc() {
  unsigned int lo, hi;
  __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
  return ((uint64_t)hi << 32) | lo;
}

int main() {
  dyden::DensityMapDaemon::startDaemon();

  static constexpr std::chrono::duration kNumSeconds =
      std::chrono::seconds(600);
  static constexpr size_t kNumPtrs = 4096;
  static constexpr int kNumThreads = 5;
  std::thread threads[kNumThreads];

  std::atomic<void*> ptrs[kNumPtrs];
  for (size_t i = 0; i < kNumPtrs; i++) {
    ptrs[i].store(nullptr, std::memory_order_release);
  }

  dyden::DynamicHistogram dynamic_histogram(
      dyden::DynamicHistogramOpts()
          .setNumBuckets(100)
          .setDecayRate(0.00001)
          .setTitle("malloc")
          .setLabel("log(cycles)")
          .setRegisterWithServer(true));
  dyden::DynamicKDE dynamic_kde(dyden::DynamicKDEOpts()
                                    .setNumKernels(100)
                                    .setDecayRate(0.00001)
                                    .setTitle("malloc")
                                    .setLabel("log(cycles)")
                                    .setRegisterWithServer(true));
  dyden::DynamicKDE2D dynamic_kde_2d(
      dyden::DynamicKDE2DOpts()
          .setNumKernels(75)
          .setDecayRate(0.00001)
          .setTitle("malloc")
          .setLabels({"log(cycles)", "log(size)"})
          .setRegisterWithServer(true));

  for (int tx = 0; tx < kNumThreads; tx++) {
    threads[tx] = std::thread(
        [&](int tx) {
          std::default_random_engine gen;
          std::uniform_int_distribution<size_t> ptr_idx_dist(0, kNumPtrs - 1);
          std::uniform_real_distribution log_size_dist(0.0, log(1024 * 1024));
          gen.seed(tx);

          const auto start_time = std::chrono::system_clock::now();
          uint64_t iteration = 0;
          while (true) {
            double log_size = log_size_dist(gen);
            size_t size = static_cast<size_t>(exp(log_size));
            uint64_t begin = rdtsc();
            void* ptr = malloc(size);
            uint64_t end = rdtsc();

            free(ptrs[ptr_idx_dist(gen)].exchange(ptr));

            double log_cycles = log(end - begin);
            dynamic_histogram.addValue(log_cycles);
            dynamic_kde.addValue(log_cycles);
            dynamic_kde_2d.addValue(log_cycles, log_size);

            if ((iteration++ & (1 << 16)) == 0) {
              if (std::chrono::system_clock::now() - start_time >=
                  kNumSeconds) {
                break;
              }
            }
          }
        },
        tx);
  }

  for (int tx = 0; tx < kNumThreads; tx++) {
    threads[tx].join();
  }
}
