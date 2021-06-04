#include <unistd.h>

#include <random>

#include "DynamicHistogram.h"
#include "DynamicKDE.h"
#include "DynamicKDE2D.h"
#include "benchmark/benchmark.h"

// Run on (8 X 2592 MHz CPU s)
// CPU Caches:
//   L1 Data 32 KiB (x8)
//   L1 Instruction 32 KiB (x8)
//   L2 Unified 256 KiB (x8)
//   L3 Unified 12288 KiB (x8)
// Load Average: 0.18, 1.33, 1.21
// ------------------------------------------------------------------------
// Benchmark                              Time             CPU   Iterations
// ------------------------------------------------------------------------
// BM_DynamicHistogramAddDecay         44.1 ns         44.1 ns     15381629
// BM_DynamicHistogramAddNoDecay       41.8 ns         41.8 ns     16371044
// BM_DynamicKDEAddDecay               47.9 ns         47.8 ns     14680331
// BM_DynamicKDEAddNoDecay             43.2 ns         43.2 ns     16316539
// BM_DynamicKDE2DAddDecay              365 ns          364 ns      1924448
// BM_DynamicKDE2DAddNoDecay            344 ns          344 ns      1985875

static void BM_DynamicHistogramAddDecay(benchmark::State &state) {
  std::vector<double> xvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
  }

  dyden::DynamicHistogram uut(
      dyden::DynamicHistogramOpts().set_num_buckets(100).set_decay_rate(
          0.0001));

  int x_idx = 0;
  int size = xvals.size();

  for (auto _ : state) {
    if (x_idx >= size) {
      x_idx -= size;
    }

    uut.addValue(xvals[x_idx]);
    x_idx += 1;
  }
}
BENCHMARK(BM_DynamicHistogramAddDecay);

void BM_DynamicHistogramAddNoDecay(benchmark::State &state) {
  std::vector<double> xvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
  }

  dyden::DynamicHistogram uut(
      dyden::DynamicHistogramOpts().set_num_buckets(100).set_decay_rate(0.0));

  int x_idx = 0;
  int size = xvals.size();

  for (auto _ : state) {
    if (x_idx >= size) {
      x_idx -= size;
    }

    uut.addValue(xvals[x_idx]);
    x_idx += 1;
  }
}
BENCHMARK(BM_DynamicHistogramAddNoDecay);

static void BM_DynamicKDEAddDecay(benchmark::State &state) {
  std::vector<double> xvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
  }

  dyden::DynamicKDE uut(
      dyden::DynamicKDEOpts().set_num_kernels(100).set_decay_rate(0.0001));

  int x_idx = 0;
  int size = xvals.size();

  for (auto _ : state) {
    if (x_idx >= size) {
      x_idx -= size;
    }

    uut.addValue(xvals[x_idx]);
    x_idx += 1;
  }
}
BENCHMARK(BM_DynamicKDEAddDecay);

void BM_DynamicKDEAddNoDecay(benchmark::State &state) {
  std::vector<double> xvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
  }

  dyden::DynamicKDE uut(
      dyden::DynamicKDEOpts().set_num_kernels(100).set_decay_rate(0.0));

  int x_idx = 0;
  int size = xvals.size();

  for (auto _ : state) {
    if (x_idx >= size) {
      x_idx -= size;
    }

    uut.addValue(xvals[x_idx]);
    x_idx += 1;
  }
}
BENCHMARK(BM_DynamicKDEAddNoDecay);

void BM_DynamicKDE2DAddDecay(benchmark::State &state) {
  std::vector<double> xvals;
  std::vector<double> yvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
    yvals.push_back(norm(gen));
  }

  dyden::DynamicKDE2D uut(
      dyden::DynamicKDE2DOpts().set_num_kernels(100).set_decay_rate(0.0001));

  int x_idx = 0;
  int y_idx = 0;
  int size = xvals.size();

  for (auto _ : state) {
    if (x_idx >= size) {
      x_idx -= size;
    }
    if (y_idx >= size) {
      y_idx -= size;
    }

    uut.addValue(xvals[x_idx], yvals[y_idx]);
    x_idx += 1;
    y_idx += 3;
  }
}
BENCHMARK(BM_DynamicKDE2DAddDecay);

void BM_DynamicKDE2DAddNoDecay(benchmark::State &state) {
  std::vector<double> xvals;
  std::vector<double> yvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
    yvals.push_back(norm(gen));
  }

  dyden::DynamicKDE2D uut(
      dyden::DynamicKDE2DOpts().set_num_kernels(100).set_decay_rate(0.0));

  int x_idx = 0;
  int y_idx = 0;
  int size = xvals.size();

  for (auto _ : state) {
    if (x_idx >= size) {
      x_idx -= size;
    }
    if (y_idx >= size) {
      y_idx -= size;
    }

    uut.addValue(xvals[x_idx], yvals[y_idx]);
    x_idx += 1;
    y_idx += 3;
  }
}
BENCHMARK(BM_DynamicKDE2DAddNoDecay);

BENCHMARK_MAIN();
