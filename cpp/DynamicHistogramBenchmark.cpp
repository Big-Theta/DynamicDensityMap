#include <unistd.h>

#include <random>

#include "DynamicHistogram.h"
#include "DynamicKDE.h"
#include "DynamicKDE2D.h"
#include "benchmark/benchmark.h"

static void BM_DynamicHistogramAddDecay(benchmark::State &state) {
  std::vector<double> xvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
  }

  dhist::DynamicHistogram uut(
      dhist::DynamicHistogramOpts().set_num_buckets(100).set_decay_rate(
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

  dhist::DynamicHistogram uut(
      dhist::DynamicHistogramOpts().set_num_buckets(100).set_decay_rate(0.0));

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

  dhist::DynamicKDE uut(
      dhist::DynamicKDEOpts().set_num_kernels(100).set_decay_rate(0.0001));

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

  dhist::DynamicKDE uut(
      dhist::DynamicKDEOpts().set_num_kernels(100).set_decay_rate(0.0));

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

void BM_DynamicKDE2DAddNoDecay(benchmark::State &state) {
  std::vector<double> xvals;
  std::vector<double> yvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
    yvals.push_back(norm(gen));
  }

  dhist::DynamicKDE2D uut(
      dhist::DynamicKDE2DOpts().set_num_kernels(100).set_decay_rate(0.0));

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
