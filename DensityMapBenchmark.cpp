#include <unistd.h>

#include <mutex>
#include <random>

#include "DynamicHistogram.h"
#include "DynamicKDE.h"
#include "DynamicKDE2D.h"
#include "benchmark/benchmark.h"
#include "submodules/digestible/include/digestible/digestible.h"

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
// BM_DynamicHistogramDecay         51.9 ns         51.9 ns     12778047
// BM_DynamicHistogramNoDecay       44.6 ns         44.6 ns     15385069
// BM_DynamicKDEDecay               49.4 ns         49.3 ns     14259306
// BM_DynamicKDENoDecay             50.0 ns         50.0 ns     13548123
// BM_DynamicKDE2DDecay              406 ns          406 ns      1669322
// BM_DynamicKDE2DNoDecay            395 ns          395 ns      1786269
// BM_TDigestThreadUnsafe           47.1 ns         47.1 ns     13581328
// BM_TDigestThreadSafe             58.8 ns         58.8 ns     10580836

static void BM_DynamicHistogramDecay(benchmark::State &state) {
  std::vector<double> xvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
  }

  dyden::DynamicHistogram uut(
      dyden::DynamicHistogramOpts().setNumBuckets(100).setDecayRate(0.0001));

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
BENCHMARK(BM_DynamicHistogramDecay);

void BM_DynamicHistogramNoDecay(benchmark::State &state) {
  std::vector<double> xvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
  }

  dyden::DynamicHistogram uut(
      dyden::DynamicHistogramOpts().setNumBuckets(100).setDecayRate(0.0));

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
BENCHMARK(BM_DynamicHistogramNoDecay);

static void BM_DynamicKDEDecay(benchmark::State &state) {
  std::vector<double> xvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
  }

  dyden::DynamicKDE uut(
      dyden::DynamicKDEOpts().setNumKernels(100).setDecayRate(0.0001));

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
BENCHMARK(BM_DynamicKDEDecay);

void BM_DynamicKDENoDecay(benchmark::State &state) {
  std::vector<double> xvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
  }

  dyden::DynamicKDE uut(
      dyden::DynamicKDEOpts().setNumKernels(100).setDecayRate(0.0));

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
BENCHMARK(BM_DynamicKDENoDecay);

void BM_DynamicKDE2DDecay(benchmark::State &state) {
  std::vector<double> xvals;
  std::vector<double> yvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
    yvals.push_back(norm(gen));
  }

  dyden::DynamicKDE2D uut(
      dyden::DynamicKDE2DOpts().setNumKernels(100).setDecayRate(0.0001));

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
BENCHMARK(BM_DynamicKDE2DDecay);

void BM_DynamicKDE2DNoDecay(benchmark::State &state) {
  std::vector<double> xvals;
  std::vector<double> yvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
    yvals.push_back(norm(gen));
  }

  dyden::DynamicKDE2D uut(
      dyden::DynamicKDE2DOpts().setNumKernels(100).setDecayRate(0.0));

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
BENCHMARK(BM_DynamicKDE2DNoDecay);

void BM_TDigestThreadUnsafe(benchmark::State &state) {
  std::vector<double> xvals;
  std::vector<double> yvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
    yvals.push_back(norm(gen));
  }

  digestible::tdigest<double> uut(100);

  int x_idx = 0;
  int size = xvals.size();

  for (auto _ : state) {
    if (x_idx >= size) {
      x_idx -= size;
    }

    uut.insert(xvals[x_idx]);
    x_idx += 1;
  }
}
BENCHMARK(BM_TDigestThreadUnsafe);

void BM_TDigestThreadSafe(benchmark::State &state) {
  std::vector<double> xvals;
  std::vector<double> yvals;
  std::default_random_engine gen;
  std::normal_distribution<double> norm(0.0, 1.0);

  for (size_t i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    xvals.push_back(norm(gen));
    yvals.push_back(norm(gen));
  }

  digestible::tdigest<double> uut(100);

  int x_idx = 0;
  int size = xvals.size();

  std::mutex mu;
  for (auto _ : state) {
    if (x_idx >= size) {
      x_idx -= size;
    }

    {
      std::scoped_lock l(mu);
      uut.insert(xvals[x_idx]);
    }

    x_idx += 1;
  }
}
BENCHMARK(BM_TDigestThreadSafe);

BENCHMARK_MAIN();
