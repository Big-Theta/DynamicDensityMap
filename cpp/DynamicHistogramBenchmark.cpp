#include <unistd.h>

#include <random>

#include "DynamicHistogram.h"
#include "DynamicHistogramReference.h"
#include "benchmark/benchmark.h"


static void BM_ReferenceAdd(benchmark::State &state) {
  std::vector<double> vals;
  std::default_random_engine gen;
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  for (int i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    vals.push_back(unif(gen));
  }

  double decay_rates[2] = {0.0, 0.0001};
  double decay_rate = decay_rates[state.range(0)];

  dhist::DynamicHistogramReference uut(
      /*max_num_buckets=*/31,
      /*decay_rate=*/decay_rate);

  int i = 0;
  int size = vals.size();

  for (auto _ : state) {
    if (i == size) {
      i = 0;
    }

    uut.addValue(vals[i]);
  }
}
BENCHMARK(BM_ReferenceAdd)->Arg(0)->Arg(1);

static void BM_AddDecay(benchmark::State &state) {
  std::vector<double> vals;
  std::default_random_engine gen;
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  for (int i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    vals.push_back(unif(gen));
  }

  double decay_rate = 0.0001;

  dhist::DynamicHistogram uut(/*max_num_buckets=*/31, /*decay_rate=*/decay_rate);

  int i = 0;
  int size = vals.size();

  for (auto _ : state) {
    if (i == size) {
      i = 0;
    }

    uut.addValue(vals[i]);
  }
}
BENCHMARK(BM_AddDecay);

void BM_AddNoDecay(benchmark::State &state) {
  std::vector<double> vals;
  std::default_random_engine gen;
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  for (int i = 0; i < sysconf(_SC_PAGE_SIZE) / sizeof(double); i++) {
    vals.push_back(unif(gen));
  }

  double decay_rate = 0.0;

  dhist::DynamicHistogram uut(/*max_num_buckets=*/31, /*decay_rate=*/decay_rate);

  int i = 0;
  int size = vals.size();

  for (auto _ : state) {
    if (i == size) {
      i = 0;
    }

    uut.addValue(vals[i]);
  }
}
BENCHMARK(BM_AddNoDecay);

BENCHMARK_MAIN();
