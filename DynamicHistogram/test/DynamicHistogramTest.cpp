#include <cstddef>
#include <string>

#include "gtest/gtest.h"

#include "DynamicHistogram.h"


using testing::Eq;

class DynamicHistogramReference {
public:
  DynamicHistogramReference(double decay_rate, size_t max_num_buckets)
      : decay_rate_(decay_rate), max_num_buckets_(max_num_buckets) {}
  virtual ~DynamicHistogramReference() {}

  void addValue(double val) {}
  void addRepeatedValue(double val) {}
  void clear() {}
  void merge(const DynamicHistogramReference &other) {}
  void copy(const DynamicHistogramReference &other) {}
  double getMin() { return 0.0; }
  double getMax() { return 0.0; }
  size_t getNumBuckets() { return 0; }
  Bucket getBucketByIndex(size_t idx) { return Bucket(0, 0, 0, 0); }
  double computeTotalCount() { return 0.0; }
  double getPercentileEstimate(double pct) { return 0.0; }
  std::string debugString() { return ""; }

protected:
  double decay_rate_;
  size_t max_num_buckets_;
};

TEST(DynamicHistogramTest, add) {
  static constexpr int kNumValues = 100;
  static constexpr double kVal = 100.0;
  DynamicHistogramReference uut(/*decay_rate=*/1.0, 10);

  for (int i = 0; i < kNumValues; i++) {
    uut.addValue(kVal);
  }

  EXPECT_EQ(uut.getMin(), kVal);
  EXPECT_EQ(uut.getMax(), kVal);
  EXPECT_EQ(uut.computeTotalCount(), kNumValues);
  EXPECT_EQ(uut.getPercentileEstimate(0.5), kVal);
}
