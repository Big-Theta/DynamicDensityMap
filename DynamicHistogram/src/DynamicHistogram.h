#ifndef DYNAMIC_HISTOGRAM_H
#define DYNAMIC_HISTOGRAM_H

class Bucket {
public:
  // Represents data in the half-open interval [min, max).
  Bucket(double min, double max, double mean, double count)
      : min_(min), max_(max), mean_(mean), count_(count) {}

private:
  double min_;
  double max_;
  double mean_;
  double count_;
};

int foo();

#endif  // DYNAMIC_HISTOGRAM_H
