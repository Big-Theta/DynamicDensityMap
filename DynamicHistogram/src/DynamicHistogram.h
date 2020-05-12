#ifndef DYNAMIC_HISTOGRAM_H
#define DYNAMIC_HISTOGRAM_H

class Bucket {
public:
  // Represents data in the half-open interval [min, max).
  Bucket(double min, double max, double count)
      : min_(min), max_(max), count_(count) {}

  double diameter() const { return max() - min(); }

  double min() const { return min_; }

  double max() const { return max_; }

  double count() const { return count_; }

private:
  double min_;
  double max_;
  double count_;
};

int foo();

#endif  // DYNAMIC_HISTOGRAM_H
