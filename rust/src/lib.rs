use std::fmt;

#[derive(Debug)]
pub struct Bucket {
    // Values in a Bucket in the range [min, max).

    min: f64,
    max: f64,
    count: f64,
}

//#[derive(Debug)]
pub struct DynamicHistogram {
    decay_rate: f64,
    generation: i64,
    ubounds: Vec<f64>,
    counts: Vec<f64>,
    // quantiles: Vec<f64>,
    // quantile_locations: Vec<f64>,
}

impl DynamicHistogram {
    pub fn new(decay_rate: f64, num_buckets: usize) -> DynamicHistogram {
        DynamicHistogram {
            decay_rate: decay_rate,
            generation: 0,
            ubounds: vec![0.0; num_buckets - 1],
            counts: vec![0.0; num_buckets],
        }
    }

    pub fn add_value(&self, _val: f64) {}

    pub fn add_repeated_value(&self, _val: f64, _count: u64) {}

    pub fn clear(&self) {}

    pub fn num_buckets(&self) -> usize { return self.counts.len() }

    pub fn bucket(&self, _index: usize) -> Bucket {
        Bucket {min: 0.0, max: 0.0, count: 0.0}
    }

    pub fn min(&self) -> f64 { 0.0 }

    pub fn max(&self) -> f64 { 0.0 }

    pub fn total_count(&self) -> f64 { 0.0 }

    pub fn estimate_quantile(&self, _quantile: f64) -> f64 { 0.0 }

    pub fn is_valid(&self) -> bool { false }

    pub fn json(&self) -> String { String::from("") }
}

impl fmt::Debug for DynamicHistogram {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let _ = write!(f, "DynamicHistogram {{\n");
        let _ = write!(f, "  decay_rate: {}\n", &self.decay_rate);
        let _ = write!(f, "  generation: {}\n", &self.generation);
        let _ = write!(f, "  total_count: {}\n", &self.total_count());
        let _ = write!(f, "  is_valid: {}\n", &self.is_valid());
        let _ = write!(f, "  buckets:\n");
        for i in 0..self.num_buckets() {
            let _ = write!(f, "    {}: {:?}\n", &i, &self.bucket(i));
        }

        write!(f, "}}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;
    use rand;
    use rand_distr::{Normal, Distribution};

    // From https://en.wikipedia.org/wiki/Geometric_series#Formula
    // Args:
    //   a: First term in series.
    //   r: Rate. Decay rate is 1 - rate in this formula.
    //   n: Number of terms in series.
    fn exponential_sum(a: f64, r: f64, n: i32) -> f64 {
      a * (1.0 - r.powi(n)) / (1.0 - r)
    }

    #[test]
    fn construct() {
        let decay_rate = 0.01;
        for s in [2, 3, 5, 100].iter() {
            let size = *s as usize;
            let hist = DynamicHistogram::new(decay_rate, size);
            assert_eq!(hist.decay_rate, decay_rate);
            assert_eq!(hist.generation, 0);
            assert_eq!(hist.ubounds.len(), size - 1);
            assert_eq!(hist.counts.len(), size);
        }
    }

    #[test]
    fn is_valid() {
        struct TestCase {
            hist: DynamicHistogram,
            expected: bool,
        }

        for tc in vec![
            TestCase {
                hist: DynamicHistogram::new(0.001, 2),
                expected: true,
            },
            // A decay_rate of 0.0 is okay.
            TestCase {
                hist: DynamicHistogram::new(0.0, 3),
                expected: true,
            },
            // decay_rate must be in [0.0, 1.0).
            TestCase {
                hist: DynamicHistogram::new(1.0, 3),
                expected: false,
            },
            // Must have at least two buckets.
            TestCase {
                hist: DynamicHistogram::new(0.001, 1),
                expected: false,
            },
            // counts.len() should equal ubounds.len() + 1
            TestCase {
                hist: DynamicHistogram {
                    decay_rate: 0.0,
                    generation: 2,
                    ubounds: vec![0.0, 1.0],
                    counts: vec![1.0, 1.0],
                },
                expected: false,
            },
            // Counts must be >= 0.0
            TestCase {
                hist: DynamicHistogram {
                    decay_rate: 0.0,
                    generation: 10,
                    ubounds: vec![1.0],
                    counts: vec![1.0, -1.0],
                },
                expected: false,
            },
            // ubounds must be sorted
            TestCase {
                hist: DynamicHistogram {
                    decay_rate: 0.0,
                    generation: 4,
                    ubounds: vec![1.0, 0.0, 2.0],
                    counts: vec![1.0, 1.0, 1.0, 1.0],
                },
                expected: false,
            },
        ] {
            assert_eq!(tc.hist.is_valid(), tc.expected,
                "tc.hist.is_valid() got {}, want {}, {:?}",
                tc.hist.is_valid(), tc.expected, tc.hist);
        }
    }

    #[test]
    fn add_with_decay() {
        let decay_rate = 0.00001;
        let num_values = 100000;
        let hist = DynamicHistogram::new(decay_rate, 31);

        let normal = Normal::new(0.0, 1.0).unwrap();

        for _ in 0..num_values {
            hist.add_value(normal.sample(&mut rand::thread_rng()));
        }

        assert_approx_eq!(
            hist.total_count(),
            exponential_sum(1.0, 1.0 - decay_rate, num_values),
            1e-6);

        /*
        Using R to compute the thresholds:
        $ R -q
        > qnorm(0.5, 0, 1)
        [1] 0
        > qnorm(0.05, 0, 1)
        [1] -1.644854
        > qnorm(0.95, 0, 1)
        [1] 1.644854
        */
        assert_approx_eq!(hist.estimate_quantile(0.5), 0.0, 1e-1);
        assert_approx_eq!(hist.estimate_quantile(0.05), -1.644854, 1e-1);
        assert_approx_eq!(hist.estimate_quantile(0.95), 1.644854, 1e-1);
    }

    #[test]
    fn add_without_decay() {
        let decay_rate = 0.0;
        let num_values = 100000;
        let hist = DynamicHistogram::new(decay_rate, 31);

        let normal = Normal::new(0.0, 1.0).unwrap();

        for _ in 0..num_values {
            hist.add_value(normal.sample(&mut rand::thread_rng()));
        }

        assert_eq!(hist.total_count(), num_values as f64);
        assert_approx_eq!(hist.estimate_quantile(0.5), 0.0, 1e-1);
        assert_approx_eq!(hist.estimate_quantile(0.05), -1.644854, 1e-1);
        assert_approx_eq!(hist.estimate_quantile(0.95), 1.644854, 1e-1);
    }
}
