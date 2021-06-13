Overview
--------

`DynamicDensityMap` is a collection of histogram-like datastructures that
autoscale to the inserted data. The `DynamicHistogram` offers the best
performance, particularly for quantile inferences. The `DynamicKDE2D`
datastructure provides an autoscaling map for two dimensional data.

Example
-------

Use the `show_density_map` binary to display images from a running server.

```console
~/DynamicDensityMap$ bazel run :show_density_map -- --show=1 --animate
```

![DynamicHistogram Animation Example](https://github.com/Big-Theta/DynamicDensityMap/blob/master/example/animated_hist.gif)

```console
~/DynamicDensityMap$ bazel run :show_density_map -- --show=3
```

![DynamicKDE2D Example](https://github.com/Big-Theta/DynamicDensityMap/blob/master/example/DynamicKDE2D.gif)

An example of building a binary in a separate workspace/repository and
including support for the `DensityMapServer` can be found in the `example`
folder. To start the example, use a `terminal` to execute

```console
~/DynamicDensityMap$ cd example
~/DynamicDenistyMap/example$ bazel run -c opt :profile_malloc
DensityMapServer listening on 0.0.0.0:50051
```

The `DynamicDensityMap` descriptions can be viewed with the `--list` option.

```console
~/DynamicDensityMap$ bazel run :show_density_map -- --list
...
list_density_maps_result {
  descriptions {
    identifier {
      identity: 1
    }
    title: "DynamicHistogram -- malloc"
    labels: "log(cycles)"
    timestamp {
      seconds: 1623444811
      nanos: 179904000
    }
    decay_rate: 1e-05
    type: HISTOGRAM
    num_containers: 100
  }
...
```

The value `identity: 1` indicates that this graph can be displayed with the
`--show=1` option.

Change the options of a `DynamicDensityMap` using the `--set_description`
option.  This can change the decay rate, the number of buckets/kernels, or
display options.  The map that is affected is identified through the
`{"identifier": {"identity": N}}` field.

```console
~/DynamicDensityMap$ bazel run :show_density_map -- --set_description='\
  {"identifier": {"identity": 1}, "title": "DynamicHistogram -- malloc",\
   "labels": ["log(cycles)"], "decay_rate": 0.0}'
```

Future Work
-----------

1. Add the ability to track specific quantiles to `DynamicHistogram`. Support
   for this would require that each tracked quantile be updated on each value
   insertion, but it would allow for high-accuracy tracking of extreme
   quantiles. See [this old code](https://github.com/Big-Theta/DynamicDensityMap/blob/d59e24844b08286c9595c0fb9a078627aff31739/cpp/DynamicHistogram.h#L515)
   for an example of the concept.
2. Support anomaly detection. A way to do this would be to store two
   `DynamicHistogram` datastructures with different decay rates. Then
   periodically compute the "distance" between these two histograms. Something
   like a Jaccard distance will probably work.  If that distance exceeds some
   threshold, emit an anomaly warning. This is similar in concept to a MACD
   trading indicator, except instead of tracking the convergence/divergence of
   the moving average, it would instead use the entire distribution.
3. The results for the `DynamicKDE2D` on `profile_malloc` look suspicious. The
   marginal density for the `log(size)` axis should be uniform, but it doesn't
   look uniform. What's up with that?
4. The `DynamicKDE2D` is about 10x slower than the `DynamicKDE`. This is
   largely due to the insertion algorithm. In `DynamicKDE`, a binary search can
   identify the two closest kernels, but with `DynamicKDE2D`, the current code
   uses a linear search to find the `N` closest kernels. There must be
   something better.  Even if an optimization is still `O(kernels)`, if only
   half of the kernels need to be inspected on each insertion, that would be an
   awesome optimization.
5. The process of querying the server and rendering histograms takes much
   longer than seems necessary. It's unclear what part of that pipeline is the
   bottleneck.
6. Plug both `DynamicHistogram` and `DynamicKDE2D` into a benchmark framework
   with support for the Linux perf event API. One of the ideas is that the
   performance of certain branches of code will be readily identifiable on a 2D
   density map, so meaningful information can be gleaned by inspecting just
   those regions. For example, a rare event in the code might require an
   expensive upkeep operation. That operation can be isolated and benchmarked
   even with the various fast-path operations happening. This happens with
   `malloc` when an allocation needs to consult the central free list or with
   `DynamicHistogram` which an insertion requires a flush of the
   `InsertionQueue`. Also, performance in benchmarks can be noisy due to
   various factors, such as thermal issues, or locality of threads. These types
   of events cause somewhat distinctive performance profiles to pop up, which
   could be classified and matched (see the anomaly detection idea in (2)).
7. The `DynamicKDE` algorithm is visually more accurate than the
   `DynamicHistogram`, but it is significantly slower to do anything that
   requires quantile inference because of a reliance on the `erfc` function to
   compute `CDF` values. The `DynamicKDE` is important because it generalizes
   to two dimensions whereas the `DynamicHistogram` does not. However, perhaps
   some distribution other than a Gaussian exists which will generalize to two
   dimensions, has a simple online update algorithm for the mean and covariance
   matrix, and also has a simple `CDF` equation. Maybe a normalized two
   dimensional quadratic function?
8. Implement this in `rust`.

