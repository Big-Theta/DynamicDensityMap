Overview
--------

`DynamicDensityMap` is a collection of histogram-like datastructures that
autoscale to the inserted data. Both `DynamicHistogram` and `DynamicKDE` work
on one dimensional data while `DynamicKDE2D` works on two dimensional data.
These data structures are highly optimized and can be used to profile streaming
data.

All of these data structures use a collection of resizable and movable
component objects -- buckets for the `DynamicHistogram` and Gaussian
distributions, or kernels, for the two `KDE` data structures. A roughly equal
proportion of all data will be stored in each of these component objects, which
produces high resolution for areas of the distribution where most of the data
falls, but also provides enough resolution to perform inference about extreme
quantiles.

A decay rate tunable dictates whether old data is slowly forgotten or if all
data will be retained. When the decay rate is non-zero, this has the effect
of producing an exponential moving distribution, similar in concept to an
exponential moving average. A non-zero decay rate produces more meaningful
animations since shocks to the distribution become visible.

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

![DynamicKDE2D Example](https://github.com/Big-Theta/DynamicDensityMap/blob/master/example/DynamicKDE2D.png)

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
option. This can change the decay rate, the number of buckets/kernels, or
display options. The map that is affected is identified through the
`{"identifier": {"identity": N}}` field.

```console
~/DynamicDensityMap$ bazel run :show_density_map -- --set_description='\
  {"identifier": {"identity": 1}, "title": "DynamicHistogram -- malloc",\
   "labels": ["log(cycles)"], "decay_rate": 0.0}'
```
