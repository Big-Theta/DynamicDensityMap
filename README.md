Example
-------

```console
~/DynamicDensityMap$ bazel run :show_density_map -- --show=1 --animate
```

```console
~/DynamicDensityMap$ bazel run :show_density_map -- --show=3
```

```console
~/DynamicDensityMap$ cd example
~/DynamicDenistyMap/example$ bazel run -c opt :profile_malloc
DensityMapServer listening on 0.0.0.0:50051
```

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

```console
~/DynamicDensityMap$ bazel run :show_density_map -- --set_description='\
  {"identifier": {"identity": 1}, "title": "DynamicHistogram -- malloc",\
   "labels": ["log(cycles)"], "decay_rate": 0.0001}'
```

