## Pull in subrepositories

This will fill in `DynamicHistogram/cpp/lib/googletest` and
`DynamicHistogram/cpp/lib/benchmark`. This is required to build the code.

$ git submodule init
$ git submodule update

## Build

### Release

```bash
DynamicHistogram/cpp $ mkdir Release
DyanmicHistogram/cpp $ cd Release
DynamicHistogram/cpp/Release $ cmake -DCMAKE_BUILD_TYPE=Release ..
DynamicHistogram/cpp/Release $ make
```

### Debug

```bash
DynamicHistogram/cpp $ mkdir Debug
DyanmicHistogram/cpp $ cd Debug
DynamicHistogram/cpp/Debug $ cmake -DCMAKE_BUILD_TYPE=Debug ..
DynamicHistogram/cpp/Debug $ make
```

