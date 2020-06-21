## Pull in subrepositories

This will fill in `DynamicHistogram/cpp/lib/googletest` and
`DynamicHistogram/cpp/lib/benchmark`. This is required to build the code.

```bash
$ git submodule init
$ git submodule update
```

## Build

### Check compiler

Before doing this, it's good to check your compilers are reasonable.

```bash
$ which clang
$ which clang++
$ export CC=`which clang`
$ export CXX=`which clang++`
```

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

## Tests

After building the `Debug` code:

```bash
DynamicHistogram/cpp/Debug $ test/DynamicHistogram_test
```

## Benchmarks

After building the `Release` code:

```bash
DynamicHistogram/cpp/Release $ test/DynamicHistogram_benchmark
```

