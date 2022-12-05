load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository", "new_git_repository")

git_repository(
    name = "protobuf",
    remote = "https://github.com/protocolbuffers/protobuf",
    tag = "v3.21.10",
)
load("@protobuf//:protobuf_deps.bzl", "protobuf_deps")
protobuf_deps()

git_repository(
    name = "com_github_grpc_grpc",
    remote = "https://github.com/grpc/grpc.git",
    tag = "v1.48.0",
)
load("@com_github_grpc_grpc//bazel:grpc_deps.bzl", "grpc_deps")
grpc_deps()

git_repository(
    name = "benchmark",
    remote = "https://github.com/google/benchmark",
    branch = "v1.7.1",
)

git_repository(
    name = "gtest",
    remote = "https://github.com/google/googletest",
    tag = "release-1.12.1",
)

load("@com_github_grpc_grpc//bazel:grpc_extra_deps.bzl", "grpc_extra_deps")
grpc_extra_deps()
