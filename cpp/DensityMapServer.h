// MIT License
//
// Copyright (c) 2021 Logan Evans
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <grpc/support/log.h>
#include <grpcpp/grpcpp.h>

#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "DynamicDensity.pb.h"

namespace dhist {

using grpc::Server;
using grpc::ServerAsyncResponseWriter;
using grpc::ServerBuilder;
using grpc::ServerCompletionQueue;
using grpc::ServerContext;
using grpc::Status;

using dhist::DynamicDensityService;
using dhist::DensityMapIdentifier;
using dhist::DensityMapDescription;
using dhist::ListDensityMapsParams;
using dhist::ListDensityMapsResult;
using dhist::DynamicHistogram;
using dhist::DynamicKDE;
using dhist::DensityMap;

class DensityMapsRegistry {
 public:
  static DensityMapsRegistry& getInstance() {
    static DensityMapsRegistry instance;
    return instance;
  }

  DensityMapsRegistry(DensityMapsRegistry const&) = delete;
  void operator=(DensityMapsRegistry const&) = delete;

  DynamicHistogram* registerDynamicHistogram(size_t num_buckets,
                                             double decay_rate = 0.0,
                                             size_t refresh_interval = 512) {
    std::scoped_lock(mutex_);
    dhists_.emplace_back(
        DynamicHistogram(num_buckets, decay_rate, refresh_interval));
    return &dhists_.back();
  }

  DynamicKDE* registerDynamicKDE(size_t num_buckets, double decay_rate = 0.0,
                                 size_t refresh_interval = 512) {
    std::scoped_lock(mutex_);
    dkdes_.emplace_back(DynamicKDE(num_buckets, decay_rate, refresh_interval));
    return &dkdes_.back();
  }

  DynamicKDE2D* registerDynamicKDE2D(size_t num_buckets,
                                     double decay_rate = 0.0,
                                     size_t refresh_interval = 512) {
    std::scoped_lock(mutex_);
    dkdes_.emplace_back(
        DynamicKDE2D(num_buckets, decay_rate, refresh_interval));
    return &dkde2ds_.back();
  }

  void ListDensityMaps(const ListDensityMapsParams* request,
                       ListDensityMapsResult* reply) {
  }

 private:
  DensityMapsRegistry() {}

  std::mutex mutex_;
  std::vector<DynamicHistogram> dhists_;
  std::vector<DynamicKDE> dkdes_;
  std::vector<DynamicKDE2D> dkde2ds_;
};

// Logic and data behind the server's behavior.
class DynamicDensityServiceImpl final : public DynamicDensityService::Service {
  Status ListDensityMaps(ServerContext* context,
                         const ListDensityMapsParams* request,
                         ListDensityMapsResult* reply) override {
    auto* description = reply->add_descriptions();

    std::string prefix("Hello ");
    reply->set_message(prefix + request->name());
    return Status::OK;
  }

  Status GetDensityMap(ServerContext* context,
                       const DensityMapIdentifier* request,
                       DensityMap* reply) override {
    std::string prefix("Hello ");
    reply->set_message(prefix + request->name());
    return Status::OK;
  }

  Status SetDensityMapOptions(ServerContext* context,
                              const DensityMapDescription* request,
                              DensityMap* reply) override {
    std::string prefix("Hello ");
    reply->set_message(prefix + request->name());
    return Status::OK;
  }
};

void RunServer() {
  std::string server_address("0.0.0.0:50051");
  DynamicDensityServiceImpl service;

  grpc::EnableDefaultHealthCheckService(true);
  grpc::reflection::InitProtoReflectionServerBuilderPlugin();
  ServerBuilder builder;
  // Listen on the given address without any authentication mechanism.
  builder.AddListeningPort(server_address, grpc::InsecureServerCredentials());
  // Register "service" as the instance through which we'll communicate with
  // clients. In this case it corresponds to an *synchronous* service.
  builder.RegisterService(&service);
  // Finally assemble the server.
  std::unique_ptr<Server> server(builder.BuildAndStart());
  std::cout << "Server listening on " << server_address << std::endl;

  // Wait for the server to shutdown. Note that some other thread must be
  // responsible for shutting down the server for this call to ever return.
  server->Wait();
}


//int main(int argc, char** argv) {
//  ServerImpl server;
//  server.Run();
//
//  return 0;
//}

}  // namespace dhist
