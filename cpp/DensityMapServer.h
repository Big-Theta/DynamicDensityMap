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
#include <grpcpp/ext/proto_server_reflection_plugin.h>
#include <grpcpp/grpcpp.h>
#include <grpcpp/health_check_service_interface.h>

#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "DynamicDensity.grpc.pb.h"
#include "DynamicDensity.pb.h"
#include "cpp/DynamicHistogram.h"
#include "cpp/DynamicKDE.h"
#include "cpp/DynamicKDE2D.h"

namespace dhist {

using grpc::Server;
using grpc::ServerAsyncResponseWriter;
using grpc::ServerBuilder;
using grpc::ServerCompletionQueue;
using grpc::ServerContext;
using grpc::Status;

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
    std::scoped_lock l(mutex_);
    dhists_.emplace_back(std::make_unique<DynamicHistogram>(
        num_buckets, decay_rate, refresh_interval, -1));
    return dhists_.back().get();
  }

  DynamicKDE* registerDynamicKDE(size_t num_buckets, double decay_rate = 0.0,
                                 size_t refresh_interval = 512) {
    std::scoped_lock l(mutex_);
    dkdes_.emplace_back(std::make_unique<DynamicKDE>(num_buckets, decay_rate,
                                                     refresh_interval, -1));
    return dkdes_.back().get();
  }

  DynamicKDE2D* registerDynamicKDE2D(size_t num_buckets,
                                     double decay_rate = 0.0,
                                     size_t refresh_interval = 512) {
    std::scoped_lock l(mutex_);
    dkde2ds_.emplace_back(std::make_unique<DynamicKDE2D>(
        num_buckets, decay_rate, refresh_interval, -1));
    return dkde2ds_.back().get();
  }

  void ListDensityMaps(const ::dynamic_density::ListDensityMapsParams* request,
                       ::dynamic_density::ListDensityMapsResult* reply) {
    for (auto& hist : dhists_) {
      Description::copyToProto(
          hist->asProto().dynamic_histogram().description(),
          reply->add_descriptions());
    }
    for (auto& kde : dkdes_) {
      Description::copyToProto(kde->asProto().dynamic_histogram().description(),
                               reply->add_descriptions());
    }
    for (auto& kde2d : dkde2ds_) {
      Description::copyToProto(
          kde2d->asProto().dynamic_histogram().description(),
          reply->add_descriptions());
    }
  }

  Status GetDensityMap(const DensityMapIdentifier* request, DensityMap* reply) {
    int32_t id = request->identity();
    for (auto& hist : dhists_) {
      if (hist->description().identifier().identifier() == id) {
        hist->toProto(reply);
        return Status::OK;
      }
    }
    for (auto& kde : dkdes_) {
      if (kde->description().identifier().identifier() == id) {
        kde->toProto(reply);
        return Status::OK;
      }
    }
    for (auto& kde2d : dkde2ds_) {
      if (kde2d->description().identifier().identifier() == id) {
        kde2d->toProto(reply);
        return Status::OK;
      }
    }
    return Status(
        grpc::StatusCode::UNKNOWN,
        "Unknown identifier: '" + std::to_string(id) + "'");
  }

  Status SetDensityMapOptions(const DensityMapDescription* request,
                              DensityMap* reply) {
    int32_t id = request->identifier().identity();
    for (auto& hist : dhists_) {
      if (hist->description().identifier().identifier() == id) {
        hist->mutable_description()->setFromProto(*request);
        hist->toProto(reply);
        return Status::OK;
      }
    }
    for (auto& kde : dkdes_) {
      if (kde->description().identifier().identifier() == id) {
        kde->mutable_description()->setFromProto(*request);
        kde->toProto(reply);
        return Status::OK;
      }
    }
    for (auto& kde2d : dkde2ds_) {
      if (kde2d->description().identifier().identifier() == id) {
        kde2d->mutable_description()->setFromProto(*request);
        kde2d->toProto(reply);
        return Status::OK;
      }
    }
    return Status(
        grpc::StatusCode::UNKNOWN,
        "Unknown identifier: '" + std::to_string(id) + "'");
  }

 private:
  DensityMapsRegistry() {}

  std::mutex mutex_;
  std::vector<std::unique_ptr<DynamicHistogram>> dhists_;
  std::vector<std::unique_ptr<DynamicKDE>> dkdes_;
  std::vector<std::unique_ptr<DynamicKDE2D>> dkde2ds_;
};

// Logic and data behind the server's behavior.
class DynamicDensityServiceImpl final
    : public ::dynamic_density::DynamicDensityService::Service {
  Status ListDensityMaps(
      ServerContext* context,
      const ::dynamic_density::ListDensityMapsParams* request,
      ::dynamic_density::ListDensityMapsResult* reply) override {
    DensityMapsRegistry::getInstance().ListDensityMaps(request, reply);
    return Status::OK;
  }

  Status GetDensityMap(ServerContext* context,
                       const DensityMapIdentifier* request,
                       DensityMap* reply) override {
    return DensityMapsRegistry::getInstance().GetDensityMap(request, reply);
  }

  Status SetDensityMapOptions(ServerContext* context,
                              const DensityMapDescription* request,
                              DensityMap* reply) override {
    return DensityMapsRegistry::getInstance().SetDensityMapOptions(request,
                                                                   reply);
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

// int main(int argc, char** argv) {
//  ServerImpl server;
//  server.Run();
//
//  return 0;
//}

}  // namespace dhist
