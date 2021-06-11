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

#include <grpcpp/grpcpp.h>

#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "DensityMapBase.h"
#include "DynamicDensity.grpc.pb.h"
#include "DynamicDensity.pb.h"

namespace dyden {

using grpc::Server;
using grpc::ServerAsyncResponseWriter;
using grpc::ServerBuilder;
using grpc::ServerCompletionQueue;
using grpc::ServerContext;
using grpc::Status;

class DensityMapRegistry {
 public:
  static DensityMapRegistry& getInstance();

  DensityMapRegistry(DensityMapRegistry const&) = delete;
  void operator=(DensityMapRegistry const&) = delete;

  void registerDensityMap(DensityMapBase* density_map);

  void ListDensityMaps(const ::dynamic_density::ListDensityMapsParams& request,
                       ::dynamic_density::ListDensityMapsResult* reply);

  Status GetDensityMap(const DensityMapIdentifier& request, DensityMap* reply);

  Status SetDensityMapOptions(const DensityMapDescription& request,
                              DensityMapDescription* reply);

 private:
  DensityMapRegistry() : next_identity_(1) {}

  std::thread server_thread_;
  std::mutex mutex_;
  std::vector<DensityMapBase*> density_maps_;
  int32_t next_identity_;
};

// Logic and data behind the server's behavior.
class DynamicDensityServiceImpl final
    : public ::dynamic_density::DynamicDensityService::Service {
 public:
  ~DynamicDensityServiceImpl();

  void Run(std::string address);

 private:
  // Class encompasing the state and logic needed to serve a request.
  class CallData {
   public:
    CallData(::dynamic_density::DynamicDensityService::AsyncService* service,
             ServerCompletionQueue* cq);
    void Proceed();

   private:
    ::dynamic_density::DynamicDensityService::AsyncService* service_;
    ServerCompletionQueue* cq_;
    ServerContext ctx_;

    ::dynamic_density::RPCQueryParams request_;
    ::dynamic_density::RPCQueryResult reply_;

    ServerAsyncResponseWriter<::dynamic_density::RPCQueryResult> responder_;

    enum CallStatus { CREATE, PROCESS, FINISH };
    CallStatus status_;  // The current serving state.
  };

  void HandleRpcs();

  std::unique_ptr<ServerCompletionQueue> cq_;
  ::dynamic_density::DynamicDensityService::AsyncService service_;
  std::unique_ptr<Server> server_;
};

class DensityMapDaemon {
 public:
  static DensityMapDaemon& startDaemon(std::string address = "0.0.0.0:50051");

  DensityMapDaemon(DensityMapDaemon const&) = delete;
  void operator=(DensityMapDaemon const&) = delete;

 private:
  DensityMapDaemon(std::string address);
};

}  // namespace dyden
