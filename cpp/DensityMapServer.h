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

#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "DynamicDensity.grpc.pb.h"
#include "DynamicDensity.pb.h"
#include "cpp/DensityMapBase.h"

namespace dyden {

using grpc::Server;
using grpc::ServerAsyncResponseWriter;
using grpc::ServerBuilder;
using grpc::ServerCompletionQueue;
using grpc::ServerContext;
using grpc::Status;

class DensityMapRegistry {
 public:
  static DensityMapRegistry& getInstance() {
    static DensityMapRegistry instance;
    return instance;
  }

  DensityMapRegistry(DensityMapRegistry const&) = delete;
  void operator=(DensityMapRegistry const&) = delete;

  void registerDensityMap(DensityMapBase* density_map) {
    std::scoped_lock l(mutex_);
    density_map->mutable_description()->set_identity(next_identity_++);
    density_maps_.emplace_back(density_map);
  }

  void ListDensityMaps(const ::dynamic_density::ListDensityMapsParams& request,
                      ::dynamic_density::ListDensityMapsResult* reply) {
    for (auto* density_map : density_maps_) {
      ::dynamic_density::DensityMap dmap = density_map->asProto();
      if (dmap.has_dynamic_histogram()) {
        Description::copyToProto(dmap.dynamic_histogram().description(),
                                 reply->add_descriptions());
      } else {
        Description::copyToProto(dmap.dynamic_kde().description(),
                                 reply->add_descriptions());
      }
    }
  }

  Status GetDensityMap(const DensityMapIdentifier& request, DensityMap* reply) {
    int32_t id = request.identity();
    for (auto* density_map : density_maps_) {
      if (density_map->description().identifier().identity() == id) {
        density_map->toProto(reply);
        return Status::OK;
      }
    }
    return Status(grpc::StatusCode::UNKNOWN,
                  "Unknown identifier: '" + std::to_string(id) + "'");
  }

  Status SetDensityMapOptions(const DensityMapDescription& request,
                              DensityMapDescription* reply) {
    int32_t id = request.identifier().identity();
    for (auto* density_map : density_maps_) {
      if (density_map->description().identifier().identity() == id) {
        density_map->mutable_description()->setFromProto(request);
        density_map->description().toProto(reply);
        return Status::OK;
      }
    }
    return Status(grpc::StatusCode::UNKNOWN,
                  "Unknown identifier: '" + std::to_string(id) + "'");
  }

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
  ~DynamicDensityServiceImpl() {
    server_->Shutdown();
    cq_->Shutdown();
  }

  void Run() {
    std::string server_address("0.0.0.0:50051");
    ServerBuilder builder;
    builder.AddListeningPort(server_address, grpc::InsecureServerCredentials());
    builder.RegisterService(&service_);
    cq_ = builder.AddCompletionQueue();
    server_ = builder.BuildAndStart();
    std::cout << "Server listening on " << server_address << std::endl;
    HandleRpcs();
  }

 private:
  // Class encompasing the state and logic needed to serve a request.
  class CallData {
   public:
    // Take in the "service" instance (in this case representing an asynchronous
    // server) and the completion queue "cq" used for asynchronous communication
    // with the gRPC runtime.
    CallData(::dynamic_density::DynamicDensityService::AsyncService* service,
             ServerCompletionQueue* cq)
        : service_(service), cq_(cq), responder_(&ctx_), status_(CREATE) {
      // Invoke the serving logic right away.
      Proceed();
    }

    void Proceed() {
      if (status_ == CREATE) {
        // Make this instance progress to the PROCESS state.
        status_ = PROCESS;

        // As part of the initial CREATE state, we *request* that the system
        // start processing SayHello requests. In this request, "this" acts are
        // the tag uniquely identifying the request (so that different CallData
        // instances can serve different requests concurrently), in this case
        // the memory address of this CallData instance.
        service_->RequestRPCQuery(&ctx_, &request_, &responder_, cq_, cq_,
                                  this);
      } else if (status_ == PROCESS) {
        // Spawn a new CallData instance to serve new clients while we process
        // the one for this CallData. The instance will deallocate itself as
        // part of its FINISH state.
        new CallData(service_, cq_);

        if (request_.has_list_density_maps_params()) {
          DensityMapRegistry::getInstance().ListDensityMaps(
              request_.list_density_maps_params(),
              reply_.mutable_list_density_maps_result());
        } else if (request_.has_get_map_with_identifier()) {
          DensityMapRegistry::getInstance().GetDensityMap(
              request_.get_map_with_identifier(),
              reply_.mutable_density_map_result());
        } else {
          assert(request_.has_set_density_map_description());
          DensityMapRegistry::getInstance().SetDensityMapOptions(
              request_.set_density_map_description(),
              reply_.mutable_density_map_description());
        }

        // And we are done! Let the gRPC runtime know we've finished, using the
        // memory address of this instance as the uniquely identifying tag for
        // the event.
        status_ = FINISH;
        responder_.Finish(reply_, Status::OK, this);
      } else {
        GPR_ASSERT(status_ == FINISH);
        // Once in the FINISH state, deallocate ourselves (CallData).
        delete this;
      }
    }

   private:
    // The means of communication with the gRPC runtime for an asynchronous
    // server.
    ::dynamic_density::DynamicDensityService::AsyncService* service_;
    // The producer-consumer queue where for asynchronous server notifications.
    ServerCompletionQueue* cq_;
    // Context for the rpc, allowing to tweak aspects of it such as the use
    // of compression, authentication, as well as to send metadata back to the
    // client.
    ServerContext ctx_;

    // What we get from the client.
    ::dynamic_density::RPCQueryParams request_;
    // What we send back to the client.
    ::dynamic_density::RPCQueryResult reply_;

    // The means to get back to the client.
    ServerAsyncResponseWriter<::dynamic_density::RPCQueryResult> responder_;

    // Let's implement a tiny state machine with the following states.
    enum CallStatus { CREATE, PROCESS, FINISH };
    CallStatus status_;  // The current serving state.
  };

  void HandleRpcs() {
    // Spawn a new CallData instance to serve new clients.
    new CallData(&service_, cq_.get());
    void* tag;  // uniquely identifies a request.
    bool ok;
    while (true) {
      // Block waiting to read the next event from the completion queue. The
      // event is uniquely identified by its tag, which in this case is the
      // memory address of a CallData instance.
      // The return value of Next should always be checked. This return value
      // tells us whether there is any kind of event or cq_ is shutting down.
      GPR_ASSERT(cq_->Next(&tag, &ok));
      GPR_ASSERT(ok);
      static_cast<CallData*>(tag)->Proceed();
    }
  }

  std::unique_ptr<ServerCompletionQueue> cq_;
  ::dynamic_density::DynamicDensityService::AsyncService service_;
  std::unique_ptr<Server> server_;
};

class DensityMapDaemon {
 public:
  static DensityMapDaemon& startDaemon() {
    static DensityMapDaemon instance;
    return instance;
  }

  DensityMapDaemon(DensityMapDaemon const&) = delete;
  void operator=(DensityMapDaemon const&) = delete;

 private:
  DensityMapDaemon() {
    static std::thread thread(
        [](int i) {
          DynamicDensityServiceImpl server;
          server.Run();
        },
        0);
    thread.detach();
  }
};

}  // namespace dyden
