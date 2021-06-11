#!/usr/bin/env python3

from matplotlib import animation
from matplotlib import pyplot
from scipy.stats import norm, multivariate_normal
from typing import Dict, List

import argparse
import bisect
import DynamicDensity_pb2
import DynamicDensity_pb2_grpc
import grpc
import json
import numpy as np
import os
import sys

parser = argparse.ArgumentParser(description='Display a dynamic density map.')
parser.add_argument("--animate", action="store_true")
parser.add_argument("--proto", type=str, default="")
parser.add_argument("--server", type=str, default="0.0.0.0:50051")
parser.add_argument(
        "--print", action="store_true",
        help="When used with --show, also print the proto.")
parser.add_argument(
        "--list", action="store_true",
        help="List all density maps found on the server.")
parser.add_argument(
        "--show", type=int, default=0,
        help="List all density maps found on the server.")
parser.add_argument(
        "--set_description", type=str, default="",
        help='A json. The "identity" must match a density map on the server. '
             'E.g., --set_description=\''
             '{"identifier": {"identity": 1}, "title": "My Title", '
             '"labels": ["MyXLabel", "MyYLabel"], "decay_rate": 0.0001}\'')
parser.add_argument(
        "--num_points", type=int, default=100,
        help="The number of points to use for the display. "
             "A small number will reduce the processing time needed to produce "
             "graphs but will also reduce resolution.")
args = parser.parse_args()


def prepare_render_hist(
        proto: DynamicDensity_pb2.DynamicHistogram, alpha: float = 1.0):
    x_min = proto.bounds[0]
    x_max = proto.bounds[-1]
    x_d = np.linspace(x_min, x_max, args.num_points)
    f = np.zeros(len(x_d))
    total_count = sum(proto.counts)

    def cdf(x):
        c = 0.0
        idx = 0

        while proto.bounds[idx + 1] < x:
            c += proto.counts[idx]
            idx += 1

        if idx < len(proto.counts):
            c += proto.counts[idx] * (
                    (x - proto.bounds[idx]) /
                    (proto.bounds[idx + 1] - proto.bounds[idx]))
        else:
            c += proto.counts[-1]

        return c / total_count

    last_cdf = 0.0
    for i in range(1, len(x_d)):
        x = x_d[i]
        this_cdf = cdf(x)
        f[i] = this_cdf - last_cdf
        last_cdf = this_cdf

    pyplot.clf()

    title = proto.description.title
    if len(proto.description.labels):
        pyplot.xlabel(proto.description.labels[0])

    pyplot.title(title)
    pyplot.fill_between(x_d, f, alpha=alpha)


def prepare_render_dkde(proto: DynamicDensity_pb2.DynamicKDE,
                        alpha: float = 1.0):
    total = sum([kernel.count for kernel in proto.kernels])
    x_min = proto.kernels[0].coord[0] - 6 * proto.kernels[0].variance[0]
    x_max = proto.kernels[-1].coord[0] + 6 * proto.kernels[-1].variance[0]
    x_d = np.linspace(x_min, x_max, args.num_points)
    f = np.zeros(len(x_d))
    norms = [norm(k.coord[0], np.sqrt(k.variance[0])) for k in proto.kernels]

    for dist_i, dist in enumerate(norms):
        weight = proto.kernels[dist_i].count / total
        low = bisect.bisect_right(x_d, dist.mean() - 5 * dist.std())
        high = bisect.bisect_right(x_d, dist.mean() + 5 * dist.std())

        last_cdf = 0.0
        for i in range(low, high):
            cdf = dist.cdf(x_d[i])
            f[i] += (cdf - last_cdf) * weight
            last_cdf = cdf

    pyplot.clf()

    title = proto.description.title
    if len(proto.description.labels):
        pyplot.xlabel(proto.description.labels[0])

    pyplot.title(title)
    pyplot.fill_between(x_d, f, alpha=alpha)


def prepare_render_dkde_2d(proto: DynamicDensity_pb2.DynamicKDE,
                           alpha: float = 1.0):
    total = sum([kernel.count for kernel in proto.kernels])

    range_x = [proto.kernels[0].coord[0], proto.kernels[0].coord[0]]
    range_y = [proto.kernels[0].coord[1], proto.kernels[0].coord[1]]
    for kernel in proto.kernels:
        if kernel.coord[0] < range_x[0]:
            range_x[0] = kernel.coord[0]
        if kernel.coord[0] > range_x[1]:
            range_x[1] = kernel.coord[0]
        if kernel.coord[1] < range_y[0]:
            range_y[0] = kernel.coord[1]
        if kernel.coord[1] > range_y[1]:
            range_y[1] = kernel.coord[1]

    diameter_x = range_x[1] - range_x[0]
    diameter_y = range_y[1] - range_y[0]
    range_x = [range_x[0] - 0.05 * diameter_x, range_x[1] + 0.05 * diameter_x]
    range_y = [range_y[0] - 0.05 * diameter_y, range_y[1] + 0.05 * diameter_y]

    xvals = np.linspace(range_x[0], range_x[1], args.num_points)
    yvals = np.linspace(range_y[0], range_y[1], args.num_points)
    xx, yy = np.meshgrid(xvals, yvals)
    f = np.zeros((len(xx), len(yy)))

    memo = {}
    def rounded(r):
        return (round(r[0], 12), round(r[1], 12))

    def memo_or_cdf(dist, r):
        r = rounded(r)
        if (dist, r) not in memo:
            memo[(dist, r)] = dist.cdf(r)
        return memo[(dist, r)]

    def volumn_at_coord(dist, x, y, diameter_x, diameter_y):
        """The cdf provides the volumn up to a point. To pull out a square
        around a point, we need 4 calls and some additions/subtractions."""
        upper_left = memo_or_cdf(dist, (x - diameter_x, y))
        upper_right = memo_or_cdf(dist, (x, y))
        lower_left = memo_or_cdf(dist, (x - diameter_x, y - diameter_y))
        lower_right = memo_or_cdf(dist, (x, y - diameter_y))
        return upper_right - lower_right - upper_left + lower_left

    diameter_x = xvals[1] - xvals[0]
    diameter_y = yvals[1] - yvals[0]

    for kernel in proto.kernels:
        dist = multivariate_normal(
                [kernel.coord[0], kernel.coord[1]],
                [[kernel.variance[0], kernel.covariance],
                 [kernel.covariance, kernel.variance[1]]], allow_singular=True)
        weight = kernel.count / total

        xlow = bisect.bisect_right(
                xvals, kernel.coord[0] - 4 * kernel.variance[0])
        xhigh = bisect.bisect_right(
                xvals, kernel.coord[0] + 4 * kernel.variance[0])
        ylow = bisect.bisect_right(
                yvals, kernel.coord[1] - 4 * kernel.variance[1])
        yhigh = bisect.bisect_right(
                yvals, kernel.coord[1] + 4 * kernel.variance[1])

        for xi in range(xlow, xhigh):
            xval = xvals[xi]
            for yi in range(ylow, yhigh):
                yval = yvals[yi]
                f[xi][yi] += (
                    weight *
                    volumn_at_coord(dist, xval, yval, diameter_x, diameter_y))

    pyplot.clf()
    ax = pyplot.axes(projection='3d')
    pyplot.title(proto.description.title, fontsize=24)

    proto_labels = proto.description.labels
    if len(proto_labels) >= 1:
        ax.set_xlabel(proto_labels[0], fontsize=18)
    if len(proto_labels) >= 2:
        ax.set_ylabel(proto_labels[1], fontsize=18)

    return ax.plot_surface(xx, yy, f, cmap='viridis', edgecolor='none')


def query_server(request):
    with grpc.insecure_channel(args.server) as channel:
        return channel.unary_unary(
            "/dynamic_density.DynamicDensityService/RPCQuery",
            request_serializer=
                DynamicDensity_pb2.RPCQueryParams.SerializeToString,
            response_deserializer=
                DynamicDensity_pb2.RPCQueryResult.FromString,
        )(request)


def gen_from_server():
    request = DynamicDensity_pb2.RPCQueryParams()
    request.get_map_with_identifier.SetInParent()
    request.get_map_with_identifier.identity = args.show

    hist_type = None

    with grpc.insecure_channel(args.server) as channel:
        while True:
            response = channel.unary_unary(
                "/dynamic_density.DynamicDensityService/RPCQuery",
                request_serializer=
                    DynamicDensity_pb2.RPCQueryParams.SerializeToString,
                response_deserializer=
                    DynamicDensity_pb2.RPCQueryResult.FromString,
            )(request)
            ddens = response.density_map_result
            if hist_type is None:
                if ddens.HasField("dynamic_histogram"):
                    hist_type = "dynamic_histogram"
                elif len(ddens.dynamic_kde.kernels[0].coord) == 1:
                    hist_type = "dynamic_kde"
                else:
                    hist_type = "dynamic_kde_2d"

            if hist_type == "dynamic_histogram":
                prepare_render_hist(ddens.dynamic_histogram)
            elif hist_type == "dynamic_kde":
                prepare_render_dkde(ddens.dynamic_kde)
            else:
                prepare_render_dkde_2d(ddens.dynamic_kde)

            yield


def interact_with_server():
    if args.list:
        request = DynamicDensity_pb2.RPCQueryParams()
        request.list_density_maps_params.SetInParent()
        print(query_server(request))
        return
    elif args.set_description:
        request = DynamicDensity_pb2.RPCQueryParams()
        request.set_density_map_description.SetInParent()

        desc = json.loads(args.set_description)
        proto_desc = request.set_density_map_description
        proto_desc.identifier.SetInParent()
        proto_desc.identifier.identity = (
                desc.get("identifier", {}).get("identity", 0))
        proto_desc.title = desc.get("title", "")
        proto_desc.labels[:] = desc.get("labels", [])
        proto_desc.decay_rate = desc.get("decay_rate", 0.0)
        proto_desc.num_containers = desc.get("num_containers", 100)
        print(query_server(request))
        return

    assert(args.show)  # 0 is an invalid identity.

    if args.animate:
        generator = gen_from_server()
        fig = pyplot.figure()
        anim = animation.FuncAnimation(
                fig, lambda i, gen: next(gen), fargs=(generator,))
        pyplot.show()
        return

    request = DynamicDensity_pb2.RPCQueryParams()
    request.get_map_with_identifier.SetInParent()
    request.get_map_with_identifier.identity = args.show
    result = query_server(request)
    ddens = result.density_map_result
    if args.print:
        print(ddens)

    if ddens.HasField("dynamic_histogram"):
        prepare_render_hist(ddens.dynamic_histogram)
    elif len(ddens.dynamic_kde.kernels[0].coord) == 1:
        prepare_render_dkde(ddens.dynamic_kde)
    elif len(ddens.dynamic_kde.kernels[0].coord) == 2:
        prepare_render_dkde_2d(ddens.dynamic_kde)
    else:
        assert(False)

    pyplot.show()


if __name__ == "__main__":
    if args.proto:
        if args.proto.startswith('/'):
            path = args.proto
        else:
            for d in [os.environ["BUILD_WORKING_DIRECTORY"],
                      os.environ["OLDPWD"],
                      os.path.abspath(os.curdir)]:
                p = os.path.join(d, args.proto)
                if os.path.exists(p):
                    path = p
                    break
            else:
                raise FileNotFoundError(f"Unable to find file '{args.proto}'")

        with open(path, "rb") as proto_in:
            serialized = proto_in.read()

        ddens = DynamicDensity_pb2.DensityMap()
        ddens.ParseFromString(serialized)

        if ddens.HasField("dynamic_histogram"):
            prepare_render_hist(ddens.dynamic_histogram)
        else:
            if ddens.dynamic_kde.kernels[0].HasField("covariance"):
                prepare_render_dkde_2d(ddens.dynamic_kde)
            else:
                prepare_render_dkde(ddens.dynamic_kde)

        pyplot.show()
    elif args.server:
        interact_with_server()
    else:
        parser.print_help()
