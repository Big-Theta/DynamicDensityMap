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
parser.add_argument("--stdin", type=bool, default=False)
parser.add_argument("--animate", type=bool, default=False)
parser.add_argument("--proto", type=str, default="")
parser.add_argument("--server", type=str, default="0.0.0.0:50051")
parser.add_argument("--list", type=bool, default=False)
parser.add_argument("--get_map", type=int, default=0)
args = parser.parse_args()


def check_histogram(data: Dict) -> bool:
    if len(data["bounds"]) == len(data["counts"]) + 1:
        return True
    return False


def extract_histograms(data: str) -> List[dict]:
    candidates = []
    open_brace = None
    nest_count = 0
    data = data.replace("\n", "")
    for i, ch in enumerate(data):
        if ch == "{":
            if nest_count == 0:
                open_brace = i
            nest_count += 1
        elif ch == "}":
            if nest_count == 0:
                continue
            elif nest_count == 1:
                candidates.append(data[open_brace:i+1])
            nest_count -= 1

    # Pick out only json.
    cleaned = []
    for candidate in candidates:
        try:
            cleaned.append(json.loads(candidate))
        except json.decoder.JSONDecodeError:
            pass

    # Make sure that only json objects with the keys "bounds" and
    # "counts" are present.
    return [c for c in cleaned if check_histogram(c)]


def dhist_to_histogram(proto: DynamicDensity_pb2.DynamicHistogram) -> dict:
    hist = {
        "bounds": proto.bounds,
        "counts": proto.counts,
    }
    hist["title"] = proto.title
    if len(proto.label):
        hist["label"] = proto.label
    return hist


def prepare_render_dkde(proto: DynamicDensity_pb2.DynamicKDE,
                        alpha: float = 1.0):
    total = sum([kernel.count for kernel in proto.kernels])
    x_min = proto.kernels[0].coord[0] - 6 * proto.kernels[0].variance[0]
    x_max = proto.kernels[-1].coord[0] + 6 * proto.kernels[-1].variance[0]
    x_d = np.linspace(x_min, x_max, 1000)
    y_d = np.zeros(1000)
    norms = [norm(k.coord[0], np.sqrt(k.variance[0])) for k in proto.kernels]
    for dist_i, dist in enumerate(norms):
        weight = proto.kernels[dist_i].count / total
        low = bisect.bisect_right(x_d, dist.mean() - 6 * dist.std())
        high = bisect.bisect_right(x_d, dist.mean() + 6 * dist.std())
        for i in range(low, high):
            y_d[i] += dist.pdf(x_d[i]) * weight

    title = proto.title
    label = ""
    if len(proto.label):
        label = proto.label[0]

    pyplot.title(title)
    pyplot.fill_between(x_d, y_d, alpha=alpha, label=label)
    pyplot.legend()


def prepare_render_dkde_2d(proto: DynamicDensity_pb2.DynamicKDE,
                           alpha: float = 1.0):
    total = sum([kernel.count for kernel in proto.kernels])

    def range_x_range_y(kernel):
        range_x = [kernel.coord[0] - 6 * kernel.variance[0],
                   kernel.coord[0] + 6 * kernel.variance[0]]
        range_y = [kernel.coord[1] - 6 * kernel.variance[1],
                   kernel.coord[1] + 6 * kernel.variance[1]]
        return range_x, range_y

    range_x, range_y = range_x_range_y(proto.kernels[0])
    for kernel in proto.kernels:
        pos_range_x, pos_range_y = range_x_range_y(kernel)
        range_x[0] = min(range_x[0], pos_range_x[0])
        range_x[1] = max(range_x[1], pos_range_x[1])
        range_y[0] = min(range_y[0], pos_range_y[0])
        range_y[1] = max(range_y[1], pos_range_y[1])

    # XXX
    range_x = [0, 4]
    range_y = [0, 1000]

    xvals = np.linspace(range_x[0], range_x[1], 100)
    yvals = np.linspace(range_y[0], range_y[1], 100)
    xx, yy = np.meshgrid(xvals, yvals)
    f = np.zeros((100, 100))

    for kernel in proto.kernels:
        dist = multivariate_normal(
                [kernel.coord[0], kernel.coord[1]],
                [[kernel.variance[0], kernel.covariance],
                 [kernel.covariance, kernel.variance[1]]], allow_singular=True)
        weight = kernel.count / total

        xlow = bisect.bisect_right(
                xvals, kernel.coord[0] - 6 * kernel.variance[0])
        xhigh = bisect.bisect_right(
                xvals, kernel.coord[0] + 6 * kernel.variance[0])
        ylow = bisect.bisect_right(
                yvals, kernel.coord[1] - 6 * kernel.variance[1])
        yhigh = bisect.bisect_right(
                yvals, kernel.coord[1] + 6 * kernel.variance[1])

        for xi in range(xlow, xhigh):
            xval = xvals[xi]
            for yi in range(ylow, yhigh):
                yval = yvals[yi]
                f[xi][yi] += weight * dist.pdf([xval, yval])

    pyplot.title(proto.description.title)

    ax = pyplot.axes(projection='3d')
    ax.plot_surface(xx, yy, f,
                    cmap='viridis', edgecolor='none')
    ax.set_zlim(0, 0.0002)


def gen_histograms(data: str):
    i = 0
    open_brace = None
    nest_count = 0
    while True:
        candidate = None
        #data = data.replace("\n", "")
        while i < len(data):
            ch = data[i]
            i += 1
            if ch == "{":
                if nest_count == 0:
                    open_brace = i - 1
                nest_count += 1
            elif ch == "}":
                if nest_count == 0:
                    continue
                elif nest_count == 1:
                    candidate = data[open_brace:i]
                    nest_count -= 1
                    break
                nest_count -= 1

        if not candidate:
            return

        # Pick out only json.
        try:
            cleaned = json.loads(candidate)
        except json.decoder.JSONDecodeError:
            continue

        if check_histogram(cleaned):
            yield cleaned


def compute_mean(histogram: Dict) -> float:
    acc = 0.0
    total = 0.0
    for i in range(len(histogram["bounds"]) - 1):
        acc += (histogram["counts"][i]
                * (histogram["bounds"][i] + histogram["bounds"][i + 1])
                / 2)
        total += histogram["counts"][i]
    return acc / total


colors = pyplot.rcParams["axes.prop_cycle"].by_key()["color"]
color_index = 0

def prepare_render(histogram: Dict, label: str = "", alpha: float = 1.0):
    global color_index

    color = colors[color_index]
    color_index += 1
    if color_index >= len(colors):
        color_index = 0

    counts = np.array(histogram["counts"])
    total_count = sum(counts)
    weights = np.array([count / total_count for count in counts])
    bins = np.array(histogram["bounds"])
    widths = bins[1:] - bins[:-1]
    heights = weights.astype(np.float) / widths
    pyplot.fill_between(
            bins.repeat(2)[1:-1], heights.repeat(2),
            color=color, alpha=alpha, label=label)
    pyplot.legend()

def next_frame(i, hist_generator):
    global color_index
    color = colors[color_index]

    pyplot.clf()

    for _ in range(25):
        histogram = next(hist_generator)

    counts = np.array(histogram["counts"])
    total_count = sum(counts)
    weights = np.array([count / total_count for count in counts])
    bins = np.array(histogram["bounds"])
    widths = bins[1:] - bins[:-1]
    heights = weights.astype(np.float) / widths

    ax = pyplot.axes()
    return ax.fill_between(
            bins.repeat(2)[1:-1], heights.repeat(2),
            color=color)

def animate(data: str):
    fig = pyplot.figure()
    anim = animation.FuncAnimation(
            fig, next_frame, fargs=(gen_histograms(data),))
    pyplot.show()

if __name__ == "__main__":
    if args.stdin:
        data = sys.stdin.read()

        if args.animate:
            animate(data)
            exit()

        hists = extract_histograms(data)
        title = None

        for hist in hists:
            if check_histogram(hist):
                label = hist.get("label", "")

                if hist.get("title"):
                    title = hist.get("title")

                prepare_render(hist, label=label, alpha=1.0 / len(hists))

        pyplot.legend()
        if title:
            pyplot.title(title)
        pyplot.show()

    elif args.proto:
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
            hist = dhist_to_histogram(ddens.dynamic_histogram)
            assert(check_histogram(hist))

            pyplot.title(hist.get("title", ""))
            label = hist.get("label", "")
            prepare_render(hist, label=label, alpha=1.0)
        else:
            if ddens.dynamic_kde.kernels[0].HasField("covariance"):
                prepare_render_dkde_2d(ddens.dynamic_kde)
            else:
                prepare_render_dkde(ddens.dynamic_kde)

        pyplot.show()
    elif args.server:
        if args.list:
            request = DynamicDensity_pb2.RPCQueryParams()
            request.list_density_maps_params.SetInParent()
        elif args.get_map:
            request = DynamicDensity_pb2.RPCQueryParams()
            request.get_map_with_identifier.SetInParent()
            request.get_map_with_identifier.identity = args.get_map
        elif args.set_map:
            assert(False)

        with grpc.insecure_channel(args.server) as channel:
            response = channel.unary_unary(
                "/dynamic_density.DynamicDensityService/RPCQuery",
                request_serializer=
                    DynamicDensity_pb2.RPCQueryParams.SerializeToString,
                response_deserializer=
                    DynamicDensity_pb2.RPCQueryResult.FromString,
            )(request)

        if response.HasField("list_density_maps_result"):
            print(response)
        elif response.HasField("density_map_result"):
            ddens = response.density_map_result
            if ddens.HasField("dynamic_histogram"):
                hist = dhist_to_histogram(ddens.dynamic_histogram)
                assert(check_histogram(hist))

                pyplot.title(hist.get("title", ""))
                label = hist.get("label", "")
                prepare_render(hist, label=label, alpha=1.0)
            else:
                if len(ddens.dynamic_kde.kernels[0].coord) == 1:
                    prepare_render_dkde(ddens.dynamic_kde)
                else:
                    prepare_render_dkde_2d(ddens.dynamic_kde)

            pyplot.show()
    else:
        parser.print_help()
