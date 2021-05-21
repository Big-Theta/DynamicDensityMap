#!/usr/bin/env python3

from matplotlib import animation
from matplotlib import pyplot
from scipy.stats import norm
from typing import Dict, List

import argparse
import bisect
import DynamicDensity_pb2
import json
import numpy as np
import os
import sys

parser = argparse.ArgumentParser(description='Display a dynamic density map.')
parser.add_argument("--stdin", type=bool, default=False)
parser.add_argument("--animate", type=bool, default=False)
parser.add_argument("--proto", type=str, default="")
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
    if proto.HasField("title"):
        hist["title"] = proto.title
    if proto.HasField("label"):
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

    title = ""
    if proto.HasField("title"):
        title = proto.title
    label = ""
    if len(proto.label):
        label = proto.label[0]

    pyplot.title(title)
    pyplot.fill_between(x_d, y_d, alpha=alpha, label=label)


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

        #print(help(DynamicDensity_pb2.DynamicDensity.ParseFromString))
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
            prepare_render_dkde(ddens.dynamic_kde)

        pyplot.legend()
        pyplot.show()
    else:
        parser.print_help()
