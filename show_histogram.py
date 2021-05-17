#!/usr/bin/env python3

from matplotlib import animation
from matplotlib import pyplot
from typing import Dict, List

import argparse
import DensityMap_pb2
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


def proto_to_histogram(proto: DensityMap_pb2.DynamicHistogram) -> dict:
    hist = {
        "bounds": proto.bounds,
        "counts": proto.counts,
    }
    if proto.HasField("title"):
        hist["title"] = proto.title
    if proto.HasField("label"):
        hist["label"] = proto.label
    return hist


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

        #print(help(DensityMap_pb2.DensityMap.ParseFromString))
        with open(path, "rb") as proto_in:
            serialized = proto_in.read()
        dhist = DensityMap_pb2.DensityMap()
        dhist.ParseFromString(serialized)
        hist = proto_to_histogram(dhist.dynamic_histogram)
        assert(check_histogram(hist))

        pyplot.title(hist.get("title", ""))
        label = hist.get("label", "")
        prepare_render(hist, label=label, alpha=1.0)
        pyplot.legend()
        pyplot.show()
    else:
        parser.print_help()

