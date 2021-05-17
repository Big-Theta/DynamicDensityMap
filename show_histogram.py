#!/usr/bin/env python3

from cpp import DensityMap_pb2
from matplotlib import animation
from matplotlib import pyplot
from typing import Dict, List

import argparse
import json
import numpy as np
import os
import sys

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("--stdin", type=bool, default=False)
parser.add_argument("--animage", type=bool, default=False)
parser.add_argument("--proto", type=str, default="")
args = parser.parse_args()


def maybe_histogram(data: Dict) -> bool:
    if (isinstance(data.get("bounds"), list)
        and isinstance(data.get("counts"), list)
        and len(data["bounds"]) == len(data["counts"]) + 1
    ):
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
    return [c for c in cleaned if maybe_histogram(c)]


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

        if maybe_histogram(cleaned):
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
        print(hists)
        title = None

        for hist in hists:
            if maybe_histogram(hist):
                label = hist.get("label", "")

                if hist.get("title"):
                    title = hist.get("title")

                prepare_render(hist, label=label, alpha=1.0 / len(hists))

        pyplot.legend()
        if title:
            pyplot.title(title)
        pyplot.show()

    elif args.proto:
        print(os.listdir(os.curdir))
        print(os.path.abspath(os.curdir))
        for a, b in os.environ.items():
            if "repos" in b:
                print()
                print(a)
                print(b)
                print()
        with open(args.proto, "rb") as proto_in:
            dhist = DensityMap_pb2.ParseFromString(proto_in.read())
            print(dhist)
