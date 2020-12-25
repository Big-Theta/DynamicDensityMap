#!/usr/bin/env python3

from matplotlib import pyplot
from typing import Dict, List
import numpy as np
import sys
import json


def maybe_histogram(data: Dict) -> bool:
    if (isinstance(data.get("bounds"), list)
        and isinstance(data.get("counts"), list)
        and len(data["bounds"]) == len(data["counts"]) + 1
    ):
        return True
    return False


def extract_histograms(data: List[str]) -> List[dict]:
    candidates = []
    open_brace = None
    nest_count = 0
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


if __name__ == "__main__":
    data = sys.stdin.read()
    print(data)

    hists = extract_histograms(data)

    for hist in hists:
        if maybe_histogram(hist):
            label = hist.get("label", "")
            if len(hists) > 1 and "title" in hist:
                label = hist["title"] + " -- " + label

            prepare_render(hist, label=label, alpha=1.0 / len(hists))

    pyplot.show()

