#!/usr/bin/env python3

from matplotlib import pyplot
from typing import Dict, List
import numpy as np
import sys
import json


def maybe_histogram(data: Dict) -> bool:
    if (isinstance(data.get("bounds"), list)
        and isinstance(data.get("counts"), list)
        and len(data["bounds"]) == len(data["counts"]) - 1
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

    # Make sure that only json objects that either have the keys "bounds" and
    # "counts" or else have children with those keys.
    filtered = []
    for candidate in cleaned:
        if maybe_histogram(candidate):
            filtered.append(candidate)
        else:
            for value in candidate.values():
                if not isinstance(value, dict) or not maybe_histogram(value):
                    break
            else:
                filtered.append(candidate)

    return filtered


def prepare_render(histogram: Dict, title: str = "", alpha: float = 1.0):
    print(histogram)
    counts = np.array(histogram["counts"])
    total_count = sum(counts)
    weights = np.array([count / total_count for count in counts])
    bins = np.array(histogram["bounds"])
    widths = bins[1:] - bins[:-1]
    print(widths)
    heights = weights.astype(np.float) / widths
    print(heights)
    pyplot.fill_between(
            bins.repeat(2)[1:-1], heights.repeat(2),
            color="steelblue", alpha=alpha)


if __name__ == "__main__":
    data = sys.stdin.read()

    hists = extract_histograms(data)

    for hist in hists:
        prepare_render(hist, alpha=1.0 / len(hists))

    pyplot.show()

