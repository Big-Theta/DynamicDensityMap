#!/usr/bin/env python3

from matplotlib import pyplot
import numpy as np
import sys
import json


data = sys.stdin.read()
histogram = json.loads(data)

counts = np.array(histogram["counts"])
total_count = sum(counts)
weights = np.array([count / total_count for count in counts])
bins = np.array(histogram["bounds"])
widths = bins[1:] - bins[:-1]
heights = weights.astype(np.float) / widths
pyplot.fill_between(
        bins.repeat(2)[1:-1], heights.repeat(2),
        color="steelblue", alpha=1.0)
pyplot.show()
