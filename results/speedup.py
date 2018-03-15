#!/usr/bin/env python3

import csv
import math
import statistics
import sys

# Computes the speedup and propagated error in the speedup, using the variance
# formula to propagate errors. r is the raw runtime, t the threaded runtime, and
# s_{r,t} are their respective standard deviations.
def speedup(r, s_r, t, s_t):
    [r, s_r, t, s_t] = map(float, [r, s_r, t, s_t])
    t1 = ((r**2) * (s_t**2)) / (t ** 4)
    t2 = (s_r**2) / (t**2)
    return (r / t, math.sqrt(t1 + t2))

if __name__ == "__main__":
    rows = []
    reader = csv.reader(sys.stdin)
    reader.__next__()
    for row in reader:
        rows.append(row)
    writer = csv.writer(sys.stdout)
    writer.writerow(['threads', 'samples', 'speedup', 'stdev'])
    raw_time, raw_stdev = rows[-1][2:4]
    for row in rows[:-1]:
        writer.writerow([*row[:2], *speedup(raw_time, raw_stdev, *row[2:4])])
