#!/usr/bin/env python3

import os
import sys

def all_times(path):
    with open(path) as f:
        all_times = {}
        for line in f:
            fields = line.split(' ')
            name = fields[0]
            times = {
                'base' : float(fields[1]),
                'integrated': float(fields[2]),
                'gpu': float(fields[3])
            }
            all_times[name] = times
        return all_times

def speedups(all_times):
    speedups = {}
    for matrix in all_times:
        base_time = all_times[matrix]['base']
        speedups[matrix] = {
            x : base_time / all_times[matrix][x] for x in all_times[matrix]
        }
    return speedups

if __name__ == "__main__":
    ats = all_times(sys.argv[1])
    speeds = speedups(ats)
    print(speeds)
