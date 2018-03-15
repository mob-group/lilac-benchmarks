#!/usr/bin/env python3

import csv
import math
import matplotlib.pyplot as plt
import statistics
import sys

def row_to_stats(row):
    new_row = []
    new_row.append(row[0])
    parsed = list(map(float, row[1:]))
    new_row.append(len(parsed))
    new_row.append(statistics.mean(parsed))
    new_row.append(statistics.stdev(parsed))
    return new_row

if __name__ == "__main__":
    in_file = sys.argv[1]
    with open(in_file, "r") as csvfile:
        reader = csv.reader(csvfile)
        writer = csv.writer(sys.stdout)
        writer.writerow(['threads', 'samples', 'time', 'stdev'])
        for row in reader:
            writer.writerow(row_to_stats(row))
