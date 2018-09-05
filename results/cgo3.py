#!/usr/bin/env python

import csv
import matplotlib.pyplot as plt
import sys

# Merge data by benchmark, then for each pair plot a pair of bars with opposing
# styles
# add label below

def merged_data(data):
    ret = {}
    for row in data:
        if row['name'] not in ret:
            ret[row['name']] = {}
        ret[row['name']][row['platform']] = row['speedup']
    return ret

if __name__ == "__main__":
    with open(sys.argv[1]) as csvfile:
        reader = csv.DictReader(csvfile)
        data = merged_data(reader)
    print(data)
