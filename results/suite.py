#!/usr/bin/env python3

import os
import sys

class Result:
    def __init__(self, line):
        items = line.split(' ')
        self.path = items[0]
        self.time = float(items[1])
        self.rows = int(items[2])
        self.cols = int(items[3])
        self.nnz = int(items[4])
        self.iters = int(items[5])

    def __repr__(self):
        return "Result({}, {}, {}, {}, {})".format(self.path, self.time,
                self.rows, self.cols, self.nnz, self.iters)

def file_results(path):
    with open(path, 'r') as f:
        return [Result(line) for line in f]

def matrix_name(result):
    return os.path.basename(result.path).split('.')[0]

if __name__ == "__main__":
    baseline = sys.argv[1]
    variants = sys.argv[2:]
    b_results = file_results(baseline)
    by_name = {}
    for r in b_results:
        by_name[matrix_name(r)] = [r]
    
    for var in variants:
        v_results = file_results(var)
        for r in v_results:
            by_name[matrix_name(r)].append(r)
