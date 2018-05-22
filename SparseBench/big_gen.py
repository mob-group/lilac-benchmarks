#!/usr/bin/env python3
import argparse
import random
import sys

def parse_fixed(row, cols, types, widths):
    if not isinstance(types, list):
        types = [types] * cols
    if not isinstance(widths, list):
        widths = [widths] * cols
    values = []
    cursor = 0
    for i in range(cols):
        values.append(types[i](row[cursor:cursor+widths[i]]))
        cursor += widths[i]
    return values

def write_fixed(values, widths):
    pass

def read_nnz(data):
    line = next(data)
    n, nnz = parse_fixed(line, 2, int, 12)
    return n, nnz

def read_row_ptr(data, n):
    row_ptr = []
    for i in range(n+1):
        line = next(data)
        item, = parse_fixed(line, 1, int, 12)
        row_ptr.append(item)
    return row_ptr

def read_a_col_ind(data, nnz):
    a = []
    col_ind = []
    for line in data:
        ind, val = parse_fixed(line, 2, [int, float], [12, 26])
        a.append(val)
        col_ind.append(ind)
    assert(len(a) == nnz)
    assert(len(col_ind) == nnz)
    return a, col_ind

def read_crs(filename):
    with open(filename, 'r') as data:
        n, nnz = read_nnz(data)
        row_ptr = read_row_ptr(data, n)
        a, col_ind = read_a_col_ind(data, nnz)
    return n, nnz, row_ptr, a, col_ind

def write_crs(f, n, nnz, row_ptr, a, col_ind):
    print("{:12}{:12}".format(n, nnz), file=f)
    for ptr in row_ptr:
        print("{:12}".format(ptr), file=f)
    for ind, val in zip(col_ind, a):
        print("{:12} {:20.17f}".format(ind, val), file=f)

def random_crs(size):
    n = size**3
    row_ptr = [1]
    row_counts = []
    for i in range(n):
        count = int(random.gauss(5, 4))
        if count < 1: count = 1
        if count > n: count = n
        row_ptr.append(count + row_ptr[-1])
        row_counts.append(count)
    nnz = row_ptr[-1] - 1
    a = []
    col_ind = []
    for (ind, i) in enumerate(row_counts):
        ind = ind+1
        sample = random.sample(range(1,n+1), i)
        if ind not in sample:
            sample.append(ind)
        cols = sorted(sample)
        for col in cols:
            col_ind.append(col)
            val = random.gauss(0, 2)
            if True: a.append(abs(val))
            else: a.append(val)
    return n, nnz, row_ptr, a, col_ind

def write_random(filename, size):
    if filename is None:
        write_crs(sys.stdout, *random_crs(size))
    else:
        with open(filename, 'w') as f:
            write_crs(f, *random_crs(size))

def read_and_dump(filename):
    crs = read_crs(filename)
    write_crs(sys.stdout, *crs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['read', 'write'])
    parser.add_argument('--filename')
    parser.add_argument('--size', type=int)
    args = parser.parse_args()
    if args.mode == 'write' and args.size is not None:
        write_random(args.filename, args.size)
    elif args.mode == 'read':
        read_and_dump(args.filename)
    else:
        print("Must specify a size if writing a random matrix", file=sys.stderr)
        exit(1)
