#!/usr/bin/env python3

import sys

def next(t):
    return {
        "native" : "opencl",
        "opencl" : "clgpu",
        "clgpu" : "native"
    }[t]

if __name__ == "__main__":
    i = 0
    with open(sys.argv[1]) as in_f, open(sys.argv[2], 'w') as big_f, open(sys.argv[3], 'w') as small_f:
        expected = 'clgpu'
        for line in in_f:
            method = line.split(' ')[-1][:-1]
            if i != 0 and method != expected:
                print(line)
                exit(1)
            if i == 0:
                big_f.write(line)
                small_f.write(line)
            elif (i - 1) % 30 < 3:
                small_f.write(line)
            else:
                big_f.write(line)
            i += 1
            expected = next(expected)
