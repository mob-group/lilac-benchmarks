#! /usr/bin/env python3

import pandas as pd
import sys

if __name__ == "__main__":
    frame = pd.read_csv(sys.argv[1])
    links = frame[frame.nnz < 10000000].link
    for row in links:
        print(row)
