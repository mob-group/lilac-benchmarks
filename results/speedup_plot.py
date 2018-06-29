#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

from sklearn import svm
from sklearn.preprocessing import StandardScaler
from sklearn_porter import Porter

def speedups(sub_f, method):
    nat = sub_f.query('method == "native"')
    meth = sub_f.query('method == @method')
    
    speedup = float(nat.time) / float(meth.time)
    return int(nat.nnz), speedup

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('train', type=str, help='training data file')
    parser.add_argument('method', type=str, help='acceleration method to compare')
    args = parser.parse_args()

    frame = pd.read_csv(args.train, sep=' ')
    data = np.array([speedups(g, args.method) for _, g in frame.groupby(by='matrix')]).T

    slows = [x for x in data[1] if x < 1]
    print(sum(slows) / len(slows))

    plt.scatter(data[0], data[1], s=20)
    plt.show()
