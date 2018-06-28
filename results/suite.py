#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

from sklearn import svm
from sklearn_porter import Porter

def name_to_output(name):
    return {
        "native" : 0,
        "clgpu" : 1,
        "opencl" : 2,
    }[name]

def output_to_name(output):
    return ["native", "opencl", "clgpu"][output]

def train_data(sub_f):
    nat = sub_f.query('method == "native"')
    out = sub_f.ix[sub_f.time.idxmin()].method
    return (int(nat.rows), int(nat.nnz)), name_to_output(out)

def get_data(frame):
    data = [train_data(sub_f) for _, sub_f in frame.groupby(by='matrix')]

    xs = [d[0] for d in data]
    ys= [d[1] for d in data]

    return xs, ys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('train', type=str, help='training data file')
    parser.add_argument('--test', type=str, help='optional test data to validate model')
    args = parser.parse_args()

    frame = pd.read_csv(args.train, sep=' ')
    train_x, train_y = get_data(frame)

    clf = svm.SVC(kernel='linear', gamma=0.001)
    clf.fit(train_x, train_y)

    if args.test is not None:
        test_frame = pd.read_csv(args.test, sep=' ')
        test_x, test_y = get_data(test_frame)
        preds = clf.predict(test_x)
        results = (np.array(preds) == np.array(test_y))
        acc = 100 * sum(results) / len(results)
        print(f"Test Accuracy: {acc:.3f}%")
    else:
        p = Porter(clf, language='C')
        print(p.export(embed_data=True))
