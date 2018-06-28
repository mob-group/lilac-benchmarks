#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

from sklearn import svm
from sklearn_porter import Porter

PLOT_NNZ_MAX=250000
PLOT_ROW_MAX=50000

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

def grid(features):
    rows, nnz = [np.array(l) for l in zip(*features)]
    r_min, r_max = rows.min(), PLOT_ROW_MAX
    n_min, n_max = nnz.min(), PLOT_NNZ_MAX
    xx, yy = np.meshgrid(np.arange(n_min, n_max, 100),
                         np.arange(r_min, r_max, 100))
    return xx, yy

def plot_contours(ax, clf, xx, yy, **params):
    """Plot the decision boundaries for a classifier.

    Parameters
    ----------
    ax: matplotlib axes object
    clf: a classifier
    xx: meshgrid ndarray
    yy: meshgrid ndarray
    params: dictionary of params to pass to contourf, optional
    """
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    out = ax.contourf(xx, yy, Z, **params)
    return out

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('train', type=str, help='training data file')
    parser.add_argument('--test', type=str, help='test data file to validate model')
    parser.add_argument('--code', action='store_true', help='generate C code for model')
    parser.add_argument('--plot', action='store_true', help='plot trained model decision contour')
    args = parser.parse_args()

    frame = pd.read_csv(args.train, sep=' ')
    train_x, train_y = get_data(frame)

    clf = svm.SVC(kernel='linear', gamma=0.001)
    clf.fit(train_x, train_y)

    if args.plot:
        fig, ax = plt.subplots()
        xx, yy = grid(train_x)
        plot_contours(ax, clf, xx, yy, cmap=plt.cm.coolwarm, alpha=0.8)
        
        rows, nnz = [np.array(l) for l in zip(*train_x)]
        ax.set_xlim(xx.min(), PLOT_NNZ_MAX)
        ax.set_ylim(yy.min(), PLOT_ROW_MAX)
        ax.scatter(nnz, rows, c=train_y, cmap=plt.cm.coolwarm, s=20, edgecolors='k')

        ax.set_xlabel("Non-zero elements")
        ax.set_ylabel("Rows")
        plt.show()

    if args.test is not None:
        test_frame = pd.read_csv(args.test, sep=' ')
        test_x, test_y = get_data(test_frame)
        preds = clf.predict(test_x)
        results = (np.array(preds) == np.array(test_y))
        acc = 100 * sum(results) / len(results)
        print(f"Test Accuracy: {acc:.3f}%")

    if args.code:
        p = Porter(clf, language='C')
        print(p.export(embed_data=True))
