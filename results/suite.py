#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

from sklearn_porter import Porter

from sklearn import svm

def name_to_output(name):
    return {
        "native" : 0,
        "opencl" : 1,
        "clgpu" : 2
    }[name]

def output_to_name(output):
    return ["native", "opencl", "clgpu"][output]

def train_data(sub_f):
    nat = sub_f.query('method == "native"')
    out = sub_f.ix[sub_f.time.idxmin()].method
    return (int(nat.rows), int(nat.nnz)), name_to_output(out)

if __name__ == "__main__":
    frame = pd.read_csv(sys.argv[1], sep=' ')
    data = [train_data(sub_f) for _, sub_f in frame.groupby(by='matrix')]

    train_x = [d[0] for d in data]
    train_y = [d[1] for d in data]

    clf = svm.SVC(gamma=0.001)
    clf.fit(train_x, train_y)
    
    preds = clf.predict(train_x)

    p = Porter(clf, language='C')
    print(p.export(embed_data=True))
