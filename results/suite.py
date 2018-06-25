#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

from sklearn import linear_model

def speedup(sub_f):
    base_time = sub_f.query('method == "native"').time
    return sub_f.assign(
        density=lambda x: x.nnz / (x.rows * x.cols),
        speedup=lambda x: float(base_time) / x.time
    )

if __name__ == "__main__":
    frame = pd.read_csv(sys.argv[1], sep=' ')
    results = pd.concat([speedup(sub_f) for _, sub_f in frame.groupby(by='matrix')])
    gpu = results.query('method == "gpu"')
    integrated = results.query('method == "integrated"')
    reg = linear_model.LinearRegression()
    reg.fit(np.array(gpu.nnz).reshape(-1,1), gpu.speedup)
    ys = reg.coef_ * gpu.nnz + reg.intercept_
    plt.plot(gpu.nnz, ys)
    plt.scatter(gpu.nnz, gpu.speedup)
    plt.show()
