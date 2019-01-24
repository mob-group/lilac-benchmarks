#!/usr/bin/env python

import numpy as np
from numpy.random import normal

MACHINES = {
    'monaco'    : 1.0,
    'michel'    : 0.4,
    'avus'      : 1.3,
    'bob'       : 0.6, 
    'spa'       : 1.4
}

IMPLS = {
    'monaco'    : ['native', 'mkl', 'gpu', 'opencl00', 'sparsex'],
    'michel'    : ['native', 'gpu', 'opencl00', 'opencl01', 'sparsex'],
    'avus'      : ['native', 'mkl', 'gpu', 'opencl10', 'sparsex'],
    'bob'       : ['native', 'opencl00', 'sparsex'],
    'spa'       : ['native', 'mkl', 'gpu', 'opencl00', 'sparsex'],
}

SPEEDUPS = {
    'native' : 1.0, 'mkl': 2.5, 'gpu': 1.2, 
    'opencl00': 0.8, 'opencl01': 0.9, 'opencl10': 1.1,
    'sparsex' : 0.4
}

BENCHES = {
    'bfs'           : [
        'roadNet-CA','kron','in-2004','ljournal-2008',
        'USA-road-d.E','higgs-twitter','Si5H12',
        'com-Youtube','erdos','eu-2005'
    ],
    'NPB'           : [
        'S', 'W', 'A', 'B', 'C', 'D'
    ],
    'PageRank'      : [
        'roadNet-CA','kron','in-2004','ljournal-2008',
        'USA-road-d.E','higgs-twitter','Si5H12',
        'com-Youtube','erdos','eu-2005'
    ],
    'parboil-spmv'  : [
        'small', 'medium', 'large'
    ],
    'Netlib-C'      : [
        '40', '60', '80', '100', '120', '140', '160'
    ],
    'Netlib-F'      : [
        '40', '60', '80', '100', '120', '140', '160'
    ],
    'pfold'    : [
        '0.small', '1.small', '2.small',
        '0.large', '1.large', '2.large',
    ],
    'ngt'      : [
        '0.small', '1.small', '2.small',
        '0.large', '1.large', '2.large',
    ]
}

def random_times(mach, impl, bench, name):
    base_time = 10.0
    perf_time = base_time * MACHINES[mach] * SPEEDUPS[impl]
    difficulty_factor = BENCHES[bench].index(name) ** 1.3
    mean_time = perf_time * difficulty_factor
    return np.abs(normal(loc=mean_time, size=5))

if __name__ == "__main__":
    for mach in MACHINES:
        for bench in BENCHES:
            for impl in IMPLS[mach]:
                for name in BENCHES[bench]:
                    ts = random_times(mach, impl, bench, name)
                    print("{},{},{},{},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f}".format(
                        mach, bench, impl, name, *ts))
