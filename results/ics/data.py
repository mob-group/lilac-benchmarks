import collections
import csv

class Row(object):
    def __init__(self, platform, bench, impl, name, times):
        self.platform = platform
        self.benchmark = bench
        self.implementation = impl
        self.name = name
        self.times = times

    @classmethod
    def from_raw(cls, raw):
        if len(raw) < 5:
            return None
        raw_times = [t for t in raw[4:] if t != '']
        if not raw_times:
            return None
        return cls(*raw[:4], [float(t) for t in raw_times])

    def mean(self):
        return sum(self.times) / float(len(self.times))

    def __repr__(self):
        return "Row(platform={}, benchmark={}, \
implementation={}, name={}, \
times={})".format(
        self.platform, self.benchmark, 
        self.implementation, self.name, self.times)

    def __str__(self):
        return "{:<16s}\t{:<16s}\t{:<16s}".format(
                    "{}[{}]".format(self.platform, self.implementation),
                    "{}:{}".format(self.benchmark, self.name),
                    "mean={:.2f}".format(self.mean())
                )

def query(dataset, platform=None, benchmark=None, 
            implementation=None, name=None):
    results = dataset

    if platform is not None:
        results = [r for r in results if r.platform == platform]

    if benchmark is not None:
        results = [r for r in results if r.benchmark == benchmark]

    if implementation is not None:
        results = [r for r in results if r.implementation == implementation]

    if name is not None:
        results = [r for r in results if r.name == name]

    return results

def list_values(dataset):
    cls = collections.namedtuple(
            'Dataset', ['platforms', 'benchmarks', 'implementations', 'names'])
    return cls(
        platforms={r.platform for r in dataset},
        benchmarks={r.benchmark for r in dataset},
        implementations={r.implementation for r in dataset},
        names={r.name for r in dataset}
    )

def read_rows(filenames):
    ret = []
    for f in filenames:
        with open(f, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                row_obj = Row.from_raw(row)
                if row_obj is not None:
                    ret.append(row_obj)
    return ret
