#!/usr/bin/env python

import matplotlib
import statistics

from typing import Dict, List

class ResultSet:
    def __init__(self, name, variants):
        self.name = name
        self.__variants = variants

    def means(self):
        ret = {}
        for var_name in self.__variants:
            ret[var_name] = statistics.mean(self.__variants[var_name])
        return ret

if __name__ == "__main__":
    x = ResultSet('xbm', {'native' : [1.0, 3.5]})
    print(x.means())
