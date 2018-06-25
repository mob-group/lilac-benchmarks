#!/usr/bin/env python3

import sys

from bs4 import BeautifulSoup

def comma_int(s):
    return int(s.replace(',', ''))

def header():
    return "name,group,kind,nnz,rows,cols,link"

class RowData:
    def __init__(self, row):
        link = row.find(text='Matrix Market ')
        self.link = link.parent['href'] if link else None

        n_rows = row.find(class_='column-num_rows')
        self.num_rows = comma_int(n_rows.text) if n_rows else None

        n_cols = row.find(class_='column-num_cols')
        self.num_cols = comma_int(n_cols.text) if n_cols else None

        n_nonzeros = row.find(class_='column-nonzeros')
        self.num_nonzeros = comma_int(n_nonzeros.text) if n_nonzeros else None

        group = row.find(class_='column-group')
        self.group = group.text if group else None

        kind = row.find(class_='column-kind')
        self.kind = kind.text if kind else None

        name = row.find(class_='column-name')
        self.name = name.text if name else None

    def __repr__(self):
        return f"RowData({self.group}, {self.kind}, {self.num_nonzeros}, {self.num_rows}, {self.num_cols}, {self.link})"

    def file_row(self):
        return f"{self.name},{self.group},{self.kind},{self.num_nonzeros},{self.num_rows},{self.num_cols},{self.link}"

if __name__ == "__main__":
    print(header())
    with open(sys.argv[1]) as fp:
        soup = BeautifulSoup(fp, 'html.parser')
        table = soup.find(id='matrices')
        for row in table.findChildren('tr'):
            print(RowData(row).file_row())
