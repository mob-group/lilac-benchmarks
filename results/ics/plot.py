#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_style('ticks')
sns.set_palette(sns.cubehelix_palette(3, rot=60/360, dark=0.2, light=0.7))

# Compute the size of a figure based on the current columnwidth
def fig_size(w, h):
    colwidth = 240.94499
    colwidth_in = colwidth / 72.26999
    return (w * colwidth_in, h * colwidth_in)

# Turn an internal representation of an implementation name into a friendly
# printable one.
def impl_map(platform, impl):
    if platform == 'michel' and impl == 'opencl10':
        return 'clSPARSE (eGPU)'
    elif platform == 'firuza' and impl == 'opencl00':
        return 'clSPARSE'
    
    return {
        'mkl': 'MKL',
        'gpu': 'cuSPARSE',
        'sparsex': 'SparseX',
        'native': 'Native'
    }[impl]

# Similarly, pretty-print the name of each benchmark / application
def bench_map(bench):
    return {
        'pfold': 'PFold',
        'ngt': 'NGT',
        'PageRank': 'PageRank',
        'bfs': 'BFS',
        'NPB': 'NPB',
        'parboil-spmv': 'Parboil SPMV',
        'Netlib-C': 'Netlib C',
        'Netlib-F': 'Netlib Fortran'
    }[bench]

# ... and the platform
def platform_map(platform):
    return {
        'monaco': 'Intel-0',
        'firuza': 'Intel-1',
        'michel': 'AMD'
    }[platform]

def get_data(infile):
    return pd.read_csv(infile)

def next_tick(val, ticks):
    return np.ceil((1/ticks) * val) / (1/ticks)

def baseline_impl(df, benches, tick_size=0.5):
    platforms = df.platform.unique()

    fig, axes = plt.subplots(1, len(benches), sharey='row', figsize=fig_size(2.1, 0.67))
    max_y = df[df['benchmark'].isin(benches)]['speedup'].max()

    for plot_num, (bench, ax) in enumerate(zip(benches, axes)):
        i = 0
        
        if plot_num == 0:
            ax.set_ylabel('Speedup (×)')
        
        results = df.query("benchmark=='{}'".format(bench))
        
        labs = []
        bars = []
        legend = []
        
        rows = [r for r in results.iterrows()]
        rows.sort(key=lambda r: platform_map(r[1].platform))

        ax.set_xlim(-0.6, 2.6)
        ax.set_ylim(0.8, next_tick(max_y, tick_size))
        ax.set_yticks(np.arange(1, next_tick(max_y, tick_size) + tick_size, tick_size))

        for row in rows:
            y_val = row[1].speedup
            bars.append(ax.bar(i, y_val, width=0.95))
            legend.append(platform_map(row[1].platform))
            labs.append(impl_map(row[1].platform, row[1].implementation))
            
            threshold = 0.85 * next_tick(max_y, tick_size)
            valign = 'top' if y_val > threshold else 'bottom'
            offset = -max_y/30 if y_val > threshold else max_y/30
            color = 'white' if y_val > threshold else 'black'

            ax.text(i, y_val + offset, impl_map(row[1].platform,
                row[1].implementation), rotation=90, verticalalignment=valign,
                horizontalalignment='center', color=color, fontsize=8)
            i += 1
                    
        ax.set_xticks([])
        ax.set_xticklabels([])
        
        ax.axhline(y=1, color='white', linestyle=':')
        
        ax.set_title(bench_map(bench))
        
    ax.legend(bars, legend, bbox_to_anchor=(1.5,0.5), loc='center')
      
    fig.tight_layout()
    plt.subplots_adjust(wspace=0.2)

    sns.despine(fig)
    return fig

def baseline(df):
    benches = ['pfold', 'ngt', 'PageRank', 'bfs']
    return baseline_impl(df, benches)

def baseline_bench(df):
    benches = ['NPB', 'parboil-spmv', 'Netlib-C', 'Netlib-F']
    return baseline_impl(df, benches, tick_size=4.0)

def marshall(df):
    def plot_bar(x, bench, platform, impl, ax):
        row = df.query('benchmark==@bench and platform==@platform').head(1)[impl]
        return ax.bar(x, row, width=0.95)

    fig, (pfold, ngt, pr, bfs) = plt.subplots(1, 4, figsize=fig_size(2.1, 0.67), sharey=True)
        
    pfold.set_title('PFold')
    ngt.set_title('NGT')
    pr.set_title('PageRank')
    bfs.set_title('BFS')
    
    groups = {
        pfold: ('pfold', [
            ('michel', 'opencl10-slow'),
            ('monaco', 'mkl-slow'),
            ('firuza', 'mkl-slow')
        ]),
        ngt: ('ngt', [
            ('michel', 'opencl10-slow'),
            ('monaco', 'mkl-slow'),
            ('firuza', 'mkl-slow')
        ]),
        pr: ('PageRank', [
            ('michel', 'gpu-slow'),
            ('monaco', 'mkl-slow'),
            ('firuza', 'opencl00-slow')
        ]),
        bfs: ('bfs', [
            ('michel', 'sparsex-slow'),
            ('monaco', 'mkl-slow'),
            ('firuza', 'mkl-slow')
        ])
    }
    
    for ax in groups:
        ax.set_ylim(0.8, 4.0)
        ax.set_xlim(-0.6, 2.6)
        
        ax.set_xticks([])
        ax.set_xticklabels([])
        
        name, bars = groups[ax]
        bar_list = []
        for i, bar in enumerate(bars):
            bar_list.append(plot_bar(i, name, *bar, ax))

        ax.axhline(y=1, color='white', linestyle=':')

    # add legend when all bars in place
    ax.legend(bar_list, '', bbox_to_anchor=(1.5,0.5), loc='center')
    
    fig.tight_layout()
    sns.despine(fig)
    return fig

def expert(df):
    sns.set_palette(sns.cubehelix_palette(2, rot=60/360, dark=0.2, light=0.7-(0.5/2)))

    fig, ax = plt.subplots(1, 1, figsize=fig_size(0.9, 0.9))
    
    b1 = ax.bar(0, df.loc[1]['openmp-expert'], width=0.95)
    b2 = ax.bar(1, df.loc[0]['openmp-expert'], width=0.95)
    
    ax.bar(2.5, df.loc[1]['opencl-expert'], width=0.95)
    ax.bar(3.5, df.loc[0]['opencl-expert'], width=0.95)

    ax.bar(5, df.loc[3]['opencl-expert'], width=0.95)
    ax.bar(6, df.loc[2]['opencl-expert'], width=0.95)

    ax.set_xticks([0.5, 3, 5.5])
    ax.set_xticklabels(['NPB\nOpenMP', 'NPB\nOpenCL', 'Parboil\nOpenCL'])
    
    ax.set_yticks([0, 0.5, 1])
    ax.set_ylabel('LiLAC Performance (×)')
    
    ax.set_title('LiLAC vs. Expert')
    ax.legend(['Intel-0', 'Intel-1'], loc='upper right')
    
    fig.tight_layout()
    sns.despine(fig)
    return fig

plot_choices = { p.__name__ : p for p in [
    baseline,
    baseline_bench,
    expert,
    marshall
]}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Baseline comparison plots for ICS')
    parser.add_argument('plot', choices=plot_choices)
    parser.add_argument('data', type=str, help='Data file to tidy')
    args = parser.parse_args()

    plot_f = plot_choices[args.plot]
    fig = plot_f(get_data(args.data))
    fig.savefig('{}.pdf'.format(args.plot))
