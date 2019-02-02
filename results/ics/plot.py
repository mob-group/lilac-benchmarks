#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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
    impl = impl.replace('-slow', '')

    if platform == 'michel' and impl == 'opencl10':
        return 'eGPU'
    elif platform == 'firuza' and impl == 'opencl00':
        return 'eGPU'
    elif platform == 'michel' and impl == 'opencl00':
        return 'iGPU'
    elif platform == 'monaco' and impl == 'opencl00':
        return 'eGPU'

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
        'NPB': 'NPB-CG',
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

def baseline_impl(df, benches, tick_size=0.5, ticks=None):
    platforms = df.platform.unique()

    fig, axes = plt.subplots(1, len(benches), sharey='row',
            figsize=fig_size(2.1, 0.5))
    max_y = df[df['benchmark'].isin(benches)]['speedup'].max()

    for plot_num, (bench, ax) in enumerate(zip(benches, axes)):
        i = 0
        
        if plot_num == 0:
            ax.set_ylabel('Speedup (×)', fontsize=8)
        
        results = df.query("benchmark=='{}'".format(bench))
        
        labs = []
        bars = []
        legend = []
        
        rows = [r for r in results.iterrows()]
        rows.sort(key=lambda r: platform_map(r[1].platform))

        ax.set_xlim(-0.6, 2.6)
        ax.set_ylim(0.8, max(ticks) if ticks is not None else next_tick(max_y, tick_size))

        if ticks is None:
            ax.set_yticks(np.arange(1, next_tick(max_y, tick_size) + tick_size, tick_size))
        else:
            ax.set_yticks(ticks)
            ax.tick_params(axis='y', labelsize=6)

        for row in rows:
            y_val = row[1].speedup
            bars.append(ax.bar(i, y_val, width=0.95))
            legend.append(platform_map(row[1].platform))
            labs.append(impl_map(row[1].platform, row[1].implementation))
            
            threshold = 0.85 * next_tick(max_y, tick_size)
            valign = 'top' if y_val > threshold else 'bottom'
            offset = -max_y/20 if y_val > threshold else max_y/20
            color = 'white' if y_val > threshold else 'black'

            ax.text(i, y_val + offset, impl_map(row[1].platform,
                row[1].implementation), rotation=90, verticalalignment=valign,
                horizontalalignment='center', color=color, fontsize=6)
            i += 1
                    
        ax.set_xticks([])
        ax.set_xticklabels([])
        
        ax.axhline(y=1, color='black', lw=0.8)
        
        ax.set_title(bench_map(bench), fontsize=10)
        
    ax.legend(bars, legend, bbox_to_anchor=(1.25,0.5), loc='center', fontsize=6)
      
    fig.tight_layout()
    plt.subplots_adjust(wspace=0.2)

    sns.despine(fig)
    return fig

def baseline(df):
    benches = ['pfold', 'ngt', 'parboil-spmv', 'bfs']
    return baseline_impl(df, benches, tick_size=0.1, ticks=np.arange(1.0, 2.1,
        0.1))

def baseline_bench(df):
    benches = ['NPB', 'PageRank', 'Netlib-C', 'Netlib-F']
    return baseline_impl(df, benches, tick_size=1, ticks=np.arange(1, 13, 1))

def marshall(df):
    def plat_color(platform):
        plats = ['michel', 'monaco', 'firuza']
        pal = sns.cubehelix_palette(3, rot=60/360, dark=0.2, light=0.7)
        return pal[plats.index(platform)]

    def plot_bar(x, bench, platform, impl, ax):
        row = df.query('benchmark==@bench and platform==@platform').head(1)[impl]
        orig = row.item()
        row = min(row.item(), 4.0)

        threshold = 0.85 * 4.0
        valign = 'top' if row > threshold else 'bottom'
        offset = -4.0/30 if row > threshold else 4.0/30
        color = 'white' if row > threshold else 'black'

        text = impl_map(platform, impl)
        ax.text(i, row + offset, text, rotation=90, verticalalignment=valign,
            horizontalalignment='center', color=color, fontsize=8)

        if orig > row:
            ax.text(i, row + 0.02, "{:.1f}×".format(orig), ha='center',
                    va='bottom', fontsize=6)
        return ax.bar(x, row, width=0.95, color=plat_color(platform))

    fig, ax = plt.subplots(1, 1, figsize=fig_size(1, 0.67), sharey=True)
    ax.set_title('Speedup vs. Naïve Memory Transfers', fontsize=10, y=1.05)
    
    groups = [
        ('pfold', [
            ('michel', 'opencl10-slow'),
            # ('monaco', 'mkl-slow'),
            # ('firuza', 'mkl-slow')
        ]),
        ('ngt', [
            ('michel', 'opencl10-slow'),
            # ('monaco', 'mkl-slow'),
            # ('firuza', 'mkl-slow')
        ]),
        ('PageRank', [
            ('michel', 'gpu-slow'),
            # ('monaco', 'mkl-slow'),
            ('firuza', 'opencl00-slow')
        ]),
        ('bfs', [
            ('michel', 'sparsex-slow'),
            # ('monaco', 'mkl-slow'),
            # ('firuza', 'mkl-slow')
        ])
    ]
    
    ax.set_ylim(0.8, 4.0)
    ax.set_yticks(np.arange(1.0, 4.1, 0.5))
    ax.set_ylabel('Speedup (×)', fontsize=8)
    
    i = 0
    ticks = []
    labels = []
    leg = {}
    for name, bars in groups:
        ticks.append(i + len(bars) / 2 - 0.5)
        labels.append(bench_map(name))

        bar_list = []
        for bar in bars:
            leg[platform_map(bar[0])] = plot_bar(i, name, *bar, ax)
            i += 1
        i += 0.5

    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=6)

    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)

    ax.axhline(y=1, color='black', lw=0.8)

    ax.legend(leg.values(), leg.keys(), bbox_to_anchor=(0.02,1), loc='upper left', fontsize=6)
    
    fig.tight_layout()
    sns.despine(fig)
    return fig

def expert(df):
    sns.set_palette(sns.cubehelix_palette(2, rot=60/360, dark=0.2, light=0.7))

    fig, ax = plt.subplots(1, 1, figsize=fig_size(0.9, 0.5))

    # ax.bar(0, df.loc[0]['native'], width=0.95)
    ax.bar(0, df.loc[0]['opencl00'], width=0.95)
    ax.bar(1, df.loc[0]['opencl-expert'], width=0.95)

    # ax.bar(3.5, df.loc[1]['native'], width=0.95)
    ax.bar(2.5, df.loc[1]['opencl00'], width=0.95)
    ax.bar(3.5, df.loc[1]['opencl-expert'], width=0.95)
    
    ax.set_xticks([0.5, 3])
    ax.set_xticklabels(['NPB-CG', 'Parboil'], fontsize=8)
    
    ax.set_yticks(np.arange(0, 1.1, 0.2))
    ax.tick_params(axis='y', labelsize=6)
    ax.set_ylabel('Performance (×)', fontsize=8)
    
    ax.set_title('LiLAC vs. Expert Implementation', fontsize=10)
    ax.legend(['LiLAC', 'Expert'], loc='upper left', fontsize=6,
            bbox_to_anchor=(0, 0.95))

    ax.axhline(y=1, color='black', lw=0.8)
    
    fig.tight_layout()
    sns.despine(fig)
    return fig

def distribution(df):
    def marker(name):
        markers = ['.', '+', 'x', '1', '*', 's']
        return dict(zip(impls, markers))[name]

    def color(name):
        pal = sns.cubehelix_palette(6, rot=60/360, dark=0.2, light=0.7)
        return pal[list(impls).index(name)]

    def new_impl_map(plat, impl):
        if plat == 'michel' and impl == 'opencl00':
            return 'CL iGPU'
        elif plat == 'michel' and impl == 'opencl10':
            return 'CL eGPU'

        return {
            'gpu': 'cuSPARSE',
            'mkl': 'MKL',
            'opencl00': 'CL eGPU',
            'sparsex': 'SparseX',
            'opencl10': 'CL iGPU'
        }[impl]
    
    sizes = ['S', 'W', 'A', 'B', 'C', 'D']

    def row_cmp(row):
        return platform_map(row[0][2]), sizes.index(row[0][3])

    speeds = df.query('benchmark=="NPB" and implementation!="native"')
    impls = ['MKL', 'cuSPARSE', 'CL eGPU', 'SparseX', 'CL iGPU']

    df = speeds.groupby(['platform', 'name'])
    data = [
        [(row['implementation'], row['speedup'], row['platform'], row['name'])
            for _, row in group.iterrows()] 
            for _, group in df
    ]
    s_data = sorted(data, key=row_cmp)

    fig, ax = plt.subplots(figsize=fig_size(1, 0.67))

    p_colors = {
        'AMD': [0.8] * 3,
        'Intel-0': [0.5] * 3,
        'Intel-1': [0.2] * 3,
    }

    leg = {}
    rects = {}

    offset = 0
    for i, row in enumerate(s_data):
        if i == 6 or i == 12:
            offset += 0.5
        base = i + offset

        p = platform_map(row[0][2])
        r = patches.Rectangle((base-0.5, 0), 1, 16, color=p_colors[p], alpha=0.2, linewidth=0)
        ax.add_patch(r)
        rects[p] = r
        for name, speed, _, _ in row:
            nn = new_impl_map(row[0][2], name)
            leg[nn] = ax.scatter(base, speed, c=[color(nn)], marker=marker(nn))

    ax.set_title('Speedup Distribution', fontsize=10)
    ax.set_ylabel('Speedup (×)', fontsize=8)
    ax.set_xticks([])
    ax.set_ylim([0, 16])
    ax.axhline(y=1, lw=0.8, c='black')
    ax.tick_params(axis='y', labelsize=6)

    vs = [*rects.values(), *leg.values()]
    ks = [*rects.keys(), *leg.keys()]
    ax.legend(vs, ks, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize=6)

    fig.tight_layout()
    sns.despine(fig)
    return fig

plot_choices = { p.__name__ : p for p in [
    baseline,
    baseline_bench,
    expert,
    marshall,
    distribution
]}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Baseline comparison plots for ICS')
    parser.add_argument('plot', choices=plot_choices)
    parser.add_argument('data', type=str, help='Data file to tidy')
    args = parser.parse_args()

    plot_f = plot_choices[args.plot]
    fig = plot_f(get_data(args.data))
    fig.savefig('{}.pdf'.format(args.plot), bbox_inches='tight')
