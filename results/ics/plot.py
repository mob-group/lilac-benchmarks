#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
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
        'bfs': 'BFS'
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

def baseline(df):
    benches = ['pfold', 'ngt', 'PageRank', 'bfs']
    platforms = df.platform.unique()

    fig, axes = plt.subplots(1, len(benches), sharey=True, figsize=fig_size(2.1, 0.67))

    for plot_num, (bench, ax) in enumerate(zip(benches, axes)):
        i = 0
        
        if plot_num == 0:
            ax.set_ylabel('Speedup (×)')
        
        results = df.query("benchmark=='{}'".format(bench))
        
        ticks = []
        labs = []
        bars = []
        legend = []
        
        rows = [r for r in results.iterrows()]
        rows.sort(key=lambda r: platform_map(r[1].platform))

        for row in rows:
            y_val = row[1].speedup
            bars.append(ax.bar(i, y_val, width=0.95))
            ticks.append(i-0.2)
            legend.append(platform_map(row[1].platform))
            labs.append(impl_map(row[1].platform, row[1].implementation))
            
            valign = 'top' if y_val > 2.5 else 'bottom'
            offset = -0.1 if y_val > 2.5 else 0.1
            color = 'white' if y_val > 2.5 else 'black'
            ax.text(i, y_val + offset, impl_map(row[1].platform,
                row[1].implementation), rotation=90, verticalalignment=valign,
                horizontalalignment='center', color=color, fontsize=8)
            i += 1
                    
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.xaxis.set_tick_params(rotation=-45)

        ax.set_ylim(0.8, 3.0)
        ax.set_xlim(-0.6, 2.6)
        
        ax.axhline(y=1, color='white', linestyle=':')
        
        ax.set_title(bench_map(bench))
        
    ax.legend(bars, legend, bbox_to_anchor=(1.5,0.5), loc='center')
      
    fig.tight_layout()
    plt.subplots_adjust(wspace=0.2)

    sns.despine(fig)
    return fig

def marshall(df):
    # TODO: fake plot for now, but adapt
    fig, axes = plt.subplots(1, 4, figsize=fig_size(2.1, 0.67))
    # end fake plotting

    fig.tight_layout()
    sns.despine(fig)
    return fig

def expert(df):
    comparisons = [
    ]
    fig, axes = plt.subplots(2, 1, figsize=fig_size(0.9, 1.8))

    # TODO fake plot for now, but adapt
    npb = axes[0]
    par = axes[1]

    npb.set_title('NPB')
    par.set_title('Parboil')

    for b in [npb, par]:
        b.set_ylabel('LiLAC Performance (×)')
        b.set_yticks([0, 0.5, 1])
    # end fake plot

    fig.tight_layout()
    sns.despine(fig)
    return fig

plot_choices = { p.__name__ : p for p in [
    baseline,
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
