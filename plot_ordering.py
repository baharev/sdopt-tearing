# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from sys import stderr
from codegen import gen_column_permutation
from py3compat import irange
from utils import has_matplotlib, DATADIR


def gen_r_c_color(rows, cols, rowp, colp, colors):
    to_str = {1: 'black', 2: 'red', 3: 'gray'}
    for k in irange(0, len(colors)):
        yield rowp[rows[k]], colp[cols[k]], to_str[colors[k]]

def plot_dm(name, rows, cols, rowp, colp, colors, sccs, show=True, msg=''):
    from matplotlib import pyplot as plt
    _, ax = setup(plt, len(rowp), len(colp))
    for r, c, color in gen_r_c_color(rows, cols, rowp, colp, colors):
        ax.add_artist( square(r, c, facecolor=color) )
    # Mark the border of SCC blocks yellow
    for r, c, size in sccs:
        ax.add_artist( square(r, c, size=size, facecolor='none', 
                              edgecolor='yellow', linewidth=3.0 ) )
    beautify_axes(ax)
    if show:
        ax.set_title(name)
        plt.tight_layout(pad=1.00)
        plt.show()
    else:
        plt.tight_layout(pad=1.00)
        plt.title(msg)
        # TODO Hard-coded path
        plt.savefig('/tmp/pics/'+name+'.pdf', bbox_inches='tight', pad_inches=0.05)
        plt.close()

#-------------------------------------------------------------------------------

def plot_hessenberg(g, rowp, colp, partitions, msg, mark_red=[ ]):
    # Compare with plot_sparsity
    from matplotlib import pyplot as plt
    _, ax = setup(plt, len(rowp), len(colp))
    for r, c in gen_entries(g, rowp, colp):
        ax.add_artist( square(r, c) )
    for r, c in partitions:
        draw_partition(ax, r, c)
    for r, c in mark_red:
        ax.add_artist( square(r, c, facecolor='r') )        
    beautify_axes(ax)
    ax.set_title(msg)
    plt.show()

def gen_entries(g, rowp, colp):
    indexof = { name : i for i, name in enumerate(colp) }
    rows = [ [indexof[c] for c in g[r]] for r in rowp ]
    return ((r,c) for r, row in enumerate(rows) for c in row)

def square(r, c, size=1, facecolor='k', edgecolor='0.7', linewidth=None):
    from matplotlib.pyplot import Rectangle
    # r and c must be swapped: row -> y axis, col -> x axis
    return Rectangle((c, r), size, size , facecolor=facecolor, 
                     edgecolor=edgecolor, linewidth=linewidth)

#-------------------------------------------------------------------------------

def plot_bipartite(g, forbidden, row_perm, col_perm):
    # The row and column identifiers are in permuted order
    indexof = { name : i for i, name in enumerate(col_perm) }
    rows = [ ]
    for r in row_perm:
        cols = g[r]
        rows.append( [ (indexof[c], (r,c) not in forbidden) for c in cols ] )
    #
    plot_sparsity(rows, len(col_perm))

#-------------------------------------------------------------------------------

def get_spiked_form_rowwise(blocks):
    col_perm = list(gen_column_permutation(blocks))
    indexof = { name : i for i, name in enumerate(col_perm) }
    rows = [ ]
    for blk in blocks:
        for eq in blk.eqs:
            row = [ (indexof[name], name in eq.elims) for name in eq.names ]
            rows.append(row)
        for y, x, _ in blk.conn_triples:
            rows.append( [ (indexof[x],True), (indexof[y],True) ] )
    return rows, col_perm

def plot_ordering(blocks):
    #from utils import serialize
    #serialize(blocks, 'blocks.pkl')
    #return
    rows, col_perm = get_spiked_form_rowwise(blocks)
    def func_draw_partitions(ax):
        draw_row_and_col_partitions(ax, blocks)
    #
    plot_sparsity(rows, len(col_perm), func_draw_partitions)

def plot_sparsity(rows, n_cols, func_draw_partitions=None):
    from matplotlib import pyplot as plt
    _, ax = setup(plt, len(rows), n_cols)
    for r, row in enumerate(rows):
        for c, allowed in row:
            # r and c must be swapped: row -> y axis, col -> x axis
            clr = get_color(r, c, allowed)
            rect = plt.Rectangle((c,r), 1,1, facecolor=clr, edgecolor='0.7')
            ax.add_artist(rect)
    if func_draw_partitions:     # used to be:
        func_draw_partitions(ax) # draw_row_and_col_partitions(ax, blocks)
    beautify_axes(ax)
    plt.show()

def get_color(r, c, allowed):
    # r and c must be swapped: row -> y axis, col -> x axis
    if c <= r: # below or on the diagonal
        clr = 'black' if allowed else 'grey'
    else:
        clr = 'red'
    return clr

def draw_row_and_col_partitions(ax, blocks):
    pos = 0
    for blk in blocks:
        if blk.eqs:
            pos += len(blk.eqs)
            draw_partition(ax, pos, pos)
        if blk.conn_triples:
            pos += len(blk.conn_triples)
            draw_partition(ax, pos, pos)

def draw_partition(ax, r, c):
    line_color, line_width = 'blue', 1
    ax.axhline(r, c=line_color, lw=line_width)        
    ax.axvline(c, c=line_color, lw=line_width)   

# Compare plotting with SDOPT
def setup(plt, nrows, ncols):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    mng = plt.get_current_fig_manager()
    mng.resize(1865,1025)
    #mng.full_screen_toggle()     
    plt.axis('scaled')
    ax.set_xlim([0, ncols])
    ax.set_ylim([0, nrows])
    return fig, ax

def beautify_axes(ax):
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])

def no_matplotlib(*args):
    stderr.write('Plotting requires a working matplotlib installation.\n')

plot_sparsity = plot_sparsity if has_matplotlib() else no_matplotlib

def main():
    from utils import deserialize
    blocks = deserialize(DATADIR + 'blocks.pkl.gz')
    plot_ordering(blocks)

if __name__ == '__main__':
    main()
