# Copyright (C) 2016 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
#
from __future__ import print_function, division
from six import iteritems

__all__ = [ 'hellerman_rarick_p_3_5' ]


def main():
    try:
        import mc33 #@UnusedImport
    except ImportError:
        import sys
        sys.stderr.write('Requires the MC33 wrapper.\n')
        return
    from testmatrices import create_difficult_pattern, create_block_pattern
    from networkx import Graph
    #
    n_eqs = 16
    g = create_difficult_pattern(n_eqs)
    solve(g, n_eqs)
    #
    g, n_eqs = create_block_pattern(12)
    solve(g, n_eqs)
    #
    g = Graph()
    n_eqs = 51
    g.add_edges_from((e, v+n_eqs) for e, v in zip(EDGELIST[::2], EDGELIST[1::2]))
    solve(g, n_eqs)


def solve(g, n_eqs):
    from plot_ordering import to_pdf
    from six.moves import xrange as range
    print('------------------------------------------------------------------')
    print('Solving problem of size', n_eqs)
    fname = '{0:03d}a'.format(n_eqs)
    to_pdf(g, list(range(n_eqs)), range(n_eqs, 2*n_eqs), fname=fname)
    #
    rowp, colp = hellerman_rarick_p_3_5(g, set(range(n_eqs)))
    #
    fname = '{0:03d}b'.format(n_eqs)
    to_pdf(g, rowp, colp, fname=fname)
    #
    print('Plot written!')


def hellerman_rarick_p_3_5(g, eqs):
    from mc33 import hr_35
    from numpy import zeros, ones, int32 as fint # Fortran INTEGER
    #
    vrs = sorted(n for n in g if n not in eqs)
    eqs = sorted(eqs)
    M = len(eqs)
    N = len(vrs)
    #
    e2i = {e: i for i, e in enumerate(eqs, 1)}
    i2e = {i: e for e, i in iteritems(e2i)}    # could be just an array
    #
    v2j = {v: j for j, v in enumerate(vrs, 1)}
    j2v = {j: v for v, j in iteritems(v2j)}    # could be just an array
    #
    NZI = g.number_of_edges() #  sum(len(adj[v]) for v in vrs) would be faster
    #
    IRN, JCN, A = zeros(NZI, fint), zeros(NZI, fint), ones(NZI, 'float64')
    #
    for k, (e, v) in enumerate(g.edges_iter(eqs)):
        IRN[k] = e2i[e]  # inefficient but leave it for now
        JCN[k] = v2j[v]
    #
    IP, IQ, IPROF, IFLAG, IW, IW1 = zeros(M, fint), zeros(N, fint), \
      zeros(N, fint), zeros(3, fint), zeros(M+N, fint), zeros(9*N+3*M, fint)
    #
    #nzo,ierr = hr_35(itype,a,irn,jcn,ip,iq,iprof,iflag,iw,iw1,m=len(ip),n=len(iq),nzi=len(a))
    NZO, IERR = hr_35(    3,A,IRN,JCN,IP,IQ,IPROF,IFLAG,IW,IW1)
    assert IERR == 0 and NZI-NZO == 0, (IERR, NZI, NZO)
    #
    rowp = [i2e[i] for i in IP]
    colp = [j2v[j] for j in IQ]
    #
    return rowp, colp


EDGELIST = [
        20, 25, 11, 11, 49, 42, 12, 12, 22, 19, 49, 45, 46, 49, 36, 34, 15,
        14, 45, 45, 43, 42, 31, 21,  5,  5, 11,  5, 36, 36, 17,  2, 47, 38,
        25, 15,  1,  1, 42, 41, 32, 28, 28, 31, 20, 20, 43, 36, 38, 38, 46,
        46, 39, 39, 30, 29, 16,  0, 31, 24, 25, 21, 16, 19, 14,  2, 13, 10,
        35, 35, 22, 25, 27, 27, 42, 42, 29, 20, 28, 28,  9,  9, 38, 35, 48,
        47, 30, 30, 35, 26, 16, 13, 47, 47, 20, 16, 16, 16, 21, 21, 26, 23,
        28, 22,  2,  1, 29, 23,  9,  4, 24, 22,  5,  1, 44, 41, 46, 43, 12,
         2, 19, 18, 12, 13, 50, 46, 22, 16, 18, 16, 36, 35,  2,  2,  8,  6,
        29, 26,  5,  4, 49, 39, 14, 19, 34, 34, 21, 20, 17,  1, 49, 48, 47,
        41, 45, 44, 27, 31, 34, 37, 32, 29, 43, 39, 29, 29, 10, 10, 34, 31,
        39, 38, 11,  3,  7,  4, 31, 27, 19, 12, 13,  9, 17, 17, 27, 26, 50,
         3,  4,  8,  3,  0,  9,  8, 32, 32, 35, 29, 15, 19, 31, 30, 20, 17,
        13, 12, 23, 17, 41, 38,  0, 50,  2, 13, 43, 33, 10,  0, 46, 40, 48,
        48, 26, 26, 25, 18, 40, 43, 18, 17, 23, 20,  8,  7, 37, 33,  6,  8,
        45, 49, 25, 24, 26, 31, 33, 37, 48, 46, 40, 40, 18, 18, 28, 25,  8,
         4, 37, 36, 34, 28,  6,  6, 38, 43, 44, 44,  7,  7, 14, 13, 41, 41,
        19, 15, 32, 37, 37, 27, 33, 32, 40, 37,  1,  2, 42, 40,  7,  5,  3,
         3,  4,  4,  6,  3, 30, 28,  8,  8, 44, 49, 24, 23, 14, 14, 19, 10,
        50, 47, 17, 14, 15, 15, 13, 11, 37, 30,  0,  0, 40, 34,  3, 13, 35,
        32, 50,  5, 41, 35, 41, 32,  1, 50, 23, 14, 38, 34, 10,  1, 39, 43,
        47, 44, 21, 25, 22, 22, 26, 22,  0, 13, 23, 23, 33, 33,  3,  8,  2,
         0, 24, 24,  5,  2, 44, 40]


if __name__ == '__main__':
    main()
