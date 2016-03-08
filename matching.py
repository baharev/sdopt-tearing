# Copyright (C) 2016 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
#
from __future__ import print_function, division

__all__ = ['maxmatch_len', 'nx_maxmatch_len', 'max_matching']

def nx_maxmatch_len(g, _eqs=None):
    from networkx import max_weight_matching
    return len(max_weight_matching(g, maxcardinality=True))//2


def maxmatch_len(g, eqs):
    _iperm, match_len = max_matching(g, eqs)
    return match_len

def has_hsl_and_numpy():
    try:
        import mc21  #@UnusedImport
        import numpy #@UnusedImport
        return True
    except ImportError:
        return False

maxmatch_len = maxmatch_len if has_hsl_and_numpy() else nx_maxmatch_len

# FIXME: Should do i2a to get back the original permutation!
#        Now IPERM is also in base 1 (Fortran!). 
#
# A code snippet, dumping the diagonalized matrix   
#    from matching import max_matching
#    iperm, match_len = max_matching(g, set(range(n_eqs)))
#    assert match_len == 51 
#    rowp = list(e-1 for e in iperm) # Fortran to base 0 indices
#    colp = list(range(n_eqs, 2*n_eqs))
#    
#    fname = '{0:03d}b'.format(n_eqs)
#    to_pdf(g, rowp, colp, fname=fname)   
    
def max_matching(g, eqs_orig):
    # See the FIXME above!
    # Precondition: eqs_orig defines a bipartite node set 
    from mc21 import maxmatch  # Compile as f2py -c -m mc21 mc21.f
    from numpy import cumsum, fromiter, zeros, int32 as fint # Fortran INTEGER
    #
    adj = g.adj
    eqs = sorted(e for e in eqs_orig if adj[e])  # non-isolated equations
    vrs = sorted(n for n in adj if n not in eqs_orig and adj[n])  # non-isolated vars
    #
    if not eqs or not vrs:
        return 0
    #
    A, B = (eqs, vrs) if len(eqs) >= len(vrs) else (vrs, eqs)
    assert len(A) >= len(B)
    #
    # Copy the bipartite graph into arrays as required by MC21, see the doc
    N = len(A)
    #i2a = { enumerate(A, 1) }
    b2j = { b: j for j, b in enumerate(B, 1) }
    cols_rowwise = fromiter((b2j[b] for a in A for b in adj[a]), fint)
    tmp_lenr = zeros(N+1, fint) # prepend the leading zero
    tmp_lenr[1:] = fromiter((len(adj[a]) for a in A), fint)
    col_start = cumsum(tmp_lenr[:-1]) + 1
    lenr = tmp_lenr[1:]
    #
    IPERM, NUMNZ, IW = zeros(N, fint), zeros(1, fint), zeros(4*N, fint) 
    # import mc21; help(mc21) in the REPL gives: 
    # maxmatch(icn,ip,lenr,iperm,numnz,iw,n=len(ip),licn=len(icn),liw=len(iw))
    assert len(cols_rowwise[(cols_rowwise < 1) | (cols_rowwise > N)]) == 0
    maxmatch(cols_rowwise, col_start, lenr, IPERM, NUMNZ, IW)
    #
    match_len = NUMNZ[0] 
    max_possible =  min(len(A), len(B))
    assert match_len > 0, match_len
    assert match_len <= max_possible, (match_len, max_possible)
    return IPERM, match_len
