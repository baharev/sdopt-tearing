# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from cPickle import dumps, loads, HIGHEST_PROTOCOL
from itertools import count
from sys import maxsize as INT_MAX
from time import time
from benchmarks import gen_benchmarks_as_undirected_bipartite
from heap_md import min_degree as heuristic_solve # no lookahead
#from min_degree import min_degree as heuristic_solve # with lookahead
from plot_ordering import plot_bipartite
#from plot_ordering import plot_hessenberg
from pqueue import PriorityQueue
from test_dense import create_test_problem # , gen_testproblems 
from utils import pairwise


def main():
    opt = { 'block size'      : 5, 
            'spike width'     : 2,
            'number of blocks': 3,
            'border width'    : 3,
          }
    OPT = opt['spike width']*opt['number of blocks'] + opt['border width']
    g, eqs = create_test_problem(opt)
    print('Problem size: {}x{}'.format(len(eqs), len(g)-len(eqs)))
    print('Optimum:', OPT)
    start = time()
    solve_problem(g, eqs, OPT=OPT)
    end = time()
    print('Overall solution time: {0:0.1f} s'.format(end-start))


def main2():
    for g, eqs, _ in gen_benchmarks_as_undirected_bipartite():
        print('Problem size: {}x{}'.format(len(eqs), len(g)-len(eqs)))
        row_perm = sorted(eqs)
        col_perm = sorted(n for n in g if n not in eqs)
        plot_bipartite(g, set(), row_perm, col_perm)
        start = time()
        solve_problem(g, eqs)
        end = time()
        print('Overall solution time: {0:0.1f} s'.format(end-start))
        return

class Elimination:
    
    def __init__(self, g, heap, parent, rowp, colp, part):
        self.g      = g      # bipartite graph of the system of equations
        self.heap   = heap   # row index - cost min-heap
        self.parent = parent # row, column permutation key of parent 
        self.rowp   = rowp   # row permutation
        self.colp   = colp   # column permutation
        self.part   = part   # (r,c) partitions, len of rowp and colp at guesses
    
    def finish_cloning(self):
        self.g    = loads(self.g)
        self.heap = loads(self.heap)
        self.rowp = self.rowp[:]
        self.colp = self.colp[:]
        self.part = self.part[:]
    
    def mark_block_boundary(self):
        self.part.append( (len(self.rowp), len(self.colp)) )


# TODO - Verify that the found solution is a valid elimination at the given cost
#      - Careful: we may not be able to improve the heuristic solution 
#        in the `while search_heap:` loop!
#      - Examine the effect of randomizing the row order
#      - UB and cutoff are different (if cutoff comes from blocks, or user 
#        defined like gap, say 10% better than UB)
#      + Play with different ordering strategies

def solve_problem(g, eqs, LB=0, UB=INT_MAX, OPT=None, SUBPROBLEM=False):
    assert LB <= UB
    assert OPT is None or LB <= OPT
    assert OPT is None or UB >= OPT
    
    _, _, _, tear_set, _ = heuristic_solve(g, eqs)
    UB = min(UB, len(tear_set))
    # TODO Get lower bound from feasible solution
    if LB==UB:
        assert OPT is None or OPT==UB, (OPT, UB)
        return UB # TODO Return the solution
    
    counter = count()
    memo = { }
    search_heap = init_search_heap(g, eqs, UB)
    
    # TODO Trivial lower bound: the lb of the first item in the min heap
    #                           + we also have a feasible solution by now
    while search_heap:
        msgfmt = 'Search heap size: {}, LB = {}, UB = {}, {}' 
        print(msgfmt.format(len(search_heap),LB,UB,'S' if SUBPROBLEM else 'M'))
        (cost_at_start, _, first_row), elim = search_heap.popitem()        
        children_ids, _ = memo.get(elim.parent, (None, None))
        if children_ids and id(elim) not in children_ids:
            continue  # parent was superseded, and children deleted in branch()
        if cost_at_start >= UB:
            assert LB < UB, (LB, UB) # ? <=
            LB = UB
            print('Aborting search at cost:', UB)
            break
        LB = max(LB, cost_at_start)
        next(counter)
        elim.finish_cloning() 
        elim.mark_block_boundary() # TODO First iter: assumes we have to guess
        solve_equation(first_row, elim) # Compulsory elimination with guesses
        while elim.heap:  # Eliminate rows that now don't require guesses
            eq, cost = elim.heap.peekitem()
            if not cost:
                solve_equation(eq, elim)
            else:
                break
        if elim.heap:  # If not done yet, we need to branch
            # Branch also drills down at the node being branched
            feas_sol = branch(search_heap, memo, UB, elim, cost_at_start, eqs)
            if feas_sol is None: # parent key of elim already in search_heap 
                continue         # with no worse cost, so discarding elim  
            elim, cost_at_start = feas_sol
        #
        elim.mark_block_boundary()
        msg = 'cost: {}, (LB = {}, UB prev = {}, OPT = {})'
        msg = msg.format(cost_at_start, LB, UB, OPT)
        print('*** Feasible solution with', msg)
        UB = min(UB, cost_at_start)
        # TODO Save best solution
        print(elim.part)
        assert cost_at_start, (cost_at_start, UB)
        assert not elim.g, len(elim.g)
        assert len(elim.rowp)>=len(elim.colp)
        assert SUBPROBLEM or len(elim.rowp)==len(elim.colp) 
        
        lb, cslc_lb, _, slcs_to_resolve = trivial_lower_bound(g,elim,OPT)
        
        #plot_hessenberg(g, elim.rowp, elim.colp, elim.part, msg+' lb = %d'% lb)
        
        if not SUBPROBLEM and slcs_to_resolve:
            # resolve those slices where the trivial lb might be improved   
            resolve_slices(g, elim, cslc_lb, slcs_to_resolve)
            new_lb = sum(v for v in cslc_lb.itervalues())
            if new_lb > lb:
                print('Overall lower bound improvement:', lb, '->', new_lb)
                lb = new_lb
        LB = max(lb, LB)
        
        #if not SUBPROBLEM:
        #    lb_with_right_nbrs(g, elim, cslc_lb, cslc_ub)
        
        assert LB <= UB
        assert OPT is None or LB <= OPT
        assert OPT is None or UB >= OPT
        if LB == UB:
            print('Optimal solution found')
            print('Search heap size:', len(search_heap))
            break
    # Done.
    if not search_heap:
        assert LB <= UB, (LB, UB)
        LB = UB
    assert LB==UB, (LB, UB)
    assert OPT is None or UB == OPT, (UB, OPT)
    print('Optimum found:', UB)
    print('Iterations:   ', next(counter))
    return UB


def trivial_lower_bound(g, elim, OPT):
    lb, cslc_lb, cslc_ub, slcs_to_resolve = 0, { }, { }, [ ]
    for (first_row_idx,start), (_,end) in pairwise(elim.part):
        # TODO We don't need the whole row_cost map, only the min element
        row_cost = get_row_cost(g, elim.colp[start:end])
        min_cost = min(row_cost.iteritems(), key=lambda t: t[1])[1]
        # spent_cost: columns in first row that are also in this column slice -1
        spent_cost = row_cost[ elim.rowp[first_row_idx] ]
        cslc_ub[(start,end)] = spent_cost
        print('Col slice:', (start,end),' lb:', min_cost,' spent:', spent_cost)
        assert min_cost <= spent_cost, (min_cost, spent_cost)
        cslc_lb[(start,end)] = min_cost
        if min_cost < spent_cost:
            slcs_to_resolve.append((start,end))
        lb += min_cost
    assert OPT is None or lb <= OPT, (lb, OPT)
    return lb, cslc_lb, cslc_ub, slcs_to_resolve

def resolve_slices(g, elim, cslc_lb, slcs_to_resolve):
    for start, end in slcs_to_resolve:
        lb = cslc_lb[(start,end)] # trivial lb
        cols = elim.colp[start:end]
        rows = set(r for c in cols for r in g[c])
        nbunch = cols[:]
        nbunch.extend(rows)
        print('<<< Recursing into subproblem')
        # TODO Deepcopy would be necessary if we had graph attributes
        new_lb = solve_problem(g.subgraph(nbunch), rows, lb, SUBPROBLEM=True)
        print('>>> Recursion finished')
        assert new_lb >= lb, (new_lb, lb)
        if new_lb > lb:
            cslc_lb[(start,end)] = new_lb
            print('Improved lb:', lb, '->', new_lb)

#-------------------------------------------------------------------------------
# Currently dead code. See notes from Apr 08, 2015: We don't know which are the
# correct right neighbors. Although it may be useful for a local search 
# improvement.

def lb_with_right_nbrs(g, elim, cslc_lb, cslc_ub):
    
    # (start,end) slices, name -> slice map
    colslcs, colname_slc = colname_slc_map(elim)
    rowslcs, rowname_slc = rowname_slc_map(elim)
    
    # We now essentially build the sparse block sparsity pattern both row- and
    # column-wise. 
    
    # row slice -> nonempty col slices
    rowslc_colslcs = { }
    for rstart, rend in rowslcs:
        rows = elim.rowp[rstart:rend]
        cols = set(c for r in rows for c in g[r])
        cslcs = set(colname_slc[c] for c in cols)
        rowslc_colslcs[(rstart,rend)] = cslcs
    
    # col slice -> nonempty row slices
    colslc_rowslcs = { }
    for cstart, cend in colslcs:
        cols = elim.colp[cstart:cend]
        rows = set(r for c in cols for r in g[c])
        rslcs = set(rowname_slc[r] for r in rows)
        colslc_rowslcs[(cstart,cend)] = rslcs
    
    print()
    print(elim.part)
    
    for cstart, cend in colslcs:
        print()
        # Non-empty row slices in (cstart,cend) column slice
        rslcs = sorted( colslc_rowslcs[(cstart,cend)] )
        # sorted so that the logging will be nice
        print('Row slices in col slice {}: {}'.format((cstart,cend), rslcs))
        
        # Strictly right neighbor column slices of (cstart,cend) column slice
        this_slc = (cstart,cend)
        nbrs = set()
        for rslc in rslcs:
            nbrs.update(slc for slc in rowslc_colslcs[rslc] if slc > this_slc)
        print('Strictly right nbrs of column slice {}: {}'.format(this_slc, \
                                                                  sorted(nbrs)))
        # Trivial lb on row slices in (cstart,cend), including right neighbors
        col_set = set(elim.colp[cstart:cend])
        for slc in nbrs:
            col_set.update(elim.colp[slice(*slc)])
        rslc_mincost = { }
        for rslc in rslcs:
            rslc_mincost[rslc] = min_rowcost(g, elim.rowp[slice(*rslc)],col_set)
        print('Col slice ub:', cslc_ub[this_slc])
        print('Row slice : lb ', sorted(rslc_mincost.iteritems()))
        
        # Recomputing lb per row slice
        rslc_lb = { }
        for rslc in rslcs:
            all_cnbrs = rowslc_colslcs[rslc]
            rhs_cnbrs = sorted(slc for slc in all_cnbrs if slc > this_slc)
            lb = rslc_mincost[rslc] - sum(cslc_ub[slc] for slc in rhs_cnbrs)
            print('r: {} c: {}  lb: {}'.format(rslc, this_slc, lb))
            rslc_lb[rslc] = lb
        # new lb on this column slice
        lb = min(rslc_lb.iteritems(), key=lambda t: t[1])[1]
        old_lb = cslc_lb[this_slc] 
        if lb > old_lb:
            print('*** Improved {}: {} -> {}'.format(this_slc, old_lb, lb))
            cslc_lb[this_slc] = lb
    
    new_lb = sum(v for v in cslc_lb.itervalues())
    print('new lb:', new_lb)
    print()
    return

def min_rowcost(g, row_set, col_set):
    # cost of r: all rows in r that are also in col_set -1
    return min( len(set(g[r]) & col_set)-1 for r in row_set )

def colname_slc_map(elim):
    'Returns: column (start,end) slices, column name -> column slice map.'
    colslices = tuple((start,end) for (_,start),(_,end) in pairwise(elim.part))
    return colslices, _name_to_slc_map(colslices, elim.colp)

def rowname_slc_map(elim):
    'Returns: row (start,end) slices, row name -> row slice map.'
    rowslices = tuple((start,end) for (start,_),(end,_) in pairwise(elim.part))
    return rowslices, _name_to_slc_map(rowslices, elim.rowp)

def _name_to_slc_map(slices, perm):
    # index_slc: row/col index in permutation -> row/col slice (start, end)
    index_slc = [ None ]*len(perm)
    counter = count()
    for start, end in slices:
        index_slc[start:end] = [ slices[next(counter)] ] * (end-start)
    # name_slc: row/col name -> row/col slice (start, end)
    name_slc = { name : index_slc[i] for i, name in enumerate(perm) }  
    return name_slc  

#-------------------------------------------------------------------------------

def get_row_cost(g, cols):
    row_cost = { }
    col_set = set(cols)
    rows = ( r for c in cols for r in g[c] if r not in row_cost )
    for r in rows:
        nbr_cols = set(g[r])
        row_cost[r] = len(nbr_cols & col_set)-1
    return row_cost

def drill_down(g_pkl, heap_pkl, rowp, colp, part, running_cost):
    # test just one of the min cost new elims
    elim = Elimination(g_pkl, heap_pkl, None, rowp, colp, part)
    elim.finish_cloning()
    while elim.heap:
        (eq, cost) = elim.heap.peekitem()
        if cost:
            elim.mark_block_boundary()
            running_cost += len(elim.g[eq])-1
        solve_equation(eq, elim)
    return elim, running_cost

def branch(search_heap, memo, UB, elim, cost_at_start, eqs):
    # Have we met this parent? Skip if yes; discard others if elim is superior.
    rowp, colp = elim.rowp, elim.colp
    key_elim = (tuple(sorted(rowp)), tuple(sorted(colp)))
    _, prev_cost = memo.get(key_elim, (None,None))
    if prev_cost is not None and prev_cost <= cost_at_start:
        return  # Already in the search heap, with no worse cost
    # TODO Consider moving all the above code up one level
    # Either not seen yet, or seen with worse cost; in the latter case discard
    ids = set() # all children that stem from the parent with the worse key
    memo[key_elim] = [ids, cost_at_start]
    # Copy as much as possible outside the loop
    g, heap, part = elim.g, elim.heap, elim.part
    g_pkl    = dumps(g,    HIGHEST_PROTOCOL)
    heap_pkl = dumps(heap, HIGHEST_PROTOCOL)
    equations = sorted( n for n in g if n in eqs )
    remaining_rows = len(equations)
    for eq in equations:
        new_cost_at_start = cost_at_start + len(g[eq])-1
        assert new_cost_at_start
        if new_cost_at_start < UB:
            new_elim = Elimination(g_pkl, heap_pkl, key_elim, rowp, colp, part)
            ids.add( id(new_elim) )
            priority = (new_cost_at_start, remaining_rows, eq)
            # TODO The heap could be updated in batch instead
            search_heap[new_elim] = priority
    # Dive at the node where we branch, return feasible solution with cost
    # TODO We could skip drill_down, if we have met this parent, however, be 
    #      careful NOT to discard new elims in the `while search_heap:` loop!
    return drill_down(g_pkl, heap_pkl, rowp, colp, part, cost_at_start)

def init_search_heap(g, eqs, UB):
    row_cost, heap = [ ], PriorityQueue()
    for eq in sorted(eqs):
        cost = len(g[eq])-1
        assert cost >= 0
        if cost < UB:
            row_cost.append((eq,cost))
            heap[eq] = cost
    assert row_cost
    search_heap = PriorityQueue()
    remaining_rows = len(heap)
    g_pkl    = dumps(g,    HIGHEST_PROTOCOL)
    heap_pkl = dumps(heap, HIGHEST_PROTOCOL)
    parent_key = None
    for first_row, cost_at_start in row_cost:
        elim_heap = Elimination(g_pkl, heap_pkl, parent_key, [ ], [ ], [ ])
        priority = (cost_at_start, remaining_rows, first_row)
        search_heap[elim_heap] = priority
    return search_heap

def solve_equation(row_index, elim):
    g = elim.g
    vars_in_row = sorted(g[row_index])
    elim.rowp.append(row_index)
    elim.colp.extend(vars_in_row)    
    nbr_rows = connected_rows(g, row_index)
    g.remove_node(row_index)
    g.remove_nodes_from(vars_in_row)
    del elim.heap[row_index]
    for eq in nbr_rows:
        ##--- TODO Remove when done with debugging
        #old_cost = elim.heap[eq]
        #new_cost = max(0, len(g[eq])-1)
        #assert new_cost < old_cost or (new_cost==0 and old_cost==0)
        ##---
        elim.heap[eq] = max(0, len(g[eq])-1) 

def connected_rows(g, eq):
    reachable_eq = set()
    for v in g[eq]:
        reachable_eq.update(g[v])    
    reachable_eq.discard(eq) # if eq has no connected vars, it is not in the set
    return sorted(reachable_eq)

if __name__ == '__main__':
    main()
