# Copyright (C) 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from collections import namedtuple
from itertools import chain
from networkx import connected_component_subgraphs
from order_util import check_nonincreasing_envelope
from pqueue import PriorityQueue
from py3compat import cPickle_loads, cPickle_dumps, cPickle_HIGHEST_PROTOCOL, izip
from test_utils import solve_test_matrices #, solve_difficult_for_ilp


def log(*args): pass
#log = print


def main():
    solve_test_matrices(solve_problem, log)
    #solve_difficult_for_ilp(solve_problem, log)
    # TODO - Implement the custom exclusion rule, then cross-check!
    #      - We need 2 heaps: one for the candidates and another one for the 
    #        remaining equations


Node = namedtuple('Node', 'cost  elims  g  heap connected')


def solve_problem(g_orig, eqs):
    # Remove empty rows and columns columns
    empty_rows = [n for n in sorted(eqs) if not g_orig[n]] 
    empty_cols = [n for n in g_orig if n not in eqs and not g_orig[n]]
    if empty_rows or empty_cols:
        g = fast_copy(g_orig)
        g.remove_nodes_from(empty_rows)
        g.remove_nodes_from(empty_cols)
        eqs_in_g = eqs - set(empty_rows)
    else:
        g = g_orig
        eqs_in_g = eqs
    #
    counter, ub, rowp = _solve_problem(g, eqs_in_g)
    #
    log('Nodes:', counter)
    log('Minimum cost:', ub)
    #
    if empty_rows:
        empty_rows.extend(rowp)
        rowp = empty_rows
    #
    rowp, colp, matches, tear_set, sink_set = get_hessenberg_form(g_orig, eqs, rowp, ub)
    record((counter, g_orig, eqs, rowp, colp, ub))
    return rowp, colp, matches, tear_set, sink_set


def _solve_problem(g_orig, eqs):
    assert isinstance(eqs, (set, dict)),'Make sure that `n in eqs` will be O(1)'
    ub, best_rowp = initial_solution(g_orig, eqs)
    #log('UB at root:', ub)
    if ub == 0:
        return 0, ub, best_rowp
    # TODO Remove zero cost rows AND columns up front
    stack = init_dfs_stack(g_orig, eqs)
    counter = 0
    while stack:
        #log('\nLevel', len(stack))
        node = stack[-1]
        if node.connected and not node.heap:
            #log('Dropping node (heap empty)')
            stack.pop()
            continue
        #
        if not node.connected:
            stack.pop()
            subcounter, subcost, rowp, cost_profiles = solve_pieces(node.g, eqs)
            # solve the pieces recursively, then increment counter by the sum,
            # set total_cost, elims appropriately 
            counter += subcounter
            total_cost = node.cost + subcost
            elims = node.elims
            elims.extend(rowp)
            #
            apply_exclusion_rule_on_pieces(stack, node.cost, cost_profiles)
            #
            if total_cost < ub:
                log('Improved UB: {} -> {}'.format(ub, total_cost))
                assert set(elims) == eqs # Only true if no preprocessing
                ub = total_cost
                best_rowp = elims
            continue
        #
        nzeros, eq = node.heap.popitem() # Why do pop- and peekitem differ in order?
        # TODO assert nzeros >= 2 # so the max below is unnecessary  
        guesses = max(0, nzeros-1)
        running_cost = node.cost + guesses
        #log('LB at node:', running_cost)
        if running_cost >= ub:
            #log('Dropping node (lb >= ub)')
            stack.pop()
            continue
        counter += 1
        elims, g = node.elims[:], fast_copy(node.g)
        #log(elims)
        rowp, heap = eliminate(g, eqs, eq)
        elims.extend(rowp)
        #log(elims)
        # TODO remove eq and all other zero cost eliminations from each 
        # parent heap if t_{r+1} there is not better than t_{s+1}
        apply_exclusion_rule(stack, running_cost, rowp)
        #
        if heap:
            stack.append(Node(running_cost, elims, g, heap, is_connected(g)))
            #log('New node added with cost', running_cost)
        elif running_cost < ub:
            log('Improved UB: {} -> {}'.format(ub, running_cost))
            assert set(elims) == eqs # Only true if no preprocessing
            ub = running_cost
            best_rowp = elims
    # TODO Append the removed zero cost rows and cols 
    #log('Nodes:', counter)
    #log('Minimum cost:', ub)
    return counter, ub, best_rowp


#@profile
def initial_solution(g_orig, eqs):
    # Duplication of heap_md.min_degree with none forbidden
    g = fast_copy(g_orig)
    eq_nzeros = PriorityQueue()
    eq_nzeros.populate((eq, len(g[eq])) for eq in sorted(eqs))
    rowp, n_tears = [ ], len(g) - len(eqs)
    while eq_nzeros:
        nzeros, eq = eq_nzeros.popitem()
        if nzeros:
            n_tears -= 1
        rowp.append(eq)
        vrs = sorted(g[eq])
        eqs_update = set(chain.from_iterable(g[v] for v in vrs))
        eqs_update.discard(eq)
        g.remove_node(eq)
        g.remove_nodes_from(vrs)
        for e in sorted(eqs_update):
            eq_nzeros[e]  = len(g[e])
    assert len(rowp) == len(eqs)
    return n_tears, rowp


def init_dfs_stack(g, eqs):
    # heap key: eq, value: nzeros
    heap = None
    connected = is_connected(g)
    if connected:
        heap = PriorityQueue()
        heap.populate((eq, len(g[eq])) for eq in sorted(eqs))
    root_node = Node(0, [ ], g, heap, connected)
    return [root_node]


def fast_copy(obj):
    pkl_str = cPickle_dumps(obj, cPickle_HIGHEST_PROTOCOL)
    return cPickle_loads(pkl_str)

#@profile
def eliminate(g, eqs, eq):
    rowp = [ eq ]
    remaining_eqs = sorted(n for n in g if n in eqs and n != eq)
    heap = PriorityQueue()
    heap.populate((r, len(g[r])) for r in remaining_eqs)
    elimination_step(g, heap, eq)
    while heap:
        eq, nzeros = heap.peekitem()
        if nzeros > 1:
            break
        heap.popitem()
        rowp.append(eq)
        elimination_step(g, heap, eq)
    return rowp, heap

#@profile
def elimination_step(g, heap, eq):
    vrs = sorted(g[eq])
    eqs_update = set(chain.from_iterable(g[v] for v in vrs))
    eqs_update.discard(eq)
    g.remove_node(eq)
    g.remove_nodes_from(vrs)
    for e in sorted(eqs_update):        # faster than the batch update
        heap[e]  = len(g[e])
    #heap.batch_update((e, len(g[e])) for e in eqs_update)


def apply_exclusion_rule(stack, running_cost, rowp):
    for node in stack:
        candids = node.heap
        for r in rowp:
            try:
                nzeros = candids[r]
            except KeyError:
                continue
            cost = node.cost + max(0, nzeros-1)
            if cost >= running_cost:
                del candids[r]


def apply_exclusion_rule_on_pieces(stack, start_cost, cost_profiles):
    profiles = [[(eq,cost+start_cost) for eq, cost in p] for p in cost_profiles]
    for node in stack:
        candids, cost = node.heap, node.cost
        for profile in profiles:
            for eq, cost_eliminated in profile:
                try:
                    nzeros = candids[eq]
                except KeyError:
                    continue
                if cost + max(0, nzeros-1) >= cost_eliminated:
                    del candids[eq]


def is_connected(g):
    assert g
    seen = set()
    nextlevel = { next(g.nodes_iter()) }
    while nextlevel:
        thislevel = nextlevel
        nextlevel = set()
        for v in thislevel:
            if v not in seen:
                seen.add(v)
                nextlevel.update(g[v])
    return len(seen) == len(g)


def argsort(lst):
    # http://stackoverflow.com/questions/3382352
    return sorted(range(len(lst)), key=lst.__getitem__)


def solve_pieces(g, eqs):
    graphs = list(connected_component_subgraphs(g))
    subeqs = [{n for n in subg if n in eqs} for subg in graphs]
    #print(subeqs)
    assert len(graphs) > 1
    sols = list(_solve_problem(*subp) for subp in izip(graphs, subeqs))
    # counter, cost, rowp = sols[i]
    total_counter = sum(t[0] for t in sols)
    total_cost    = sum(t[1] for t in sols)
    c_minus_r = [len(subg)-2*len(seqs) for subg, seqs in izip(graphs, subeqs)]
    keys = [(c_m_r, cost, min(seqs)) for c_m_r, (_, cost, _), seqs in izip(c_minus_r, sols, subeqs)]
    elims = list(chain.from_iterable(sols[i][2] for i in argsort(keys)))
    return total_counter, total_cost, elims, get_t_profiles(graphs, sols)


def get_t_profiles(graphs, sols):
    # counter, cost, rowp = sols[i]
    list_g_rowp = list(izip(graphs, (s[2] for s in sols)))
    t_profiles = [t_profile(subg, rowp) for subg, rowp in list_g_rowp]
    assert len(graphs) == len(t_profiles)
    eq_costs = [profile[-1] for profile in t_profiles]
    assert all(s[1] == eq_cost[1] for s, eq_cost in izip(sols, eq_costs))
    return t_profiles


def t_profile(g, rowp):
    profile, seen_cols, cost = [ ], set(), 0
    for eq in rowp:
        new_cols = set(g[eq]) - seen_cols
        cost += max(0, len(new_cols)-1)
        seen_cols.update(new_cols)
        profile.append((eq, cost))
    return profile


def get_hessenberg_form(g, eqs, rowp, ub):
    assert set(rowp) == eqs, (rowp, eqs)
    colp, seen_cols, matches = [ ], set(), { }
    adj = g.adj
    for eq in rowp:
        vrs = sorted(set(adj[eq]) - seen_cols)
        if vrs:
            colp.extend(vrs)
            seen_cols.update(vrs)
            var = vrs[0] # or [-1] for last
            assert eq  not in matches
            assert var not in matches
            matches[eq]  = var
            matches[var] = eq
    isolated_cols = sorted(n for n in g if n not in eqs and len(g[n])==0)
    colp.extend(isolated_cols)
    assert len(colp) == len(g) - len(eqs)
    check_nonincreasing_envelope(g, rowp, colp)
    sink_set = { n for n in rowp if n not in matches }
    tear_set = { n for n in g if n not in eqs and n not in matches }
    #_plot_bipartite(g, set(), rowp, colp)
    assert len(tear_set) == ub, (len(tear_set), ub)
    return rowp, colp, matches, tear_set, sink_set

#-------------------------------------------------------------------------------
# FIXME Clean up and remove call to record in the search algorithm

_worst_cases = [ ] # Store (counter, g_orig, eqs, best_rowp), max K_WORST length
K_WORST = 10       # and ordered in descending order by counter

def record(counter_data):
    if len(_worst_cases) < K_WORST:
        _worst_cases.append(counter_data)
        _worst_cases.sort(key=lambda tup: tup[0], reverse=True)
    elif _worst_cases[-1][0] < counter_data[0]:
        _worst_cases[-1] = counter_data
        _worst_cases.sort(key=lambda tup: tup[0], reverse=True)

def __record(counter_data):
    if counter_data[0] == 1:
        _worst_cases.append(counter_data)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
