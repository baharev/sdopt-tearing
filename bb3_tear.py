# Copyright (C) 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from itertools import chain
from time import time
from networkx import connected_component_subgraphs
from namedlist import namedlist
from order_util import check_nonincreasing_envelope
from plot_ordering import to_pdf
from pqueue import PriorityQueue
from py3compat import cPickle_loads, cPickle_dumps, cPickle_HIGHEST_PROTOCOL, \
                      irange, izip
from testmatrices import create_difficult_pattern, create_block_pattern
from utils import print_timestamp
#from test_utils import solve_test_matrices#, solve_difficult_for_ilp

TIME_LIMIT = 3600

class TimeLimitReached(Exception):
    pass

def log(*args, **kwargs):
    pass
  
def log_indent(*args, **kwargs):
    pass

# log = print
#  
# def log_indent(*args, **kwargs):
#     level = kwargs.pop('level')
#     log('    '*level, end='')
#     log(*args, **kwargs)

#-------------------------------------------------------------------------------

def main():
    print_timestamp()
    n_eqs = 16
    g = create_difficult_pattern(n_eqs)
    solve(g, n_eqs, solve_problem)
    #
    g, n_eqs = create_block_pattern(12)
    solve(g, n_eqs, solve_problem)
    #
    #solve_test_matrices(solve_problem, log)
    #solve_difficult_for_ilp(solve_problem, log)

def solve(g, n_eqs, func):
    print('------------------------------------------------------------------')
    print('Solving problem of size', n_eqs)
    msg   = 'Size: ' + str(n_eqs)
    fname = '{0:03d}a'.format(n_eqs)
    to_pdf(g, list(irange(n_eqs)), irange(n_eqs, 2*n_eqs), msg, fname)
    #
    res = func(g, set(irange(n_eqs)))
    #
    print('Explored', res.explored, 'nodes')
    msg   = 'OPT = {}, BT: {}'.format(res.ub, res.explored)
    fname = '{0:03d}b'.format(n_eqs)
    to_pdf(g, res.rowp, res.colp, msg, fname)
    print_timestamp()

#-------------------------------------------------------------------------------

# counter: to generate unique node IDs; history: map of node id -> [ events ];
# asm: key -> node id where the asm was seen 

History = namedlist('History', 'counter  history  asm  incumbent  start')

_H = None

def reset_history():
    global _H
    _H = History(counter=0, history={}, asm={}, incumbent=None, start=time())

def register(who, recursion_level):
    assert _H.counter not in _H.history, 'Forgot to call reset_history?'
    _H.history[_H.counter] = [_H.counter, who, 'rlevel %d' % recursion_level]
    _H.counter += 1

def get_id():
    return _H.counter

def explored_so_far():
    return _H.counter

def add_event(node_id, event):
    _H.history[node_id].append(event)

def history(node_id_or_key):
    if isinstance(node_id_or_key, int):
        return _H.history[node_id_or_key]
    else:
        assert isinstance(node_id_or_key, tuple), type(node_id_or_key)
        node_id = _H.asm[node_id_or_key]
        return _H.history[node_id]

def store_asm(key, node_id):
    assert key not in _H.asm
    _H.asm[key] = node_id

def update_asm(key, node_id):
    assert key in _H.asm
    _H.asm[key] = node_id

def set_incumbent(indent, ub, best_rowp):
    if not indent:
        _H.incumbent = (ub, best_rowp)

def check_time_limit():
    elapsed = time() - _H.start
    if elapsed > TIME_LIMIT:
        raise TimeLimitReached()        

#-------------------------------------------------------------------------------

TearingResult = namedlist('TearingResult', '''rowp  colp  matches  tear_set  
                           sink_set  ub  explored  optimal  time''')

def solve_problem(g_orig, eqs, log=log):
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
    reset_history()
    seen_asm = { }
    start = time()
    try:
        ub, rowp = _solve_problem(g, eqs_in_g, seen_asm)
        optimal = True
    except TimeLimitReached:
        assert _H.incumbent is not None
        log('Time limit of {}s reached!'.format(TIME_LIMIT))
        ub, rowp = _H.incumbent 
        optimal = False
    #
    explored = explored_so_far()
    log('Nodes explored:', explored)
    log('Minimum cost:  ', ub)
    #
    if empty_rows:
        empty_rows.extend(rowp)
        rowp = empty_rows
    #
    ub += len(empty_cols)
    t = get_hessenberg_form(g_orig, eqs, rowp, ub)
    elapsed = time() - start
    return TearingResult(*t, ub=ub, explored=explored, optimal=optimal, time=elapsed)


Node = namedlist('Node', '''cost  elims  g  heap  connected  remaining_eqs  
                            rowp_to_parent  asm_key  fin_cost  fin_seq  dbg_id''')

def _solve_problem(g_orig, eqs, seen_asm, indent=0):
    assert isinstance(eqs, (set, dict)),'Make sure that `n in eqs` will be O(1)'
    ub, best_rowp = initial_solution(g_orig, eqs)
    log_indent('UB at root:', ub, level=indent)
    #ub = 1000000
    if ub == 0:
        return ub, best_rowp
    set_incumbent(indent, ub, best_rowp)
    stack = init_dfs_stack(g_orig, eqs, indent)
    while stack:
        log()
        log_indent('Level ', len(stack), level=indent)
        node = stack[-1]
        if node.connected and not node.heap:
            log_indent('Dropping node (heap empty)', level=indent)
            stack.pop()
            propagate_finishing_sequence(stack, node, ub, seen_asm)
            continue
        if not node.connected:
            stack.pop()
            ub, best_rowp, guesses, rowp = \
                _solve_disconnected(stack, node, eqs, seen_asm, indent, ub, best_rowp)
            # FIXME Chunk up the rowp at least if the permutations are skipped
            store_finishing_sequence(node, guesses, rowp)
            propagate_finishing_sequence(stack, node, ub, seen_asm)
            continue
        if node.asm_key in seen_asm:
            log_indent('Seen ASM', level=indent)
            fin_lb = seen_asm[node.asm_key][0]
            cost_lb = node.cost + fin_lb
            if cost_lb >= ub:
                log_indent('ASM discarding: {} >= {}'.format(cost_lb, ub), level=indent)
                stack.pop()
                store_lower_bound(node, fin_lb)
                propagate_finishing_sequence(stack, node, ub, seen_asm)
                continue
            log_indent('Solving ASM: {} < {}'.format(cost_lb, ub), level=indent)     
        running_cost, eq, guesses = _next_eq(node, indent)
        if running_cost >= ub:
            log_indent('Dropping node (lb >= ub)', level=indent)
            stack.pop()
            store_lower_bound(node, guesses)
            propagate_finishing_sequence(stack, node, ub, seen_asm)
            continue
        elims, g, heap, remaining_eqs, rowp = _solve_eq(stack, node, indent, 
                                                        running_cost, eq)
        if heap:
            asm_key = to_asm_key(remaining_eqs)
            stack.append(Node(running_cost, elims, g, heap, is_connected(g), 
                              remaining_eqs, rowp, asm_key, None, None, get_id()))
            register('dfs down', indent)
            log_indent('New node added with cost', running_cost, level=indent)
        else: # running_cost < ub guaranteed by the above check!
            assert running_cost < ub, (running_cost, ub)
            log_indent('Improved UB: {} -> {}'.format(ub, running_cost), level=indent)
            assert set(elims) == eqs
            ub = running_cost
            best_rowp = elims
            set_incumbent(indent, ub, best_rowp)
            store_finishing_sequence(node, guesses, rowp)
        check_time_limit()
    log_indent('Nodes:', explored_so_far(), level=indent)
    log_indent('Minimum cost:', ub, level=indent)
    return ub, best_rowp

def to_asm_key(iterable):
    return tuple(sorted(iterable))

def store_lower_bound(node, guesses):
    if node.fin_cost is None or guesses < node.fin_cost:
        node.fin_cost = guesses
        node.fin_seq  = None

def store_finishing_sequence(node, guesses, rowp):
    if node.fin_cost is None or guesses < node.fin_cost or \
      (guesses == node.fin_cost and node.fin_seq is None):
        node.fin_cost = guesses
        node.fin_seq  = [ rowp ]

def update_parent(node, parent):
    fin_cost = (node.cost - parent.cost) + node.fin_cost
    if parent.fin_cost is None or fin_cost < parent.fin_cost or \
      (fin_cost == parent.fin_cost and parent.fin_seq is None):
        parent.fin_cost = fin_cost
        if node.fin_seq is None:
            parent.fin_seq = None
        else:
            parent.fin_seq = [ node.rowp_to_parent ]
            parent.fin_seq.extend(node.fin_seq)

# FIXME Simply take the node as argument below:
def _upsert(seen_asm, key, node_fin_cost, node_fin_seq, node_id, node_rem_eqs):
    assert node_fin_cost is not None, history(node_id)
    assert node_fin_seq is None or \
      set(chain.from_iterable(node_fin_seq)) == node_rem_eqs, history(node_id)
    # seen_asm value: (cost, fin_seq) 
    log('key =', key)
    if key not in seen_asm:
        seen_asm[key] = (node_fin_cost, node_fin_seq)
        store_asm(key, node_id)
        return
    #---------------------------------------------------
    # key in seen_asm:
    cost, fin_seq = seen_asm[key]
    # stored lb - node lb  
    if fin_seq is None and node_fin_seq is None:
        if cost < node_fin_cost: # update only if the new lower bound is better
            seen_asm[key] = (node_fin_cost, node_fin_seq)
            update_asm(key, node_id)
    # stored lb - node finishing seq
    elif fin_seq is None and node_fin_seq is not None:
        assert cost <= node_fin_cost, (cost, node_fin_cost)
        seen_asm[key] = (node_fin_cost, node_fin_seq)
        update_asm(key, node_id)
    # stored fin - node lb
    elif fin_seq is not None and node_fin_seq is None:
        assert cost >= node_fin_cost, (cost, node_fin_cost)
    # stored fin - node fin
    else:
        assert fin_seq is not None and node_fin_seq is not None
        assert cost == node_fin_cost, (cost, node_fin_cost)


def propagate_finishing_sequence(stack, node, ub, seen_asm):
    if not stack:
        log('Root node reached!')
        assert node.fin_cost == ub, (node.fin_cost, ub, node.fin_seq)
        log(node.fin_seq)
        if node.fin_seq is not None:
            rowp = list(chain.from_iterable(node.fin_seq))
            log(rowp)
            _colp = get_hessenberg_form(node.g, node.remaining_eqs, rowp, ub)[1]
            #to_pdf(node.g, rowp, colp, '', 'aaa')
        return
    #---------------------------------------------------
    _upsert(seen_asm, node.asm_key, node.fin_cost, node.fin_seq, node.dbg_id, node.remaining_eqs)
    #---------------------------------------------------
    parent = stack[-1]
    update_parent(node, parent)

def _solve_disconnected(stack, node, eqs, seen_asm, indent, ub, best_rowp):
    log_indent('Disconnected', level=indent)
    subcost, rowp, cost_profiles = \
        solve_pieces(node.g, node.remaining_eqs, seen_asm, indent+1)
    log_indent('Disconnected, solved', len(cost_profiles), 'pieces recursively', level=indent)
    total_cost = node.cost + subcost
    elims = node.elims
    elims.extend(rowp)
    apply_exclusion_rule_on_pieces(stack, node.cost, cost_profiles)
    if total_cost < ub:
        log_indent('Improved UB: {} -> {}'.format(ub, total_cost), level=indent)
        assert set(elims) == eqs
        ub = total_cost
        best_rowp = elims
        set_incumbent(indent, ub, best_rowp)
    return ub, best_rowp, subcost, rowp


def _next_eq(node, indent):
    nzeros, eq = node.heap.popitem()
    guesses = max(0, nzeros-1)
    running_cost = node.cost + guesses
    log_indent('LB at node:', running_cost, level=indent)
    return running_cost, eq, guesses


def _solve_eq(stack, node, indent, running_cost, eq):
    elims, g, remaining_eqs = node.elims[:], fast_copy(node.g), set(node.remaining_eqs)
    log_indent(elims, level=indent)
    rowp, heap = eliminate(g, remaining_eqs, eq)
    elims.extend(rowp)
    log_indent(elims, level=indent)
    apply_exclusion_rule(stack, running_cost, rowp)
    return elims, g, heap, remaining_eqs, rowp


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


def init_dfs_stack(g, eqs, indent):
    # heap key: eq, value: nzeros
    heap = None
    connected = is_connected(g)
    if connected:
        heap = PriorityQueue()
        heap.populate((eq, len(g[eq])) for eq in sorted(eqs))
    remaining_eqs = set(eqs)
    asm_key = to_asm_key(remaining_eqs)
    root_node = Node(0, [ ], g, heap, connected, remaining_eqs, None, asm_key,
                     None, None, get_id())
    register('init', indent)
    return [root_node]


def fast_copy(obj):
    pkl_str = cPickle_dumps(obj, cPickle_HIGHEST_PROTOCOL)
    return cPickle_loads(pkl_str)

#@profile
def eliminate(g, remaining_eqs, eq):
    rowp = [ eq ]
    remaining_eqs.remove(eq)
    heap = PriorityQueue()
    heap.populate((r, len(g[r])) for r in sorted(remaining_eqs))
    elimination_step(g, heap, eq)
    while heap:
        eq, nzeros = heap.peekitem()
        if nzeros > 1:
            break
        remaining_eqs.remove(eq)
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


def solve_pieces(g, remaining_eqs, seen_asm, indent):
    graphs = list(connected_component_subgraphs(g))
    subeqs = [{n for n in subg if n in remaining_eqs} for subg in graphs]
    #print(subeqs)
    assert len(graphs) > 1
    sols = list(_solve_problem(subg, seqs, seen_asm, indent=indent) 
                               for subg, seqs in izip(graphs, subeqs))
    # cost, rowp = sols[i]
    total_cost = sum(t[0] for t in sols)
    c_minus_r = [len(subg)-2*len(seqs) for subg, seqs in izip(graphs, subeqs)]
    keys = [(c_m_r, cost, min(seqs)) for c_m_r, (cost, _), seqs in izip(c_minus_r, sols, subeqs)]
    elims = list(chain.from_iterable(sols[i][1] for i in argsort(keys)))
    return total_cost, elims, get_t_profiles(graphs, sols)


def get_t_profiles(graphs, sols):
    # cost, rowp = sols[i]
    list_g_rowp = list(izip(graphs, (s[1] for s in sols)))
    t_profiles = [t_profile(subg, rowp) for subg, rowp in list_g_rowp]
    assert len(graphs) == len(t_profiles)
    eq_costs = [profile[-1] for profile in t_profiles]
    assert all(s[0] == eq_cost[1] for s, eq_cost in izip(sols, eq_costs))
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

if __name__ == '__main__':
    main()
