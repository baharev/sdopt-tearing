# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from itertools import count
import six
from networkx import DiGraph, Graph, is_connected
from networkx.algorithms.bipartite import is_bipartite_node_set
#from plot import plot
from py3compat import irange

# From: Partitioning and Tearing of Networks Applied to Process Flowsheeting;
# Truls Gundersen and Terje Hertzberg; Modeling, Identification and Control; 
# 1983, Vol 4, No 3, pp. 139-165. DOI: 10.4173/mic.1983.3.2  
# http://www.mic-journal.no/PDF/1983/MIC-1983-3-2.pdf


def gen_benchmarks_as_undirected_bipartite():
    # The system is artificially made square
    for dig in gen_benchmark_digraphs():
        if dig.number_of_selfloops():
            print('Skipping it, as this digraph has self-loops\n')
            continue # Or we could just remove the self-loops...
        yield digraph_to_undirected_bipartite(dig)

def digraph_to_undirected_bipartite(dig):
    # The system is artificially made square
    #plot(dig, prog='sfdp')
    # We will introduce fake nodes (equations) to make the system square
    nodeid = count()
    eq_id = { n : [next(nodeid) for _ in irange(dig.node[n]['weight'])] \
              for n in sorted(dig) }
    #
    edgeid = count()
    var_id = { e : 'x%d' % next(edgeid) for e in sorted(dig.edges_iter()) }
    # Build the bipartite graph: edges of dig become the var node set, and 
    # each var is connected to its equation(s).
    g = Graph(name=dig.graph['name'])
    for u,v in dig.edges_iter():
        var = var_id[(u,v)]
        append_eqs_for_var(g, eq_id[u], var)
        append_eqs_for_var(g, eq_id[v], var)
    #plot(g, prog='sfdp')
    eqs = set()
    for eq_list in six.itervalues(eq_id):
        eqs.update( eq_list )
    forbidden = set()
    return g, eqs, forbidden

def append_eqs_for_var(g, eqs, var):
    for eq in eqs:
        assert not g.has_edge(eq, var), (eq, var)
        g.add_edge(eq, var)    

#-------------------------------------------------------------------------------

def digraph_as_rectangular_bipartite(dig):
    nodeid = count()
    eq_id = { n : next(nodeid) for n in sorted(dig) }
    edgeid = count()
    var_id = { e : 'x%d' % next(edgeid) for e in sorted(dig.edges_iter()) }
    # Build the bipartite graph: edges of dig become the var node set, and 
    # each var is connected to its equation(s).
    g = Graph(name=dig.graph['name'])
    for u, v in dig.edges_iter():
        eq1 =  eq_id[u]
        eq2 =  eq_id[v]
        var = var_id[(u,v)]
        assert not g.has_edge(eq1, var), (eq1, var)
        assert not g.has_edge(eq2, var), (eq2, var)        
        g.add_edge(eq1, var)
        g.add_edge(eq2, var)
    eqs = set( six.itervalues(eq_id) )
    assert is_connected(g)
    assert is_bipartite_node_set(g, eqs)   
    set_lower_bound(g, eqs)
    return g, eqs

def set_lower_bound(g, eqs):
    ncons = len(eqs)
    nvars = len(g)-ncons
    dof   = nvars-ncons
    assert dof > 0, dof # expected an underdetermined system
    min_col_count = min(len(g[n]) for n in g if n not in eqs)
    assert min_col_count >= 2, min_col_count
    lb = dof + min_col_count-1
    g.graph['trivial_lb'] = lb
    print('Size: {}x{}, dof = {}, lb = {}'.format(ncons, nvars, dof, lb))

def gen_digraphs_as_rectangular_bipartite():
    for dig in gen_benchmark_digraphs():
        if dig.number_of_selfloops():
            print('Skipping it, as this digraph has self-loops\n')
            continue # Or we could just remove the self-loops...
        yield digraph_as_rectangular_bipartite(dig)

#-------------------------------------------------------------------------------

def gen_benchmark_digraphs():
    for name in sorted(TESTPROBLEMS):
        yield create_testproblem(name)
        print('Finished:', name)
        print() 

def create_testproblem(name):
    print('===================================================================')
    print('Creating:', name)
    return create_graph(name, TESTPROBLEMS[name])

def create_graph(name, edges_str):
    edges_by_src = [ ]
    for line in edges_str.splitlines():
        line = line.strip()
        if line:
            edges_by_src.append([int(i) for i in line.split()])
    g = DiGraph(name=name)
    for edges in edges_by_src:
        src = edges[0]
        for target in edges[1:]:
            assert not g.has_edge(src, target), (src, target)
            g.add_edge(src, target, {'weight':1, 'orig_edges':[(src,target)]})
    # setup node weights
    for n, d in g.nodes_iter(data=True):
        d['weight'] = g.in_degree(n) 
    return g

#-------------------------------------------------------------------------------

TESTPROBLEMS = {
                
'Subproblem 8 (opt=3)' : '''
2  4  6  24  
4  5  6  
5  6  7  
6  27  
7  8  18  
8  5  14  27  
14  7  
16  14  
17  16  22  
18  17  21  
20  17  
21  20  22  
22  2  
24  16  20  26  
26  21  
27  7  26  
''',

'Heavy water subgraph (opt=6)' : '''
1  16  20  
4  1  
11  4  
12  4  29  
13  12  25  
16  13  
17  11  16  
18  11  17  
20  18  
21  16  38  
22  21  34  
25  22  
26  20  25  
27  20  26  
29  27  
30  25  48  
31  30  43  
34  31  
35  29  34  
36  29  35  
38  36  
39  34  57  
40  39  53  
43  40  
44  38  43  
45  38  44  
48  45  
49  43  57  
53  49  
54  48  53  
55  48  54  
57  55
''',

'Problem 1 (opt=15)' : '''
1       2   3   4   5   6
2   1       3   4   5   6
3   1   2       4   5   6
4   1   2   3       5   6
5   1   2   3   4       6
6   1   2   3   4   5    
''',

# FIXME - Self-loops are not handled correctly (grb_* for example)
#       - (opt=number) is not checked
#       - Document the IIS method 0: it is not only faster, it is also better
#       - Increased cutoff in get all cycles. Can it backfire?

# 'Complete graph n=10 (opt=45)' : '''
# 1       2   3   4   5   6   7   8   9   10
# 2   1       3   4   5   6   7   8   9   10
# 3   1   2       4   5   6   7   8   9   10
# 4   1   2   3       5   6   7   8   9   10
# 5   1   2   3   4       6   7   8   9   10
# 6   1   2   3   4   5       7   8   9   10
# 7   1   2   3   4   5   6       8   9   10
# 8   1   2   3   4   5   6   7       9   10
# 9   1   2   3   4   5   6   7   8       10
# 10  1   2   3   4   5   6   7   8   9   
# ''',

'Problem 2 (opt=2)' : '''
 1   2   4   6
 2   3
 3   4   8
 4   5
 5   1  10
 6   4   7
 7   2   8  11
 8   9
 9   4  11
10   6
11  12
12   4  10
''',

'Problem 3 (opt=6)' : '''
 1  2  3   
 2  3  6   
 3  4  5   
 4  1  2  5   
 5  1  2  7   
 6  7  8   
 7  8   
 8  9   10  13   
 9  10   6   7  12   
10   6   7  11
11  12  13  
12  13  
13  14  15  
14  11  12 15  
15  11  12  
''',

'Problem 4 (opt=6)' : '''
 1   5
 2   5
 3   7 
 4   1   5   8
 5   6
 6   2   3   7  10  12
 7  13  19
 8  14
 9   4
10   9
11   6
12   7  17
13   7
14  15
15  16
16  10
17  11  16  18
18   7  17
19  13  18
''',

'Problem 5 (opt=3)' : '''
 1   2 
 2   3
 3   4  11
 4   5 
 5   6
 6   7  12
 7   8
 8   9
 9  10
10  17  20
11  21
12  11
13  15  23
14  13
15  18
16   8  15
17  16
18  19
19  20
20   1
21  22  20
22  23
23  24
24  25
25  14  20
''',

'Problem 6 (opt=5)' : '''
 1   2
 2   3   7
 3   4   5   6  10 
 4  25
 5   3
 6  24
 7   1   9
 8   1
 9  10  14
10  11  15
11  12
12  26
13   9
14  27  28
15  13  16  17  18
16  26
17  15
18  19
19  20  26
20  18  21
21  29
22   1
23   3
26   8
''',

'Problem 7 (opt=3)'  : '''
 1   2   5   
 2   3   
 3   4    
 4   5  14   
 5   6  
 6   8  
 7   8  
 8   9  
 9  10   
10  11  20
11   1  12  
12  13  19  
13   2  15  23  
14  15  19    
15  16  
16  17  
17  18  19  
18  19  23
19   7  
20  29  
21  12  
22  21  
23  24  
24  25  
25  26  19  
26  22  
27  29  
28  27  
29  28  30  
30  22  
''',

'Problem 8 (opt=5)' : '''
 1   2  33
 2   3   6  24
 3   4 
 4   5   6 
 5   6   7
 6  25
 7   8  18
 8   5   9  27
 9  10
10  11
11  12
12  13
13  14 
14   7  29
15  14
16  15
17  16  22
18  17  21
19  17
20  19
21  20  22
22  23
23   1
24  16  20  26
25  27
26  28
27   7  26
28  21
29  30  41
30  31
31  33
32  31
33  34
34  32  35
35  36
36   1  37
37  40
38  33  37
39  38
40  30  39
41  40
''',

'Problem 9 (opt=8)' : '''
 1   9   
 2   8  
 3  12   
 4  13  42    
 5  17  36   
 6  17   
 7   8  13  16
 8   9  12  
 9  10   
10  11  
11  13  
12  16  25  
13  14  
14  15  
15  17  
16  17  
17  18  21  47  
18  19  20  
19  47  49  
20  50  
21  22  25  
22  23  
23  24  
24  25  47    
26  28  29  50  
27  29  32  
28  26  39  
29  26  27  30  31  
30  28  29  
31  29  32  
32  31  35  
33  30  
34  35  38   
35  25  34  37  40  
36  25  35  
37  42  
38  30  36    
40  41  
41  42  47
42  44  
43  46  47  
44  43  
45  46  
46  47    
48  13  42  
49  47  
50  26  27
''',

'Problem 10 (opt=12)' : '''
 1  16  20
 2   1 
 3   2 
 4   3
 5   4   8
 6   5   7 
 7   8  47
 8   9
 9  10
11   6
12   4  29
13  12  25
14  13 
15  14 
16  15 
17  16  11
18  17  19 
19  11  47
20  18
21  16  38
22  34  21
23  22
24  23
25  24
26  25  20
27  26  28
28  20  47
29  27
30  25  48
31  30  43 
32  31 
33  32 
34  33 
35  34  29
36  35  37
37  29  47
38  36
39  34  57
40  39  53
41  40
42  41
43  42
44  43  38
45  44  46 
46  38  47
47  65
48  45
49  43  67
50  49  61
51  50
52  51
53  52
54  53  48 
55  54  56
56  48  65
57  55
58  73
59  58 
60  59
61  60
62  61  57
63  62  64
64  57  65
65  76
66  53  89
67  63
68  66  74
69  75  76
70  84
71  70
72  71  68
73  72
74  73  67
75  74
76  77
77  67  98
78  75  76 
79  86  104
80  72
81  79  87
82  81
83  82
84  83  85
85  80  86
86  101
87  84  80
88  89
89  78  69
90  88  99 
91  90  94
92  91
93  92
94  93  95
95  94  88
96  95  97
97  89  98 
98  107
99  100  102
100  89  101
101  96 
102  109
103  100 102
104  103 
105  104  106
106  101
107  108
108   89
109  108
'''

}
