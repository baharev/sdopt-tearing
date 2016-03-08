# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>

def dummy(_, prog=None): pass

def plot(directed_graph, prog='dot'):
    from matplotlib import pyplot as plt
    from networkx import graphviz_layout, draw_networkx
    positions = graphviz_layout(directed_graph, prog)
    draw_networkx(directed_graph, pos=positions, node_size=800)
    mng = plt.get_current_fig_manager()
    # TODO Post a wrapper to Code Review?
    mng.full_screen_toggle()
    #mng.resize(1865,1025)
    plt.show()

def plot_tearing_result(G, elims, prog='neato', font_size=12, node_size=800):
    from matplotlib import pyplot as plt
    import networkx as nx
    # graphviz_layout bug: http://stackoverflow.com/a/35280794/341970
    from networkx.drawing.nx_agraph import graphviz_layout
    from sys import stderr
    
    ax = plt.figure().add_subplot(111)
    try:
        pos = graphviz_layout(G, prog)
    except IndexError:
        stderr.write('Known bug with pydot, see ' 
                     'https://github.com/networkx/networkx/issues/1836\n')
        return 
    nodes = G.nodes()
    
    nx.draw_networkx_nodes(G, pos, nodes, node_color='r', node_size=node_size)
    labels = { n : str(n) for n in nodes }
    nx.draw_networkx_labels(G, pos, labels, font_size=font_size)
    
    torn_set = set(elims) 
    kept = [ e for e in G.edges_iter() if e not in torn_set ]
    
    nx.draw_networkx_edges(G, pos, edgelist=elims, edge_color='red')
    
    nx.draw_networkx_edges(G, pos, edgelist=kept, edge_color='black')
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    mng = plt.get_current_fig_manager()
    #mng.full_screen_toggle()
    mng.resize(1865,1025)  
    plt.show()
