import matplotlib.pyplot as plt
import networkx as nx


def plot_graph_components(nodes, edges, output_file_name='labels_and_colors.png'):
    G = nx.Graph()
    labels = {}
    for node in nodes:
        labels[node] = r'$  %s$' % str(node)
        G.add_node(node)
    for u, v in edges:
        G.add_edge(u, v)

    try:
        from networkx import graphviz_layout
        layout = nx.graphviz_layout
    except ImportError:
        print("PyGraphviz not found; drawing with spring layout; will be slow.")
        layout = nx.spring_layout
    n = len(nodes)
    p = 0.001
    pos = layout(G)
    nx.draw(G, pos,
            with_labels=False,
            node_size=100
            )
    Gcc = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)
    G0 = Gcc[0]
    nx.draw_networkx_edges(G0, pos,
                           with_labels=False,
                           edge_color='black',
                           width=1.0
                           )
    # show other connected components
    for Gi in Gcc[1:]:
        if len(Gi) > 1:
            nx.draw_networkx_edges(Gi, pos,
                                   with_labels=False,
                                   edge_color='black',
                                   alpha=0.9,
                                   width=1.0
                                   )

    # nx.draw_networkx_nodes(G,pos, node_size=2, alpha=0.6)
    # nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.9)
    # nx.draw_networkx_labels(G,pos,labels,font_size=14)
    plt.axis('off')
    plt.savefig(output_file_name)
