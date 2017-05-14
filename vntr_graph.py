import matplotlib.pyplot as plt
import networkx as nx
from Bio import pairwise2


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


def get_nodes_and_edges_of_vntr_graph(vntrs):
    edges = []
    nodes = []
    for i in range(len(vntrs)):
        nodes.append(vntrs[i].id)
        for j in range(i + 1, len(vntrs)):
            if vntrs[i].is_homologous_vntr(vntrs[j]):
                edges.append((vntrs[i].id, vntrs[j].id))

    nodes = [1, 2, 3, 4, 5, 8, 9, 10, 12, 16, 17, 18, 19, 21, 22, 24, 25, 28, 29, 30, 31, 32, 33, 34, 38, 40, 47, 53,
             57, 59, 66, 67, 68, 69, 70, 71]
    edges = [(1, 8), (1, 16), (2, 17), (4, 18), (8, 16), (30, 32), (30, 33), (32, 33), (34, 40), (34, 47), (38, 57),
             (38, 59), (38, 67), (40, 47), (57, 59), (57, 67), (59, 67)]

    return nodes, edges
