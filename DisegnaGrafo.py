import os

import networkx as nx
import matplotlib.pyplot as plt


def draw_labeled_graph(labeled_graph, graph_name, labels):
    """
    Disegna un grafo etichettato. Ogni nodo e ogni spigolo sono visualizzati,
    e le etichette degli spigoli vengono mostrate sopra di essi.

    :param labeled_graph: Dizionario che rappresenta il grafo etichettato,
                          nella forma {v: {w: label}} dove v e w sono i nodi e label è l'etichetta dello spigolo (classe di equivalenza).
    """
    G = nx.Graph()
    max_length = max(len(bin(label)[2:]) for label in labels.values())
    label_dict = {v: bin(labels[v])[2:].zfill(max_length) for v in labels}
    bit_difference = 0
    for v in labeled_graph:
        for w in labeled_graph[v]:
            if not G.has_edge(v, w):
                v_label = label_dict[v]
                w_label = label_dict[w]
                for i in range(max_length - 1, -1, -1):
                    if v_label[max_length - i - 1] != w_label[max_length - i - 1]:
                        bit_difference = i
                G.add_edge(v, w, label=bit_difference)
                bit_difference = 0
    pos = nx.kamada_kawai_layout(G)


    composed_labels = {}
    for v in G:
        composed_labels[v] = str(str(v) + "\n" + str(label_dict[v]))
    plt.figure()
    nx.draw(G, pos, with_labels=False, node_color='lightblue', edge_color='gray', node_size=850, font_size=10)
    nx.draw_networkx_labels(G, pos, labels=composed_labels, font_size=7)
    edge_labels = nx.get_edge_attributes(G, "label")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='red', font_size=7, label_pos=0.5)

    plt.draw()
    if not os.path.exists("Grafi"):
        os.mkdir("Grafi")
    plt.savefig("Grafi/" + graph_name)
    plt.show()

    return max_length