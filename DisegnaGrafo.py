import networkx as nx
import matplotlib.pyplot as plt


def draw_labeled_graph(labeled_graph):
    """
    Disegna un grafo etichettato. Ogni nodo e ogni spigolo sono visualizzati,
    e le etichette degli spigoli vengono mostrate sopra di essi.

    :param labeled_graph: Dizionario che rappresenta il grafo etichettato,
                          nella forma {v: {w: label}} dove v e w sono i nodi e label è l'etichetta dello spigolo (classe di equivalenza).
    """
    G = nx.Graph()
    for v in labeled_graph:
        for w, label in labeled_graph[v].items():
            if not G.has_edge(v, w):
                G.add_edge(v, w, label=label)

    pos = nx.kamada_kawai_layout(G)

    nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=350, font_size=7,
            font_weight='bold')

    edge_labels = nx.get_edge_attributes(G, 'label')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='red', font_size=5, label_pos=0.5)

    plt.title("Grafo etichettato")
    plt.show()
