import networkx as nx
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps


def normalized_data(raw_data):
    """Normalized Data from raw response of String.

    returns proper data about Protein A, Protein B and with their score.
    """
    lines = raw_data.text.split("\n")
    data = [l.split("\t") for l in lines]
    df = pd.DataFrame(data[1:-1], columns=data[0])
    # Protein A | Protein B | Score
    return df[["preferredName_A", "preferredName_B", "score"]]


def graph_info(G, proteins):
    # Network info about edges (Interactions) and Nodes (Proteins)
    print(f"Number of nodes: {G.number_of_nodes()}")
    print(f"Number of edges: {G.number_of_edges()}")
    print("Avg degree of each protein:")
    for protein in proteins:
        print(
            f"Protein {protein}: {round(nx.average_neighbor_degree(G).get(protein))} - {G.degree(protein)}"
        )


def rescale(l, new_min, new_max):
    # Returns degree of each |V|
    arr = list(l)
    return [
        (x - min(arr)) / (max(arr) - min(arr)) * (new_max - new_min) + new_min
        for x in arr
    ]


def draw_graph(data, proteins, option="normal", mst=False):
    # A graph instance to show interactions of Protein 'A' to Protein 'B'
    G = nx.Graph(name="Protein Interactions Graph")
    interactions = np.array(normalized_data(data))

    for i in range(len(interactions)):
        # ["Pro_A", "Pro_B", "Score"]
        protein_a = interactions[i][0]
        protein_b = interactions[i][1]
        # Score |E| between |V1| and |V2|
        weight = float(interactions[i][2])
        # Connect nodes (Protein A and B) with simple weighted |E|
        G.add_weighted_edges_from([(protein_a, protein_b, weight)])

    graph_info(G, proteins)

    if option == "normal":
        pos = nx.spring_layout(G)
        plt.figure(figsize=(11, 11))
        nx.draw_networkx(G)
        plt.axis("off")
        plt.show()
    elif option == "betweenness_centrality":
        # Matplotlib plasma colormap
        graph_colormap = colormaps.get_cmap("plasma")
        # Each node color config to varies with degree
        node_color = rescale([G.degree(v) for v in G], 0.0, 0.9)
        node_color = [graph_colormap(i) for i in node_color]
        # Node size varies with betweenness centrality - map to range [10, 100]
        bc = nx.betweenness_centrality(G)
        node_size = rescale([v for v in bc.values()], 1500, 7000)
        # Edge width shows 1-weight to convert cost back to strength of
        # interaction
        edge_width = rescale([float(G[u][v]["weight"]) for u, v in G.edges], 0.1, 4)
        # Edge color also demonstrates weight
        edge_color = rescale([float(G[u][v]["weight"]) for v, u in G.edges], 0.1, 1)
        edge_color = [graph_colormap(i) for i in edge_color]

        pos = None
        T = None

        if mst == False:
            pos = nx.spring_layout(G)
        else:
            T = nx.minimum_spanning_tree(G)
            pos = nx.spring_layout(T)

        figure = plt.figure(figsize=(19, 6), facecolor=[0.7, 0.7, 0.7, 0.1])
        figure.canvas.manager.set_window_title(option + f"_{str(mst).lower()}")

        nx.draw_networkx(
            G if mst == False else T,
            pos=pos,
            with_labels=True,
            node_color=node_color,
            node_size=node_size,
            edge_color=edge_color,
            width=edge_width,
            font_color="white",
            font_weight="bold",
            font_size="9",
        )
        plt.axis("off")
        plt.show()
    else:
        print("Option type error")


def main():
    # Top 20 list of proteins selected from a larger set of human proteins
    # related to serotonin.
    proteins_list: [str] = [
        "TPH1",
        "COMT",
        "SLC18A2",
        "HTR1B",
        "HTR2C",
        "HTR2A",
        "MAOA",
        "TPH2",
        "HTR1A",
        "HTR7",
        "SLC6A4",
        "GABBR2",
        "POMC",
        "GNAI3",
        "NPY",
        "ADCY1",
        "PDYN",
        "GRM2",
        "GRM3",
        "GABBR1",
    ]
    # Convert to `.tsv` or Tab Separated Values like `.csv`
    proteins = "%0d".join(proteins_list)

    # Species is kind of humans = 9606
    url = f"https://string-db.org/api/tsv/network?identifiers={proteins}&species=9606"
    try:
        r = requests.get(url)
        draw_graph(r, proteins_list, option="betweenness_centrality", mst=False)
    except requests.exceptions.RequestException as _:
        raise SystemExit("Network error, Please check your internet connection.")


if __name__ == "__main__":
    main()
