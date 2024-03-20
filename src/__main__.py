import networkx as nx
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import _cm


def normalized_data(raw_data):
    """Normalized Data from raw response of String.

    returns proper data about Protein A, Protein B and with their score.
    """
    lines = raw_data.text.split("\n")
    data = [l.split("\t") for l in lines]
    df = pd.DataFrame(data[1:-1], columns=data[0])
    # Protein A | Protein B | Score
    return df[["preferredName_A", "preferredName_B", "score"]]


def draw_graph(data, proteins):
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

    # Network info about edges (Interactions) and Nodes (Proteins)
    print(f"Number of nodes: {G.number_of_nodes()}")
    print(f"Number of edges: {G.number_of_edges()}")
    print("Avg degree of each protein:")
    for protein in proteins:
        print(f"Protein {protein}: {nx.average_neighbor_degree(G).get(protein)}")

    pos = nx.spring_layout(G)
    plt.figure(figsize=(11, 11))
    nx.draw_networkx(G)
    plt.axis("off")
    plt.show()


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
        draw_graph(r, proteins_list)
    except requests.exceptions.RequestException as e:
        raise SystemExit("Network error, Please check your internet connection.")


if __name__ == "__main__":
    main()
