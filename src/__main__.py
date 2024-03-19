import networkx as nx
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import _cm

# Top 20 list of proteins selected from a larger set of human proteins related
# to serotonin.
protein_list: [str] = [
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

proteins = "%0d".join(protein_list)

# Species is kind of humans = 9606
url = f"https://string-db.org/api/tsv/network?identifiers={proteins}&species=9606"
r = requests.get(url)

lines = r.text.split("\n")
# Makes component based on tabs
data = [l.split("\t") for l in lines]
df = pd.DataFrame(data[1:-1], columns=data[0])
interactions = df[["preferredName_A", "preferredName_B", "score"]]

# Make a graph instance to show interactions of Protein 'A' to Protein 'B'
G = nx.Graph(name="Protein Interactions Graph")

interactions = np.array(interactions)
for i in range(len(interactions)):
    interaction = interactions[i]
    # Protein 'A' node
    a = interaction[0]
    # Protein 'B' node
    b = interaction[1]
    # Weight
    w = float(interaction[2])
    # Add weighted edge to graph
    G.add_weighted_edges_from([(a, b, w)])


# Network info about edges (Interactions) and Nodes (Proteins)
# print(f"Number of nodes: {G.number_of_nodes()}")
# print(f"Number of edges: {G.number_of_edges()}")
# print(f"Average degree: {G.degree}")

# Graph shape of NetworkX
pos = nx.spring_layout(G)
plt.figure(figsize=(11, 11))

nx.draw_networkx(G)
plt.axis("off")
plt.show()
