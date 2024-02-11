import networkx as nx

G = nx.Graph()

G.add_edge(1, 2)
G.add_edge(1, 3)
G.add_edge(2, 3)

# Graph shape of NetworkX
pos = nx.spring_layout(G)
nx.draw_networkx(G)
plt.axis("off")
plt.show()
