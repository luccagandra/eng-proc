import networkx as nx
import matplotlib.pyplot as plt

def generate_image():
  # Cria um grafo de grade 2 por 3 por 4
  G = nx.grid_graph(dim=(2, 2))
  # Desenha o grafo com os rótulos dos nós
  nx.draw(G, with_labels=True)
  # Mostra o diagrama
  plt.show()

generate_image()