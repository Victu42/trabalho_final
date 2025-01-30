import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_polar_data(file_path, color, label):
    """
    Função para plotar dados em coordenadas polares.

    Parâmetros:
        file_path (str): Caminho do arquivo de dados.
        color (str): Cor da linha do gráfico.
        label (str): Legenda do gráfico.
    """
    # Lê o arquivo, ignorando a primeira linha
    data = pd.read_csv(file_path, skiprows=1, delim_whitespace=True, header=None)

    # Assume que a primeira coluna é phi (ângulo) e a segunda é r (raio)
    phi = data[0].values  # Ângulo em radianos
    r = data[1].values    # Raio

    # Plota os dados
    ax.plot(phi, r, label=label, color=color, linewidth=2)

# Cria uma figura com um único gráfico polar
plt.figure(figsize=(8, 8))
ax = plt.subplot(111, projection='polar')

# Plota a órbita (orbit_data.txt)
plot_polar_data("orbit_data.txt", "blue", "Órbita")

# Plota a trajetória da luz (eye_data.txt)
plot_polar_data("eye_data.txt", "red", "Trajetória da Luz")

# Configurações do gráfico
ax.set_title("Órbita e Trajetória da Luz em Coordenadas Polares (r em função de $\phi$)", pad=20)
ax.grid(True)
ax.legend(loc='upper right')

# Exibe o gráfico
plt.show()