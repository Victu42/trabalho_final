import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_polar_data(file_path, title, color, label):

    # Lê o arquivo, ignorando a primeira linha
    data = pd.read_csv(file_path, skiprows=1, delim_whitespace=True, header=None)

    # Assume que a primeira coluna é phi (ângulo) e a segunda é r (raio)
    phi = 180*data[0].values/(2*pi)  # Ângulo em radianos
    r = data[1].values    # Raio

    # Cria o gráfico em coordenadas polares
    ax.plot(phi, r, label=label, color=color, linewidth=2)
    ax.set_title(title, pad=20)
    ax.grid(True)
    ax.legend(loc='upper right')
    ax.set_rmax(np.max(r)+0.1*np.max(r))

# Cria uma figura com dois subplots lado a lado
plt.figure(figsize=(16, 8))

# Primeiro gráfico: Órbita (orbit_data.txt)
ax = plt.subplot(121, projection='polar')
plot_polar_data("orbit_data.txt", "Órbita em Coordenadas Polares (r em função de $\phi$)", "blue", "Órbita")

# Segundo gráfico: Trajetória da luz (eye_data.txt)
ax = plt.subplot(122, projection='polar')
plot_polar_data("eye_data.txt", "Trajetória da Luz em Coordenadas Polares (r em função de $\phi$)", "red", "Trajetória da Luz")

# Ajusta o layout e exibe os gráficos
plt.tight_layout()
plt.show()