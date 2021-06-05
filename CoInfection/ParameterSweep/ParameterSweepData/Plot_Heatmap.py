import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

data_to_plot = 3
production_multiplier = [0.01,0.10,1.0,10.0,100.0]
diffusion_multiplier = [0.01,0.10,1.0,10.0,100.0]

heatmap_matrix_cells = np.genfromtxt('heatmap_matrix_cells.txt')
heatmap_matrix_max = np.genfromtxt('heatmap_matrix_max.txt')
heatmap_matrix_total = np.genfromtxt('heatmap_matrix_total.txt')

if data_to_plot == 1:
    sns.heatmap(heatmap_matrix_cells, cmap="jet", xticklabels=production_multiplier, yticklabels = diffusion_multiplier)
    plt.title('Cells infected ratio (A/B)')
    plt.ylabel('Virus diffusion ratio (A/B)')
    plt.xlabel('Virus production ratio (A/B)')
    plt.savefig("Cell_Ratio.pdf")

if data_to_plot == 2:
    sns.heatmap(heatmap_matrix_max, cmap="jet", xticklabels=production_multiplier, yticklabels = diffusion_multiplier)
    plt.title('Max virus ratio (A/B)')
    plt.ylabel('Virus diffusion ratio (A/B)')
    plt.xlabel('Virus production ratio (A/B)')
    plt.savefig("Max_Virus_Ratio.pdf")

if data_to_plot == 3:
    sns.heatmap(heatmap_matrix_total, cmap="jet", xticklabels=production_multiplier, yticklabels = diffusion_multiplier)
    plt.title('Total virus ratio (A/B)')
    plt.ylabel('Virus diffusion ratio (A/B)')
    plt.xlabel('Virus production ratio (A/B)')
    plt.savefig("Total_Virus_Ratio.pdf")

plt.show()