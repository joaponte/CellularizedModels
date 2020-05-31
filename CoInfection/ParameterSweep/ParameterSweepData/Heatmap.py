import numpy as np

total_replicates = 20
names = ['time', 'U', 'AI1', 'AI2', 'AD', 'BI1', 'BI2', 'BD', 'VA', 'VB']
initial_number_of_cells = 7396.0

production_multiplier = [0.01,0.10,1.0,10.0,100.0]
diffusion_multiplier = [0.01,0.10,1.0,10.0,100.0]

heatmap_matrix_cells = np.zeros((5,5))
heatmap_matrix_max = np.zeros((5,5))
heatmap_matrix_total = np.zeros((5,5))

for i in range(len(production_multiplier)):
    a = production_multiplier[i] * 100
    for j in range(len(diffusion_multiplier)):
        b = diffusion_multiplier[j] * 100
        m = np.zeros((20,3))
        for k in range(1,21):
            #print('cellularizedmodel_%.5d_%.5d_%i.txt' % (a,b,k))
            file = np.genfromtxt('cellularizedmodel_%.5d_%.5d_%i.txt' % (a,b,k), delimiter=',', names=names, skip_header=1, max_rows=865)

            cells_infected_by_virusA = (file['AI1'][-1] + file['AI2'][-1] + file['AD'][-1])
            cells_infected_by_virusB = (file['BI1'][-1] + file['BI2'][-1] + file['BD'][-1])
            ratio_cells_A_to_B = np.log10(cells_infected_by_virusA / cells_infected_by_virusB)

            max_A = max(file['VA'])
            max_B = max(file['VB'])
            ratio_max_A_to_B = np.log10(max_A / max_B)

            total_A = np.sum(file['VA'])
            total_B = np.sum(file['VB'])
            ratio_total_A_to_B = np.log10(total_A / total_B)

            m[i-1] = (ratio_cells_A_to_B,ratio_max_A_to_B,ratio_total_A_to_B)
        heatmap_matrix_cells[i][j] = np.mean(m,0)[0]
        heatmap_matrix_max[i][j] = np.mean(m,0)[1]
        heatmap_matrix_total[i][j] = np.mean(m,0)[2]

np.savetxt('heatmap_matrix_cells.txt',heatmap_matrix_cells)
np.savetxt('heatmap_matrix_max.txt',heatmap_matrix_max)
np.savetxt('heatmap_matrix_total.txt',heatmap_matrix_total)