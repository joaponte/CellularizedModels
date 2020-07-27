# Data control options
plot_pop_data_freq = 10  # Plot population data frequency (disable with 0)
write_pop_data_freq = 0  # Write population data to simulation directory frequency (disable with 0)
plot_med_diff_data_freq = 10  # Plot total diffusive field amount frequency (disable with 0)
write_med_diff_data_freq = 0  # Write total diffusive field amount frequency (disable with 0)
plot_spat_data_freq = 0  # Plot spatial data frequency (disable with 0)
write_spat_data_freq = 0  # Write spatial data to simulation directory frequency (disable with 0)
plot_death_data_freq = 10  # Plot death data frequency (disable with 0)
write_death_data_freq = 0  # Write death data to simulation directory frequency (disable with 0)

plot_ode_sol = True  # Set to True to instantiate ODE model and plot, for comparison with spatial model

# Conversion Factors
s_to_mcs = 5 * 60  # s/mcs
um_to_lat_width = 2.0  # um/lattice_length

# Experimental Parameters
exp_cell_diameter = 10.0  # um

exp_virus_dc = 0.025  # um^2/s

exp_cytokine_dc_cyto = 0.1  # um^2/s; estimated diffusion constant in cytoplasm (B)
# ^ the  /10 is not experimental; added because of relative small area and because virus D is (or was) slowed down

# =============================
# CompuCell Parameters
# cell
cell_diameter = exp_cell_diameter / um_to_lat_width
cell_volume = cell_diameter ** 2
volume_lm = 9

# virus diffusion
virus_dc = exp_virus_dc * s_to_mcs / (um_to_lat_width ** 2)  # virus diffusion constant

cytokine_dc = exp_cytokine_dc_cyto * s_to_mcs / (um_to_lat_width ** 2)  # CK diff cst

# Lambda Chemotaxis
lamda_chemotaxis = 1E5

# Antimony/SBML model step size
sbml_step_size = 1.0

# Regulation shape parameter
shape_param = 10

# Initial fraction of infected cells
init_infect = 0.05
