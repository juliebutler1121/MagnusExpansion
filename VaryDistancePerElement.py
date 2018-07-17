########################################################################
# VaryDistancePerElement.py
# Julie Butler
# July 17th, 2018
# Version 1.0
#
# Explores the variation between H(s) calcualted using the Magnus Expansion and H(s)
# calculated using the derivative method, and how the variation varies as the energy level
# spacing varies# Explores the variation between H(s) calcualted using the Magnus Expansion and H(s)
# calculated using the derivative method, and how the variation varies as the energy level
# spacing varies
########################################################################

########################################################################
# Function Outline:
########################################################################

#################################################
#
#                              IMPORTS
#
#################################################
# Third-Party Imports
from numpy import array, arange
from pylab import *
import matplotlib as mpl
# Local Imports
from MagnusExpansionButler import main as magnus_expansion_main
from srg_pairing import main as srg_main
from MatrixComparisons import compare_square_matrices_by_element as diff
from GraphingCapabilities import set_up_graph
from GraphingCapabilities import make_4_subplots_per_figure as make_single_figure

#################################################
#
#                        Discrepancy calculations
#
#################################################
# flow parameters for snapshot images
flowparams = array([0.,0.001,0.01,0.05,0.1, 1., 5., 10.])

# Range of energy level spacings from 0 to 2 by 0.2
d_range = arange (0.0, 2.2, 0.2)

# Fixed interaction value
# Values taken from srg_pairing.
fixed_g = 0.5

# Fixed value of the threshold to be used
threshold_value = 0.001

# To hold the final values
discrepancy_vary_d = []

# Collects the data for varying the energy spacing
for d in d_range:
	srgHs = srg_main (flowparams, d, fixed_g)
	meHs = magnus_expansion_main (flowparams, threshold_value, d, fixed_g)
	discrepancy = []
	for i, j in zip (meHs, srgHs):
		discrepancy.append (diff (i, j, 6))
	discrepancy_vary_d.append (discrepancy)

# Gathers the data 
distance_0_2 = discrepancy_vary_d[1][3:7]
distance_0_6 = discrepancy_vary_d[3][3:7]
distance_1 = discrepancy_vary_d[5][3:7]
distance_1_6 = discrepancy_vary_d[8][3:7]
distance_2 = discrepancy_vary_d[10][3:7]

#################################################
#
#                         GRAPHING
#
#################################################
# Create discrete colormap
cmap = mpl.colors.ListedColormap(['0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7',
							'0.8', '0.9', '1.0'])
bounds = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
# Titles for the subplots
subplot_titles = ["Flow Parameter: 0.05", "Flow Parameter: 0.1", "Flow Parameter: 1.0",
                            "Flow Parameter: 5.0"]
index_list = [0, 1, 2, 3]

#################################################
#                              Figure 1
#                   Energy Level Spacing: 10.0
#################################################
make_single_figure (1, distance_0_2, "Energy Level Spacing: 0.2", 
    subplot_titles, index_list, cmap, norm)

#################################################
#                                Figure 2
#                      Energy Level Spacing: 1.0
#################################################
make_single_figure (2, distance_0_6, "Energy Level Spacing: 0.6", 
    subplot_titles, index_list, cmap, norm)

#################################################
#                                Figure 3
#                     Energy Level Spacing: 0.1
#################################################
make_single_figure (3, distance_1, "Energy Level Spacing: 1.0", 
    subplot_titles, index_list, cmap, norm)

#################################################
#                                Figure 4
#                     Energy Level Spacing: 0.01
#################################################
make_single_figure (4, distance_1_6, "Energy Level Spacing: 1.6", 
    subplot_titles, index_list, cmap, norm)

#################################################
#                                Figure 5
#                     Energy Level Spacing: 0.001
#################################################
make_single_figure (5, distance_2, "Energy Level Spacing: 2.0",
     subplot_titles, index_list, cmap, norm)

#################################################
#                                         Show
#################################################
# Show the plot
show()
