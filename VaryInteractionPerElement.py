########################################################################
# VaryInteractionPerElement.py
# Julie Butler
# July 17th, 2018
# Version 1.0
#
# Explores the variation between H(s) calcualted using the Magnus Expansion and H(s)
# calculated using the derivative method, and how the variation varies as the interaction
# strength varies
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

# Range of energy level spacings from -1.0 to 1.0 by 0.01
g_range = arange (-1.0, 1.2, 0.2)

# Fixed energy level spacing
# Values taken from srg_pairing.
fixed_d = 1.0

# Fixed value of the threshold to be used
threshold_value = 0.001

# To hold the final values
discrepancy_vary_g = []

# Collects the data for varying the energy spacing
for g in g_range:
	srgHs = srg_main (flowparams, fixed_d, g)
	meHs = magnus_expansion_main (flowparams, threshold_value, fixed_d, g)
	discrepancy = []
	for i, j in zip (meHs, srgHs):
		discrepancy.append (diff (i, j, 6))
	discrepancy_vary_g.append (discrepancy)

# Gathers the data 
interaction_negative_1 = discrepancy_vary_g[0][3:7] 
interaction_negative_0_2 = discrepancy_vary_g[4][3:7] 
interaction_0_2 = discrepancy_vary_g[6][3:7] 
interaction_1 = discrepancy_vary_g[10][3:7] 

#################################################
#
#                         GRAPHING
#
#################################################
# Imports here to prevent error

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
#                  Interaction Strength: -1.0
#################################################
make_single_figure (1, interaction_negative_1, "Interaction Strength: -1.0",
     subplot_titles, index_list, cmap, norm)

#################################################
#                                Figure 2
#                       Interaction Strength: -0.2
#################################################
make_single_figure (2, interaction_negative_0_2, "Interaction Strength: -0.2", 
    subplot_titles, index_list, cmap, norm)

#################################################
#                                Figure 3
#                        Interaction Strength: 0.2
#################################################
make_single_figure (3, interaction_0_2, "Interaction Strength: 0.2", 
    subplot_titles, index_list, cmap, norm)

#################################################
#                                Figure 4
#                       Interaction Strength: 1.0
#################################################
make_single_figure (4, interaction_1, "Interaction Strength: 1.0", 
    subplot_titles, index_list, cmap, norm)

#################################################
#                                         Show
#################################################
# Show the plot
show()
