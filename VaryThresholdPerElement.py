########################################################################
# VaryThresholdPerElement.py
# Julie Butler
# July 11th, 2018
# Version 1.0
#
# A supporting code which runs the SRG code from srg_pairing.py and MagnusExpansionButler.py
# Calculates the matrix H(s) at a variety of different flow parameters.  For the magnus expansion, it
# also calculates the matrices using a variety of threshold values for loop termination.
########################################################################

########################################################################
# Function Outline:
# ratio_matrix (matrix): Return a matrix such that each element has been divided by the global maximum of 
#   the matrix
# make_single_figure (figure_number, discrepancy_list, graph_title, cmap, norm): Makes the four
#   subplots for a single threshold value
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

#################################################
#
#                         ratio_matrix
#
#################################################
def ratio_matrix (matrix):
    """
    Return a matrix such that each element has been divided by the global maximum of 
    the matrix.
    Input:
        matrix (a matrix): The starting matrix.
    Output:
        matrix (a matrix): The matrix after the ratio calcualtions have been performed
    """
    # Find the maximum value
    local_maxes = [max(sub_array) for sub_array in matrix]
    global_max = max (local_maxes)
    # Convert the elements of the matrix to ratios
    for i in range (0, 6):
    	for j in range (0, 6):
    		matrix[i][j] = matrix[i][j] / global_max
    return matrix

#################################################
#
#                      make_single_figure
#
#################################################
def make_single_figure (figure_number, discrepancy_list, graph_title, cmap, norm):
    """
    Makes the four subplots for a single threshold value
    Input:
        figure_number (an integer): the number of the figure
        discrepancy_list (an array): Contains the data to be graphed
        graph_title (a string): The title to be printed at the top of the graph
        cmap and norm (matplotlib objects): handle the color mapping 
    """
    figure (figure_number)
    # Set the figure title
    suptitle (graph_title, fontsize=40, fontweight='bold')
    # Titles for the subplots
    subplot_titles = ["Flow Parameter: 0.05", "Flow Parameter: 0.3", "Flow Parameter: 1.0",
                            "Flow Parameter: 5.0"]
    index_list = [1, 6, 20, 100]
    # Setting up the four subplots
    gs = plt.GridSpec(2, 2)
    for i in range (0, 2):
        iter_number = 2*i
        for j in range (0, 2):
            # Select the correct subplot
            plt.subplot (gs [i, j])
            # Set the subplot title
            title (subplot_titles [iter_number], fontsize=20, fontweight='bold')
            # Collect the data
            matrix = discrepancy_list [index_list [iter_number]]
            matrix = ratio_matrix (matrix)
            # Plot the matrix using the created color map
            imshow(matrix, interpolation='nearest', 
                     extent=[0.5, 0.5+6, 0.5, 0.5+6],
                     cmap=cmap,
            		 norm=norm)
            # Create the color bar legend
            colorbar()
            set_up_graph ("", "", False)
            iter_number += 1


#################################################
#
#                        Discrepancy calculations
#
#################################################
# Gets a list of numbers to be used for flow parameters to solve the differential equation.  current
# values range from 0 to 10 with a spacing of 0.05, leading to 200 total flow parameters used to 
# solve the differential equations
flowparams_new = arange (0, 10, 0.05)

# Getting the H(s) matrices using the Magnus expansion for a variety of threshold values, given in 
# the parameters.
# Naming Convention: meHs10 and meHs1 refer directly to the threshold used.  meHsex refers to a 
# threshold of 10e-x.
meHs10 = magnus_expansion_main (flowparams_new, 10, 1, 0.5)
meHs1 = magnus_expansion_main (flowparams_new, 1, 1, 0.5)
meHse1 = magnus_expansion_main (flowparams_new, 0.1, 1, 0.5)
meHse2 = magnus_expansion_main (flowparams_new, 0.01, 1, 0.5)
meHse3 = magnus_expansion_main (flowparams_new, 0.001, 1, 0.5)
meHse4 = magnus_expansion_main (flowparams_new, 0.0001, 1, 0.5)
meHse5 = magnus_expansion_main (flowparams_new, 0.00001, 1, 0.5)
meHse6 = magnus_expansion_main (flowparams_new, 0.000001, 1, 0.5)
meHse7 = magnus_expansion_main (flowparams_new, 0.0000001, 1, 0.5)
meHse8 = magnus_expansion_main (flowparams_new, 0.00000001, 1, 0.5)
meHse9 = magnus_expansion_main (flowparams_new, 0.000000001, 1, 0.5)

# Gets the H(s) matrices using the derivative method 
srgHs = srg_main (flowparams_new, 1, 0.5)

# Used to store the discrepancies between the Magnus matrices and the derivative matrix.  A variety of 
# thresholds are used.  Same naming convention as described above.
discrepancy10 = []
discrepancy1 = []
discrepancye1 = []
discrepancye2 = []
discrepancye3 = []
discrepancye4 = []
discrepancye5 = []
discrepancye6 = []
discrepancye7 = []
discrepancye8 = []
discrepancye9 = []

for a, b, c, d, e, f, g, h, i, j, k,  s in zip (meHs10, meHs1, meHse1, meHse2, meHse3, 
	meHse4, meHse5, meHse6, meHse7, meHse8, meHse9, srgHs):

    discrepancy10.append (diff (a, s, 6))
    discrepancy1.append (diff (b, s, 6))
    discrepancye1.append (diff (c, s, 6))
    discrepancye2.append (diff (d, s, 6))
    discrepancye3.append (diff (e, s, 6))
    discrepancye4.append (diff (f, s, 6))
    discrepancye5.append (diff (g, s, 6))
    discrepancye6.append (diff (h, s, 6))
    discrepancye7.append (diff (i, s, 6))
    discrepancye8.append (diff (j, s, 6))	
    discrepancye9.append (diff (k, s, 6))

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

#################################################
#                              Figure 1
#                       Threshold: 10.0
#################################################
make_single_figure (1, discrepancy10, "Threshold: 10.0", cmap, norm)

#################################################
#                                Figure 2
#                           Threshold: 1.0
#################################################
make_single_figure (2, discrepancy1, "Threshold: 1.0", cmap, norm)

#################################################
#                                Figure 3
#                           Threshold: 0.1
#################################################
make_single_figure (3, discrepancye1, "Threshold: 0.1", cmap, norm)

#################################################
#                                Figure 4
#                           Threshold: 0.01
#################################################
make_single_figure (4, discrepancye2, "Threshold: 0.01", cmap, norm)

#################################################
#                                Figure 5
#                           Threshold: 0.001
#################################################
make_single_figure (5, discrepancye3, "Threshold: 0.001", cmap, norm)

#################################################
#                                         Show
#################################################
# Show the plot
show()
