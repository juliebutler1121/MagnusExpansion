########################################################################
# GraphingCapabilities.py
# Julie Butler
# July 17th, 2018
# Version 1.0
#
# A collection of functions for setting up graphs and formatting them in a way that is useful for 
# insertion into papers.
########################################################################

########################################################################
# Function Outline:
# set_up_graph (x_label, y_label, is_legend): Format the axes of the graph and inserts 
#   legend if wanted.
# make_4_subplots_per_figure (figure_number, data_list, graph_title, subplot_titles, 
#   index_list, cmap, norm): Makes the four subplots for a single threshold value
########################################################################

#################################################
#
#                                 IMPORTS
#
#################################################
# Third-Party Imports
from pylab import *
# Local Imports
from MatrixComparisons import ratio_matrix

def set_up_graph (x_label, y_label, is_legend):
    """
    Format the axes of the graph and inserts legend if wanted.
    Input:
    x_label (a string): The label of the x axis.
    y_label (a string): The label of the y axis.
    is_legend (a boolean): True if a legend is to be added to the graph
    """
    # Change size and font of tick labels
    fontsize = 12
    ax = gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    # Set the axis labels
    xlabel(x_label, fontsize=16, fontweight='bold')
    ylabel(y_label, fontsize=16, fontweight='bold')
    if is_legend:
        # Make the legend
        legend (fontsize = 14)

#################################################
#
#                      make_single_figure
#
#################################################
def make_4_subplots_per_figure (figure_number, data_list, graph_title, 
                subplot_titles, index_list, cmap, norm):
    """
    Makes the four subplots for a single threshold value
    Input:
        figure_number (an integer): the number of the figure
        data_list (an array): Contains the data to be graphed
        graph_title (a string): The title to be printed at the top of the graph
        subplot_titles (an array of strings):  A list of four strings to be used as the
            titles for the four subplots
        index_list (an array of integers): A list of four integers to be used as indices
            to access the data
        cmap and norm (matplotlib objects): handle the color mapping 
    """
    figure (figure_number)
    # Set the figure title
    suptitle (graph_title, fontsize=40, fontweight='bold')
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
            matrix = data_list [index_list [iter_number]]
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
