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
########################################################################
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
