########################################################################
# VaryDistanceAndInteraction.py
# Julie Butler
# July 12th, 2018
# Version 1.0
#
# Compares the effect of varying the energy level spacing and the interaction strength of the 
# Hamiltonians that result from the Magnus Expansion.  The results from the Mangus Expansion will
# be compared to the results of solving the SRG pairing problems using the derivative scheme.
#########################################################################

#################################################
#
#                                     IMPORTS
#
#################################################
import numpy as np
from numpy import array, arange
from MatrixComparisons import compareSquareMatrices as diff
from srg_pairing import main as srgMain
from MagnusExpansionButler import main as meMain

#################################################
#
#                                  Calculations
#
#################################################
# flow parameters for snapshot images
flowparams = array([0.,0.001,0.01,0.05,0.1, 1., 5., 10.])


# Range of energy level spacings from 0 to 2 by 0.2
d_range = arange (0.0, 2.0, 0.2)

# Range of interaction strengths from -1.0 to 1.0 by 0.2
g_range = arange (-1.0, 1.0, 0.2)

# Fixed values to be used when the other one is varying
# Values taken from srg_pairing.
fixed_d = 1.0
fixed_g = 0.5

# Fixed value of the threshold to be used
threshold_value = 0.001

# To hold the final values
discrepancy_vary_d = []
discrepancy_vary_g = []

# Collects the data for varying the energy spacing
for d in d_range:
	srgHs = srgMain (flowparams, d, fixed_g)
	meHs = meMain (flowparams, threshold_value, d, fixed_g)
	discrepancy = []
	for i, j in zip (meHs, srgHs):
		discrepancy.append (diff (i, j, 6))
	discrepancy_vary_d.append (discrepancy)

# Collects the data for varying the interaction strengrh
for g in g_range:
	srgHs = srgMain (flowparams, fixed_d, g)
	meHs = meMain (flowparams, threshold_value, fixed_d, g)
	discrepancy = []
	for i, j in zip (meHs, srgHs):
		discrepancy.append (diff (i, j, 6))
	discrepancy_vary_g.append (discrepancy)

# Collecting the data for making the energy level spacing graphs
discrepancy_d_0 = []
discrepancy_d_0_001 = []
discrepancy_d_0_01 = []
discrepancy_d_0_1 = []
discrepancy_d_5 = []
for i in discrepancy_vary_d:
	discrepancy_d_0.append (i[0])
	discrepancy_d_0_001.append (i[1])
	discrepancy_d_0_01.append (i[2])
	discrepancy_d_0_1.append (i[4])
	discrepancy_d_5.append (i[5])

# Collecting the data for making the interaction strength graphs
discrepancy_g_0 = []
discrepancy_g_0_001 = []
discrepancy_g_0_01 = []
discrepancy_g_0_1 = []
discrepancy_g_5 = []
for i in discrepancy_vary_g:
	discrepancy_g_0.append (i[0])
	discrepancy_g_0_001.append (i[1])
	discrepancy_g_0_01.append (i[2])
	discrepancy_g_0_1.append (i[4])
	discrepancy_g_5.append (i[5])

#################################################
#
#                              Graphing
#
#################################################
from pylab import *
dot_size = 4

#################################################
#                                   Figure 1
#                Varying the energy level spacing 
#################################################
figure (1)

#################################################
#                                 Subplot 1.1
# Energy level spacing vs discrepancy for energy level spacing values
# between 0.0 and 2.0 with a spacing of 0.2 and flow parameters of 
# 0.0, 0.001, 0.01, 0.1, and 5.0
#################################################
subplot (211)
# Plot the data
plot (d_range, discrepancy_d_0, 'b-o', linewidth=dot_size, label='Flow Parameter: 0.0')
plot (d_range, discrepancy_d_0_001, 'g-o', linewidth=dot_size, label='Flow Parameter: 0.001')
plot (d_range, discrepancy_d_0_01, 'r-o', linewidth=dot_size, label='Flow Parameter: 0.01')
plot (d_range, discrepancy_d_0_1, 'm-o', linewidth=dot_size, label='Flow Parameter: 0.1')
plot (d_range, discrepancy_d_5, 'c-o', linewidth=dot_size, label='Flow Parameter: 5.0')
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
xlabel('Energy Level Spacing', fontsize=16, fontweight='bold')
ylabel('Discrepancy', fontsize=16, fontweight='bold')
# Make the legend, but only for this plot
legend (fontsize = 14)

#################################################
#                                 Subplot 1.2
# Energy level spacing vs discrepancy for energy level spacing values
# between 0.0 and 2.0 with a spacing of 0.2 and flow parameters of 
# 0.0, 0.001, 0.01, 0.1, and 5.0
#################################################
subplot (212)
# Plot the data
plot (d_range, discrepancy_d_0, 'b-o', linewidth=dot_size, label='Flow Parameter: 0.0')
plot (d_range, discrepancy_d_0_001, 'g-o', linewidth=dot_size, label='Flow Parameter: 0.001')
plot (d_range, discrepancy_d_0_01, 'r-o', linewidth=dot_size, label='Flow Parameter: 0.01')
plot (d_range, discrepancy_d_0_1, 'm-o', linewidth=dot_size, label='Flow Parameter: 0.1')
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
xlabel('Energy Level Spacing', fontsize=16, fontweight='bold')
ylabel('Discrepancy', fontsize=16, fontweight='bold')

#################################################
#                                   Figure 2
#                    Varying the interaction strength 
#################################################
figure (2)

#################################################
#                                 Subplot 2.1
# Interaction strength vs discrepancy for energy level spacing values
# between -1.0 and 1.0 with a spacing of 0.2 and flow parameters of 
# 0.0, 0.001, 0.01, 0.1, and 5.0
#################################################
subplot (311)
# Plot the data
plot (g_range, discrepancy_g_0, 'b-o', linewidth=dot_size, label='Flow Parameter: 0.0')
plot (g_range, discrepancy_g_0_001, 'g-o', linewidth=dot_size, label='Flow Parameter: 0.001')
plot (g_range, discrepancy_g_0_01, 'r-o', linewidth=dot_size, label='Flow Parameter: 0.01')
plot (g_range, discrepancy_g_0_1, 'm-o', linewidth=dot_size, label='Flow Parameter: 0.1')
plot (g_range, discrepancy_g_5, 'c-o', linewidth=dot_size, label='Flow Parameter: 5.0')
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
xlabel('Interaction Strength', fontsize=16, fontweight='bold')
ylabel('Discrepancy', fontsize=16, fontweight='bold')
# Make the legend, but only for this plot
legend (fontsize = 14)

#################################################
#                                 Subplot 2.2
# Interaction strength vs discrepancy for energy level spacing values
# between -1.0 and 1.0 with a spacing of 0.2 and flow parameters of 
# 0.0, 0.001, 0.01, 0.1, and 5.0
#################################################
subplot (312)
# Plot the data
plot (g_range, discrepancy_g_0, 'b-o', linewidth=dot_size, label='Flow Parameter: 0.0')
plot (g_range, discrepancy_g_0_001, 'g-o', linewidth=dot_size, label='Flow Parameter: 0.001')
plot (g_range, discrepancy_g_0_01, 'r-o', linewidth=dot_size, label='Flow Parameter: 0.01')
plot (g_range, discrepancy_g_0_1, 'm-o', linewidth=dot_size, label='Flow Parameter: 0.1')
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
xlabel('Interaction Strength', fontsize=16, fontweight='bold')
ylabel('Discrepancy', fontsize=16, fontweight='bold')

#################################################
#                                 Subplot 2.3
# Interaction strength vs discrepancy for energy level spacing values
# between -1.0 and 1.0 with a spacing of 0.2 and flow parameters of 
# 0.0, 0.001, 0.01
#################################################
subplot (313)
# Plot the data
plot (g_range, discrepancy_g_0, 'b-o', linewidth=dot_size, label='Flow Parameter: 0.0')
plot (g_range, discrepancy_g_0_001, 'g-o', linewidth=dot_size, label='Flow Parameter: 0.001')
plot (g_range, discrepancy_g_0_01, 'r-o', linewidth=dot_size, label='Flow Parameter: 0.01')
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
xlabel('Interaction Strength', fontsize=16, fontweight='bold')
ylabel('Discrepancy', fontsize=16, fontweight='bold')

#################################################
#                                        Show
#################################################
# Show the plot
show()
