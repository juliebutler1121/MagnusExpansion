########################################################################
# VaryThreshold.py
# Julie Butler
# July 11th, 2018
# Version 1.0
#
# A supporting code which runs the SRG code from srg_pairing.py and MagnusExpansionButler.py
# Calculates the matrix H(s) at a variety of different flow parameters.  For the magnus expansion, it
# also calculates the matrices using a variety of threshold values for loop termination.
#
# To-do:
# 1. Document the code
# 2. Possibly change the cosmetics of the graph (B&W vs color)
# 3. change the naming scheme of meHs and discrepancy
# 4. Make three new graphs threshold vs discrepancy for 5 values of s
########################################################################
# The first two imports run the two different methods for solving the SRG problem 
from MagnusExpansionButler import main as meMain
from srg_pairing import main as srgMain
# A useful method for finding the difference between two square matrices (see below)
from MatrixComparisons import compareSquareMatricesByElement as diff
# Two needed methods from the numpy collection
from numpy import array, arange
# Nice printing method for debugging purposes
from HelpfulPrintingJunk import print2DMatricesByRow as printer

#################################################
# MagnusExpansionButler.compareSquareMatrices
# A method of comparing the differences between two square matrices
# but it returns a numbers, which is good from loop comparisons and 
# graphing.  For the two square matrices x and y of dimension d, the 
# difference between the two matrices, delta, is defined as 
# delta = SUM_{i = 0}^d SUM_{j=0}^d ||xij|-|yij||
#################################################

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
meHs10 = meMain (flowparams_new, 10, 1, 0.5)
meHs1 = meMain (flowparams_new, 1, 1, 0.5)
meHse1 = meMain (flowparams_new, 0.1, 1, 0.5)
meHse2 = meMain (flowparams_new, 0.01, 1, 0.5)
meHse3 = meMain (flowparams_new, 0.001, 1, 0.5)
meHse4 = meMain (flowparams_new, 0.0001, 1, 0.5)
meHse5 = meMain (flowparams_new, 0.00001, 1, 0.5)
meHse6 = meMain (flowparams_new, 0.000001, 1, 0.5)
meHse7 = meMain (flowparams_new, 0.0000001, 1, 0.5)
meHse8 = meMain (flowparams_new, 0.00000001, 1, 0.5)
meHse9 = meMain (flowparams_new, 0.000000001, 1, 0.5)

# Gets the H(s) matrices using the derivative method 
srgHs = srgMain (flowparams_new, 1, 0.5)

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
from pylab import *
import matplotlib as mpl

#################################################
#                                Figure 1
#                           Threshold: 10
#################################################
figure (1)
# Set the figure title
suptitle ('Threshold: 10.0', fontsize=40, fontweight='bold')

#################################################
#                                   Subplot 1.1
#                           Flow Parameter: 0.05
#################################################
subplot (221)
# Set the subplot title
title ('Flow Parameter: 0.05', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancy10 [1]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Create discrete colormap
cmap = mpl.colors.ListedColormap(['0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7',
							'0.8', '0.9', '1.0'])
bounds = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
fontsize = 12
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 1.2
#                           Flow Parameter: 0.3
#################################################
subplot (222)
# Set the subplot title
title ('Flow Parameter: 0.3', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancy10 [6]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 1.3
#                           Flow Parameter: 1.0
#################################################
subplot (223)
# Set the subplot title
title ('Flow Parameter: 1.0', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancy10 [100]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 1.4
#                           Flow Parameter: 5.0
#################################################
subplot (224)
# Set the subplot title
title ('Flow Parameter: 5.0', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancy10 [100]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                Figure 2
#                           Threshold: 1.0
#################################################
figure (2)
# Set the figure title
suptitle ('Threshold: 1.0', fontsize=40, fontweight='bold')

#################################################
#                                   Subplot 2.1
#                           Flow Parameter: 0.05
#################################################
subplot (221)
# Set the subplot title
title ('Flow Parameter: 0.05', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancy1 [1]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Create discrete colormap
cmap = mpl.colors.ListedColormap(['0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7',
							'0.8', '0.9', '1.0'])
bounds = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
fontsize = 12
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 2.2
#                           Flow Parameter: 0.3
#################################################
subplot (222)
# Set the subplot title
title ('Flow Parameter: 0.3', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancy1 [6]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 2.3
#                           Flow Parameter: 1.0
#################################################
subplot (223)
# Set the subplot title
title ('Flow Parameter: 1.0', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancy1 [100]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 2.4
#                           Flow Parameter: 5.0
#################################################
subplot (224)
# Set the subplot title
title ('Flow Parameter: 5.0', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancy1 [100]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                Figure 3
#                           Threshold: 10
#################################################
figure (3)
# Set the figure title
suptitle ('Threshold: 0.1', fontsize=40, fontweight='bold')

#################################################
#                                   Subplot 3.1
#                           Flow Parameter: 0.05
#################################################
subplot (221)
# Set the subplot title
title ('Flow Parameter: 0.05', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancye1 [1]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
fontsize = 12
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 3.2
#                           Flow Parameter: 0.3
#################################################
subplot (222)
# Set the subplot title
title ('Flow Parameter: 0.3', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancye1 [6]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 3.3
#                           Flow Parameter: 1.0
#################################################
subplot (223)
# Set the subplot title
title ('Flow Parameter: 1.0', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancye1 [100]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 3.4
#                           Flow Parameter: 5.0
#################################################
subplot (224)
# Set the subplot title
title ('Flow Parameter: 5.0', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancye1 [100]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                Figure 4
#                           Threshold: 10
#################################################
figure (4)
# Set the figure title
suptitle ('Threshold: 0.01', fontsize=40, fontweight='bold')

#################################################
#                                   Subplot 4.1
#                           Flow Parameter: 0.05
#################################################
subplot (221)
# Set the subplot title
title ('Flow Parameter: 0.05', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancye2 [1]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
fontsize = 12
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 4.2
#                           Flow Parameter: 0.3
#################################################
subplot (222)
# Set the subplot title
title ('Flow Parameter: 0.3', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancye2 [6]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 4.3
#                           Flow Parameter: 1.0
#################################################
subplot (223)
# Set the subplot title
title ('Flow Parameter: 1.0', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancye2 [100]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 4.4
#                           Flow Parameter: 5.0
#################################################
subplot (224)
# Set the subplot title
title ('Flow Parameter: 5.0', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancye2 [100]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                Figure 5
#                           Threshold: 10
#################################################
figure (5)
# Set the figure title
suptitle ('Threshold: 0.001', fontsize=40, fontweight='bold')

#################################################
#                                   Subplot 5.1
#                           Flow Parameter: 0.05
#################################################
subplot (221)
# Set the subplot title
title ('Flow Parameter: 0.05', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancye3 [1]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
fontsize = 12
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 5.2
#                           Flow Parameter: 0.3
#################################################
subplot (222)
# Set the subplot title
title ('Flow Parameter: 0.3', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancye3 [6]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 5.3
#                           Flow Parameter: 1.0
#################################################
subplot (223)
# Set the subplot title
title ('Flow Parameter: 1.0', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancye3 [100]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

#################################################
#                                   Subplot 5.4
#                           Flow Parameter: 5.0
#################################################
subplot (224)
# Set the subplot title
title ('Flow Parameter: 5.0', fontsize=20, fontweight='bold')
# Collect the data
matrix = discrepancye3 [100]
# Find the maximum value
local_maxes = [max(sub_array) for sub_array in matrix]
global_max = max (local_maxes)
# Convert the elements of the matrix to ratios
for i in range (0, 6):
	for j in range (0, 6):
		matrix[i][j] = matrix[i][j] / global_max
# Plot the matrix using the created color map
imshow(matrix, interpolation='nearest', 
                 extent=[0.5, 0.5+6, 0.5, 0.5+6],
                 cmap=cmap,
		 norm=norm)
# Create the color bar legend
colorbar()
# Change size and font of tick labels
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')


#################################################
#                                         Show
#################################################
# Show the plot
show()
