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
########################################################################
# The first two imports run the two different methods for solving the SRG problem 
from MagnusExpansionButler import main as meMain
from srg_pairing import main as srgMain
# A useful method for finding the difference between two square matrices (see below)
# Maybe moved to a more general math file later
from MagnusExpansionButler import compareSquareMatrices as diff 
# Two needed methods from the numpy collection
from numpy import array, arange

#################################################
# MagnusExpansionButler.compareSquareMatrices
# A method of comparing the differences between two square matrices
# but it returns a numbers, which is good from loop comparisons and 
# graphing.  For the two square matrices x and y of dimension d, the 
# difference between the two matrices, delta, is defined as 
# delta = SUM_{i = 0}^d SUM_{j=0}^d ||xij|-|yij||
#################################################

# Gets a list of numbers to be used for flow parameters to solve the differential equation.  current
# values range from 0 to 10 with a spacing of 0.05, leading to 200 total flow parameters used to 
# solve the differential equations
flowparams_new = arange (0, 10, 0.05)

# Getting the H(s) matrices using the Magnus expansion for a variety of threshold values, given in 
# the parameters.
meHs = meMain (flowparams_new, 10)
meHs0 = meMain (flowparams_new, 1)
meHs1 = meMain (flowparams_new, 0.1)
meHs2 = meMain (flowparams_new, 0.01)
meHs3 = meMain (flowparams_new, 0.001)
meHs4 = meMain (flowparams_new, 0.0001)
meHs5 = meMain (flowparams_new, 0.00001)
meHs6 = meMain (flowparams_new, 0.000001)
meHs7 = meMain (flowparams_new, 0.0000001)
meHs8 = meMain (flowparams_new, 0.00000001)
meHs9 = meMain (flowparams_new, 0.000000001)

# Gets the H(s) matrices using the derivative method 
srgHs = srgMain (flowparams_new)

# Are g
discrepancy = []
discrepancy0 = []
discrepancy1 = []
discrepancy2 = []
discrepancy3 = []
discrepancy4 = []
discrepancy5 = []
discrepancy6 = []
discrepancy7 = []
discrepancy8 = []
discrepancy9 = []

for a, b, c, d, e, f, g, h, i, j, k,  s in zip (meHs, meHs0, meHs1, meHs2, meHs3, meHs4, meHs5, meHs6, meHs7, meHs8, meHs9, srgHs):
	discrepancy.append (diff (a, s, 6))
	discrepancy0.append (diff (b, s, 6))
	discrepancy1.append (diff (c, s, 6))
	discrepancy2.append (diff (d, s, 6))
	discrepancy3.append (diff (e, s, 6))
	discrepancy4.append (diff (f, s, 6))
	discrepancy5.append (diff (g, s, 6))
	discrepancy6.append (diff (h, s, 6))
	discrepancy7.append (diff (i, s, 6))
	discrepancy8.append (diff (j, s, 6))	
	discrepancy9.append (diff (k, s, 6))



from pylab import *
rc('axes', linewidth=2)
figure (1)
subplot (321)
# Make a dummy plot
plot (flowparams_new, discrepancy,'bo',  linewidth = 1, label='Threshold: 10')
plot (flowparams_new, discrepancy1, 'go',linewidth = 1, label='Threshold: 0.1')
plot (flowparams_new, discrepancy3, 'ro',linewidth = 1, label='Threshold: 0.001')
plot (flowparams_new, discrepancy5, 'co',linewidth = 1, label='Threshold: 0.00001')
plot (flowparams_new, discrepancy7, 'mo',linewidth = 1, label='Threshold: 0.0000001')
# Change size and font of tick labels
# Again, this doesn't work in interactive mode.
fontsize = 12
ax = gca()

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

xlabel('Flow Parameter', fontsize=16, fontweight='bold')
ylabel('Discrepancy', fontsize=16, fontweight='bold')

legend (fontsize = 14)
subplot (323)
plot (flowparams_new[79:], discrepancy1[79:], 'go',linewidth = 1, label='Threshold: 0.1')
plot (flowparams_new[79:], discrepancy3[79:], 'ro',linewidth = 1, label='Threshold: 0.001')
plot (flowparams_new[79:], discrepancy5[79:], 'co',linewidth = 1, label='Threshold: 0.00001')
plot (flowparams_new[79:], discrepancy7[79:], 'mo',linewidth = 1, label='Threshold: 0.0000001')

# Change size and font of tick labels
# Again, this doesn't work in interactive mode.
fontsize = 12
ax = gca()

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

xlabel('Flow Parameter', fontsize=16, fontweight='bold')
ylabel('Discrepancy', fontsize=16, fontweight='bold')

#legend (fontsize = 14)

subplot (325)
plot (flowparams_new[79:], discrepancy3[79:], 'ro',linewidth = 1, label='Threshold: 0.001')
plot (flowparams_new[79:], discrepancy5[79:], 'co',linewidth = 1, label='Threshold: 0.00001')
plot (flowparams_new[79:], discrepancy7[79:], 'mo',linewidth = 1, label='Threshold: 0.0000001')

# Change size and font of tick labels
# Again, this doesn't work in interactive mode.
fontsize = 12
ax = gca()

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

xlabel('Flow Parameter', fontsize=16, fontweight='bold')
ylabel('Discrepancy', fontsize=16, fontweight='bold')

#legend (fontsize = 14)

# Save figure
savefig('magnus_derivative_comparison.png')

show()
