########################################################################
# VaryThreshold.py
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
########################################################################

#################################################
#
#                      Imports
#
#################################################
# Third-Party Imports
from numpy import array, arange
# Local Imports
from MagnusExpansionButler import main as magnus_expansion_main
from srg_pairing import main as srg_main
from MatrixComparisons import compare_square_matrices as diff
from GraphingCapabilities import set_up_graph

#################################################
#
#                        Discrepancy calculations
#
#################################################
# Gets a list of numbers to be used for flow parameters to solve the differential equation. 
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

# Used to store the discrepancies between the Magnus matrices and the derivative matrix.  
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
#                             Graphing
#
#################################################
# Import here to prevent an error
from pylab import *
rc('axes', linewidth=2)

#################################################
#
#                                       Figure 1
#
#################################################
figure (1)

#################################################
#                                      Subplot 1.1
# Flow parameter vs. Discrepancy for threshold values of 10, 10e-1, 
# 10e-3, 10e-5, 10e-7 and flow paramter values of 0 - 10.
#################################################
subplot (311)
# Plot the data
plot (flowparams_new, discrepancy10,'bo',  linewidth = 1, label='Threshold: 10')
plot (flowparams_new, discrepancye1, 'go',linewidth = 1, label='Threshold: 0.1')
plot (flowparams_new, discrepancye3, 'ro',linewidth = 1, label='Threshold: 0.001')
plot (flowparams_new, discrepancye5, 'co',linewidth = 1, label='Threshold: 0.00001')
plot (flowparams_new, discrepancye7, 'mo',linewidth = 1, label='Threshold: 0.0000001')
set_up_graph ("Flow Parameter", "Discrepancy", True)

#################################################
#                                 Subplot 1.2
# Flow parameter vs. Discrepancy for threshold values 10e-1, 10e-3,
# 10e-5, and 10e-7 and for flow parameter values 4-10
#################################################
subplot (312)
# Plot the data
plot (flowparams_new[80:], discrepancye1[80:], 'go',linewidth = 1, label='Threshold: 0.1')
plot (flowparams_new[80:], discrepancye3[80:], 'ro',linewidth = 1, label='Threshold: 0.001')
plot (flowparams_new[80:], discrepancye5[80:], 'co',linewidth = 1, label='Threshold: 0.00001')
plot (flowparams_new[80:], discrepancye7[80:], 'mo',linewidth = 1, label='Threshold: 0.0000001')
set_up_graph ("Flow Parameter", "Discrepancy", False)

#################################################
#                                    Subplot 1.3
# Flow parameter vs. Discrepancy for threshold values 10e-3, 10e-5, 
# and 10e-7 and for flow parameters 4 - 10.
#################################################
subplot (313)
# Plot the data
plot (flowparams_new[80:], discrepancye3[80:], 'ro',linewidth = 1, label='Threshold: 0.001')
plot (flowparams_new[80:], discrepancye5[80:], 'co',linewidth = 1, label='Threshold: 0.00001')
plot (flowparams_new[80:], discrepancye7[80:], 'mo',linewidth = 1, label='Threshold: 0.0000001')
set_up_graph ("Flow Parameter", "Discrepancy", False)

#################################################
#
#                                    Figure 2
#
#################################################
figure (2)

#################################################
#                                    Subplot 1.1
# Threshold vs Discrepancy for threshold values 10, 1, 10e-1, 10e-2,
# 10e-3, 10e-4. 1-e05, 10e-6, 10e-7, 10e-8, and 10e-9 and for 
# flow parameter values of 0, 0.3, 1, and 5
#################################################
subplot (311)
# Collect the data
threshold_values = [10, 1, 10e-1, 10e-2, 10e-3, 10e-4, 10e-5, 10e-6, 10e-7, 10e-8, 10e-9]
discrepancy_s_0 = [discrepancy10[0], discrepancy1[0], discrepancye1[0], discrepancye2[0],
        discrepancye3[0], discrepancye4[0], discrepancye5[0], discrepancye6[0], 
        discrepancye7[0], discrepancye8[0], discrepancye9[0]]
discrepancy_s_0_3 = [discrepancy10[6], discrepancy1[6], discrepancye1[6], discrepancye2[6],
        discrepancye3[6], discrepancye4[6], discrepancye5[6], discrepancye6[6], 
        discrepancye7[6], discrepancye8[6], discrepancye9[6]]
discrepancy_s_1 = [discrepancy10[20], discrepancy1[20], discrepancye1[20], discrepancye2[20],
        discrepancye3[20], discrepancye4[20], discrepancye5[20], discrepancye6[20], 
        discrepancye7[20], discrepancye8[20], discrepancye9[20]]
discrepancy_s_5 = [discrepancy10[100], discrepancy1[100], discrepancye1[100], discrepancye2[100],
        discrepancye3[100], discrepancye4[100], discrepancye5[100], discrepancye6[100], 
        discrepancye7[100], discrepancye8[100], discrepancye9[100]]
# Plot the data
plot (threshold_values, discrepancy_s_0, 'b--', linewidth = 4, label='Flow Parameter: 0.0')
plot (threshold_values, discrepancy_s_0_3, 'g--', linewidth = 4, label='Flow Parameter: 0.3')
plot (threshold_values, discrepancy_s_1, 'r--', linewidth = 4, label='Flow Parameter: 1.0')
plot (threshold_values, discrepancy_s_5, 'm--', linewidth = 4, label='Flow Parameter: 5.0')
set_up_graph ("Threshold Value", "Discrepancy", True)

#################################################
#                                    Subplot 2.2
# Threshold vs Discrepancy for threshold values 10e-1, 10e-2,
# 10e-3, 10e-4. 1-e05, 10e-6, 10e-7, 10e-8, and 10e-9 and for 
# flow parameter values of 0, 0.3, 1, and 5
#################################################
subplot (312)
# Plot the data
plot (threshold_values[2:], discrepancy_s_0[2:], 'b--', linewidth = 4, label='Flow Parameter: 0.0')
plot (threshold_values[2:], discrepancy_s_0_3[2:], 'g--', linewidth = 4, label='Flow Parameter: 0.3')
plot (threshold_values[2:], discrepancy_s_1[2:], 'r--', linewidth = 4, label='Flow Parameter: 1.0')
plot (threshold_values[2:], discrepancy_s_5[2:], 'm--', linewidth = 4, label='Flow Parameter: 5.0')
set_up_graph ("Threshold Value", "Discrepancy", False)

#################################################
#                                    Subplot 2.3
# Threshold vs Discrepancy for threshold values 10e-6, 10e-7, 10e-8,
# and 10e-9 and for flow parameter values of 0, 0.3, 1, and 5
#################################################
subplot (313)
# Plot the data
plot (threshold_values[7:], discrepancy_s_0[7:], 'b--', linewidth = 4, label='Flow Parameter: 0.0')
plot (threshold_values[7:], discrepancy_s_0_3[7:], 'g--', linewidth = 4, label='Flow Parameter: 0.3')
plot (threshold_values[7:], discrepancy_s_1[7:], 'r--', linewidth = 4, label='Flow Parameter: 1.0')
plot (threshold_values[7:], discrepancy_s_5[7:], 'm--', linewidth = 4, label='Flow Parameter: 5.0')
set_up_graph ("Threshold Value", "Discrepancy", False)

#################################################
#                                     Show
#################################################
# Show the plot
show()
