##########################################################################################
# Magnus-Butler.py
# Julie Butler
# Version 1.0
# July 9th, 2018
#
# A code that applies the Magnus Expansion to the In-Medium Similarity Renormalization
# Group Theory for nuclear physics.  Some of the code needed to perform calculations will 
# be written in the file linearAlgebra.py, which was created to hold useful snippets of 
# code that oertain to linear algebra calculations.
#
# To-Do
# 1. Document the code
# 2. Run a time test
##########################################################################################

##################################################
#
#                 IMPORTS
#
##################################################
# Methods for performing commutators and nested commutators
from LinearAlgebra import quantumMechanics 
# Allows for symbolic representation of variables
from sympy import *
import numpy as np
from numpy import array, dot, diag, reshape
from scipy.integrate import odeint
from fractions import Fraction as Fr
from sys import exit 
from math import exp
from MatrixComparisons import compareSquareMatrices



##################################################
#
#           GLOBAL VARIABLES
#
##################################################
d = 1
g = 0.5

H0 =  array ([
			 [2*d - g, -g/2, -g/2, -g/2, -g/2, 0],
			 [-g/2, 4*d - g, -g/2, -g/2, 0, -g/2],
			 [-g/2, -g/2, 6*d - g, 0, -g/2, -g/2],
			 [-g/2, -g/2, 0, 6*d - g, -g/2, -g/2],
			 [-g/2, 0, -g/2, -g/2, 8*d - g,  -g/2],
			 [0, -g/2, -g/2, -g/2, -g/2, 10*d - g]
			])
max_n = 1000
#threshold = 0.0000000001

qm = quantumMechanics ()

def bernoulli(n):
	A= [0] * (n+1)
	for m in range(n+1):
		A[m] = Fr(1, m+1)
		for j in range(m, 0, -1):
			A[j-1] = j*(A[j-1] - A[j])
	return A[0] # (which is Bn)

def factorial (n):
	if n == 0:
		return 1
	else:
		return n * factorial (n-1)



def magnusExpansion (omega, s, threshold):
	# need a 6x6 matrix 
	omega = reshape (omega, (6, 6))
	# failsafes so there is no infinite loop
	# Needed for commutator (rewrite this part)
	# Calculates H(s) using the formula H(s) = SUM 1/k! [omega, H(0)]^(k)
	previous_term = (1/factorial (0)) * qm.nestedCommutator (omega, H0, 0)
	Hs = previous_term
	k = 1
	delta = 10000
	while (delta > threshold):
		if (k == max_n):
			print ("Summation failed to converge to required threshold.")
			print ("Program will now terminate")
			exit ()
		current_term = (1/factorial (k)) * qm.nestedCommutator (omega, H0, k)
		Hs = Hs + current_term
		delta = compareSquareMatrices (current_term, previous_term, 6)
		previous_term = current_term
		k = k + 1
	# Calculates eta from H(s) using eta = [Hd, Hod]
	Hd = diag (diag (Hs))
	Hod = Hs - Hd
	eta = qm.commutator(Hd, Hod)
	# Computes the sum domega/dt = SUM Bn/n! [omega, eta]^(n)
	previous_term = (bernoulli (0) / factorial (0)) * qm.nestedCommutator (omega, eta, 0)
	domega_ds = previous_term
	delta = 10000
	n = 1
	while (delta > threshold):
		if (n == max_n):
			print ("Summation failed to converge to required threshold.")
			print ("Program will now terminate")
			exit ()
		current_term = (bernoulli (n) / factorial (n)) * qm.nestedCommutator (omega, eta, n)
		domega_ds = domega_ds + current_term
		delta = compareSquareMatrices (current_term, previous_term, 6)
		previous_term = current_term 
		n = n + 1
	# needs to be returned as a single dimension array for ODE solver
	domega_ds = reshape (domega_ds, -1)
	return domega_ds.tolist()
		




def main (flowparams, threshold):
	omega0 = array ([
					[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
				])
	omega0 = reshape (omega0, -1)
	#flowparams = array([0.,0.001,0.01,0.05,0.1, 1., 5., 10.])

	omegas = odeint (magnusExpansion, omega0, flowparams, args=(threshold,))
	Hs_list = []
	for omega, s in zip(omegas, flowparams):
		omega = reshape (omega, (6, 6))
		# Needed for commutator (rewrite this part)
		# Calculates H(s) using the formula H(s) = SUM 1/k! [omega, H(0)]^(k)
		previous_term = (1/factorial (0)) * qm.nestedCommutator (omega, H0, 0)
		Hs = previous_term
		k = 1
		delta = 10000
		while (delta > threshold):
			if (k == max_n):
				print ("Summation failed to converge to required threshold.")
				print ("Program will now terminate")
				exit ()
			current_term = (1/factorial (k)) * qm.nestedCommutator (omega, H0, k)
			Hs = Hs + current_term
			delta = compareSquareMatrices (current_term, previous_term, 6)
			previous_term = current_term
			k = k + 1
		Hs_list.append(Hs)
	return Hs_list
	





