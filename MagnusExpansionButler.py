#########################################################################
# Magnus-Butler.py
# Julie Butler
# Version 1.0
# July 9th, 2018
#
# A code that applies the Magnus Expansion to the In-Medium Similarity Renormalization
# Group Theory for nuclear physics. The specific case this Magnus Expansion is written for is the 
# pairing model containing eight single particle states and four particles.
#
# To-Do
# 1. Document the code
#########################################################################

#########################################################################
# The Magnus Expansion and its Algorithm:
#
#
#########################################################################


##################################################
#
#                 IMPORTS
#
##################################################
# Methods for performing commutators and nested commutators
from LinearAlgebra import commutator, nestedCommutator 
# Allows for matrices to be compared in a way that will satisfy loop conditions.
# See below.
from MatrixComparisons import compareSquareMatrices
# Numpy imports which are needed for creating and manipulating the matrices
import numpy as np
from numpy import array, dot, diag, reshape
# The ODE solver, see below
from scipy.integrate import odeint
# Needed to calculated the Bernoulli numbers
from fractions import Fraction as Fr
# Allows for the program to terminate itself incase one of the loops does not converge
from sys import exit 

##################################################
# compareSquareMatrices
#
#
##################################################

##################################################
# odeint:
#
#
##################################################

##################################################
#
#                      GLOBAL VARIABLES
#
##################################################
# Sets the maximum number of iterations that a loop can perform.  If the number of iterations reaches
# the value of man_n, the program terminates with an error message.
max_iterations = 100


##################################################
#
#                                 bernoulli
#
##################################################
def bernoulli(n):
	"""
	Input:
		n (an integer):  The index of the Bernoulli number to be calculated
	Output:
		A[0] (a fraction): The nth Bernoulli number
	Algorithm taken from the following website, unmodified:
	https://rosettacode.org/wiki/Bernoulli_numbers#Python
	"""
	A= [0] * (n+1)
	for m in range(n+1):
		A[m] = Fr(1, m+1)
		for j in range(m, 0, -1):
			A[j-1] = j*(A[j-1] - A[j])
	return A[0] # (which is Bn)

##################################################
#
#                                 factorial 
#
##################################################
def factorial (n):
	"""
	Input:
		n (an integer): The number which is to be used as the starting point for 
			calculating the factorial.
	Output:
		Unnamed (an integer): The nth factorial.
	Calculates the factorial of n.  The factorial of n is denoted as n! and is defined as
	n! = n*(n-1)*(n-2)*...*2*1, where 0! = 1.  The algorithm is implemented using 
	recursion with the base case being zero.
	"""
	if n == 0:
		return 1
	else:
		return n * factorial (n-1)

##################################################
#
#                                 hamiltonian
#
##################################################
def hamiltonian (d, g):
	"""
	Input:
		d (a double): The energy  level spacing 
		g (a double): The interaction strength
	Output:
		Unnamed (a 6x6 numpy matrix): The result of calculating the Hamiltonian matrix, 
			H0, with the specified values of d and g
	Calculates the 6x6 Hamiltonian matrix, H0, which describes the pairing model with eight states
	and four particles using the specified values of d and g.
	"""
	return array([
			[2*d-g, -g/2, -g/2, -g/2, -g/2, 0],
			[-g/2, 4*d-g, -g/2, -g/2, 0, -g/2],
			[-g/2, -g/2, 6*d-g, 0, -g/2, -g/2],
			[-g/2, -g/2, 0, 6*d-g, -g/2, -g/2],
			[-g/2, 0, -g/2, -g/2, 8*d-g, -g/2],
			[0, -g/2, -g/2, -g/2, -g/2, 10*d-g]
		])


##################################################
#
#                              mangusExpansion
#
##################################################
def magnusExpansion (omega, s, threshold, H0):
	# Reshapes the inputed omega array from a single dimension to a 6x6 matrix
	omega = reshape (omega, (6, 6))

	# Calculates H(s) using the formula H(s) = SUM 1/k! [omega, H(0)]^(k)
	previous_term = (1/factorial (0)) * nestedCommutator (omega, H0, 0)
	Hs = previous_term
	k = 1
	delta = 10000
	while (delta > threshold):
		if (k == max_iterations):
			print ("Summation failed to converge to required threshold.")
			print ("Program will now terminate")
			exit ()
		current_term = (1/factorial (k)) * nestedCommutator (omega, H0, k)
		Hs = Hs + current_term
		delta = compareSquareMatrices (current_term, previous_term, 6)
		previous_term = current_term
		k = k + 1

	# Calculates eta from H(s) using eta = [Hd, Hod]
	Hd = diag (diag (Hs))
	Hod = Hs - Hd
	eta = commutator(Hd, Hod)

	# Computes the sum domega/dt = SUM Bn/n! [omega, eta]^(n)
	previous_term = (bernoulli (0) / factorial (0)) * nestedCommutator (omega, eta, 0)
	domega_ds = previous_term
	delta = 10000
	n = 1
	while (delta > threshold):
		if (n == max_iterations):
			print ("Summation failed to converge to required threshold.")
			print ("Program will now terminate")
			exit ()
		current_term = (bernoulli (n) / factorial (n)) * nestedCommutator (omega, eta, n)
		domega_ds = domega_ds + current_term
		delta = compareSquareMatrices (current_term, previous_term, 6)
		previous_term = current_term 
		n = n + 1

	# needs to be returned as a single dimension array for ODE solver
	domega_ds = reshape (domega_ds, -1)
	return domega_ds.tolist()
		




def main (flowparams, threshold, d, g):
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
	H0 = hamiltonian (d, g)
	omegas = odeint (magnusExpansion, omega0, flowparams, args=(threshold,H0,))
	Hs_list = []
	for omega, s in zip(omegas, flowparams):
		omega = reshape (omega, (6, 6))
		# Needed for commutator (rewrite this part)
		# Calculates H(s) using the formula H(s) = SUM 1/k! [omega, H(0)]^(k)
		previous_term = (1/factorial (0)) * nestedCommutator (omega, H0, 0)
		Hs = previous_term
		k = 1
		delta = 10000
		while (delta > threshold):
			if (k == max_iterations):
				print ("Summation failed to converge to required threshold.")
				print ("Program will now terminate")
				exit ()
			current_term = (1/factorial (k)) * nestedCommutator (omega, H0, k)
			Hs = Hs + current_term
			delta = compareSquareMatrices (current_term, previous_term, 6)
			previous_term = current_term
			k = k + 1	
		Hs_list.append(Hs)
	return Hs_list
	





