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
# 1. Write the code
# 2. Document the code
# 3. Run a time test
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
from scipy.linalg import eigvalsh
from scipy.integrate import odeint

##################################################
#
#           GLOBAL VARIABLES
#
##################################################
d = Symbolic ("d")
g = Symbolic ("g")

H0 = [ [2*d - g, -g/2, -g/2, -g/2, -g/2, 0],
       [-g/2, 4*d - g, -g/2, -g/2, 0, -g/2],
       [-g/2, -g/2, 6*d - g, 0, -g/2, -g/2],
       [-g/2, -g/2, 0, 6*d - g, -g/2, -g/2],
       [-g/2, 0, -g/2, -g/2, 8*d - g, -g/2],
       [0, -g/2, -g/2, -g/2, -g/2, 10*d - g]
     ]

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
		return n * facorial (n-1)

def magnusExpansion ():
	omega0 = 0
	Hd  = diag(diag(H))
	Hod = H-Hd
	qm1 = quamtumMechanics (Hd, Hod)
	eta = qm1.commutator ()
	qm2 = quantumMechanics (omega0, eta)
	previous_term = (bernoulli (0) / factorial (0)) * qm2.nestedCommutator (0)
	threshhold = 0.001
	delta = 1000
	total_sum = previous_term
	n = 1
	while (threshold < delta):
		current_term = (bernoulli (n) / factorial (n)) * qm2.nestedCommutator (n)
		total_sum = total_sum + current_term
		delta = current_term - previous_term
		previous_term = current_term
		n = n + 1
	# turn initial Hamiltonian into a linear array
	y0  = reshape(H0, -1)                 
	# flow parameters for snapshot images
	flowparams = array([0.,0.001,0.01,0.05,0.1, 1., 5., 10.])
	# integrate flow equations - odeint returns an array of solutions,
	# which are 1d arrays themselves
	ys  = odeint(derivative, y0, flowparams, args=(dim,))
	# reshape individual solution vectors into dim x dim Hamiltonian
	# matrices
	Hs  = reshape(ys, (-1, dim,dim))
	



