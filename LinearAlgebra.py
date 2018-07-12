#########################################################################
# LinearAlgebra.py
# Julie Butler
# Version 1.1
# July 9th, 2018
#
# A colection of methods that perform linear algebra related calculations.  The methods are organized
# alphabetically, with an outline of all methods given below.
#########################################################################

#########################################################################
# METHODS:
# commutator (x, y): computes the commutator of the numpy arrays x and y.
# nestedCommutator (x, y, n): computes the nested commutator to the nth degree of the numpy 
#	arrays x and y
#########################################################################
##################################################
#
#                 IMPORTS
#
##################################################
# Allows for the use of effecient arrays and must basic linear algebra calculations
import numpy as np 
# Commonly used methods from the numpy library 
from numpy import dot

##################################################
#
#                  commutator
#
##################################################
def commutator (x, y):
	"""	
	Input:
		x (a numpy array): The first array that the nested commutator is to be performed on.
		y (a numpy array): Thesecond array that the nested commutator is to be performed on.
	Output:
		Unnamed (depends on input):  The results of performing the commutator on x and y

	 Performs a commutator on the numpy arrays x and y.  The commutator of A and B is defined 
	as [A, B] = AB - BA.  The order of the inputs could matter, depending on the whether the inputs
	are matrices or just vectors.  
	"""
	return dot (x, y) - dot (y, x)

##################################################
#
#              nestedCommutator
#
##################################################
def nestedCommutator (x, y, n):
	"""
	Input:
		x (a numpy array):  The first array that the nested commutator is to be performed
			on.
		y (a numpy array): The second array that the nested commutator is to be performed
			on.
		n (an integer): the number of times the nested commutator routine
			is to be performed.
	Output:
		Unnamed (depends on input): The result of performing the nested commutator n
			times on class level numpy arrays a and b.
			
	Performs a nested commutator with n iterations on the numpy arrays x and y.  For the 
	nested commutator [a, b]^(n) = [a, [a, b]^(n-1)] and [a, b]^(0) = b.  The order of x and y
	is important as it will affect the results.
	"""
	if n == 0:
		return y
	else:
		return commutator (x, nestedCommutator (x, y, n-1))
		
