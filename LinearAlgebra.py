##########################################################################################
# LinearAlgebra.py
# Julie Butler
# Version 1.0
# July 9th, 2018
#
# A colection of methods that perform linear algebra related calculations.  The code is
# split into different classes by topic/subject.
##########################################################################################

##########################################################################################
# CLASSES:
# quantumMechanics: methods that are especially perternate to quantum mechanics 
#	calculations 
# linearAlgebra: basic linear algebra calculations that are not included in numpy or that
#	needed to be modified from the numpy versions
##########################################################################################

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
#              quantumMechanics
#
# Methods:
# commutator: takes in two numpy arrays and returns the commutator of them, which is
#	defined as AB-BA.
# nestedCommutator: takes in two numpy arrays and an integer and performs the nested
#	commutator [a, b]^(n), where [a, b]^(n) = [a, [a, b]^(n-1)] and 
#	[a, b]^(0) = b
##################################################

class quantumMechanics:
	##################################################
	#
	#                   __init__
	#
	##################################################
	def __init__ (self):
		"""
		Input:
			None.
		Output:
			None.
			
		Initializes an object of the class, but does nothing else.
		"""
	
	##################################################
	#
	#                  commutator
	#
	##################################################
	def commutator (self, x, y):
		return dot (x, y) - dot (y, x)
	
	##################################################
	#
	#              nestedCommutator
	#
	##################################################
	def nestedCommutator (self, x, y, n):
		"""
		Input:
			n (an integer): the number of times the nested commutator routine
				is to be performed.
		Output:
			Unnamed (a double): The result of performing the nested commutator n
				times on class level numpy arrays a and b.
				
		Performs a nested commutator with n iterations on teh class level numpy arrays
			a and b.  For the nested commutator [a, b]^(n) = [a, [a, b]^(n-1)] and
			[a, b]^(0) = b.
		"""
		if n == 0:
			return y
		else:
			return self.commutator (x, self.nestedCommutator (x, y, n-1))
		
