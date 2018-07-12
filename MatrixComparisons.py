#########################################################################
# Matrix Comparisons.py
# Julie Butler
# July 12th, 2018
# Version 1.0
#
# A collection of methods that allow for the values contained in two matrices to be eassily compared.
# Currently only implemented for square matrices.
#
# TO-DO:
# 1. Document the code
#########################################################################

#########################################################################
# Methods:
# compareSquareMatrices (x, y, dim): Takes the two square matrices x and y, each of dimension dim,
#	and computes the equation delta = SUM_{i = 0}^d SUM_{j=0}^d ||xij|-|yij||, to give a 
#	numerival answer to determine the difference between the two matrices.
# compareSquareMatricesByElement (x, y, dim): Takes the two square matrices x and y, each of 
#	dimension dim, and computes the difference by matrix element and returns a square matrix of 
#	dimension dim.  For example, if the returned matrix is called delta, 
#	then delta_{ii} = ||x_{ii}|-|y_{ii}||
#########################################################################

#################################################
#
#                      compareSquareMatrices
#
#################################################
def compareSquareMatrices (x, y, dim):
	"""
	Input:
		x (a matrix): a square matrix of dimension dim
		y (a matrix): a square matrix of dimension dim.  Must be the same size as x.
		dim (an integer):  The size of one side of one of the matrices.
	Output:
		total_difference (a double): The result of performing the calculation
			SUM_{i = 0}^d SUM_{j=0}^d ||xij|-|yij||
	Computes the difference between two matrices and converting this difference into a single 
	number using the equation delta = SUM_{i = 0}^d SUM_{j=0}^d ||xij|-|yij||.  This method
	is useful when making comparisons for loops and when graphing.
	"""
	total_difference = 0
	for i in range (0, dim):
		for j in range (0, dim):
			total_difference = total_difference + abs(abs(x[i][j]) - abs(y[i][j]))
	return total_difference

#################################################
#
#               compareSquareMatricesByElement
#
#################################################
def compareSquareMatricesByElement (x, y, dim):
	"""
	Input:
		x (a matrix): a square matrix of dimension dim
		y (a matrix): a square matrix of dimension dim.  Must be the same size as x.
		dim (an integer):  The size of one side of one of the matrices.
	Output:
		total_difference (a matrix): Each element of total_difference, delta_{ij} is calcualted as 
			follows: delta_{ij} = ||x_{ij}| - |y_{ij}||.
	Calcualtes the difference between two square matrices by element.   Useful when looking into 	
	patterns between two sets of matrix data.
	"""
	total_difference = []
	for i in range (0, dim):
		row_difference = []
		for j in range (0, dim):
			row_difference.append (abs (abs (x[i][j]) - abs (y[i][j])))
		total_difference.append (row_difference)
	return total_difference
				
