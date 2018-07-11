def compareSquareMatrices (x, y, dim):
	total_difference = 0
	for i in range (0, dim):
		for j in range (0, dim):
			total_difference = total_difference + abs(abs(x[i][j]) - abs(y[i][j]))
	return abs(total_difference)
