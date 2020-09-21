import sys

infile = open(sys.argv[1],'r')

with open(sys.argv[1]) as infile:

	line = infile.readline().split()
	
	while line:

		min_prec = min([line[0][::-1].find('.'),line[1][::-1].find('.')])
		
		val_1 = round(float(line[0]), min_prec)
		val_2 = round(float(line[1]), min_prec)
		
		if abs( val_1 - val_2) >= 1E-3:
			print("Values do not match:")
			print("min prec", min_prec)
			print( "vals:", val_1, val_2)
			print("	",line[0], line[1])
			print("")

		line = infile.readline().split()
		
