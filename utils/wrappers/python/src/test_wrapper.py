import wrapper_python

#param_file = "/Users/Lindsey11/My_codes/chimes_calculator/utils/wrappers/cpp_and_c/tests/liquid_carbon/params.2b_liqcarb.txt"
param_file = "/Users/Lindsey11/My_codes/chimes_calculator/utils/wrappers/cpp_and_c/tests/uranium_nitride/params_UN_2B-12_3B-8_4B-4"


wrapper_python.set_chimes()
wrapper_python.read_params(param_file)

rij = 2.26
dr  = [1.0, 2.5, 0.6]
typ = ["N", "U"]
frc = [None]*3
frc   [0] = [1.0,2.0,3.0]
frc   [1] = [4.0,5.0,6.0]
prs = [0.0]*9
enr = 300.0

for i in range(9):
	prs[i] = i*i

wrapper_python.chimes_compute_2b_props(rij, dr, typ, frc, prs, enr)


rij = [2.26, 2.26, 3.25]
vec = [[3.0, 3.0, 3.0],[3.0, 3.0, 3.0],[3.0, 3.0, 3.0]]
typ = ["N", "U", "U"]
frc = [None]*3
frc   [0] = [1.0,2.0,3.0]
frc   [1] = [4.0,5.0,6.0]
frc   [2] = [7.0,8.0,9.0]
prs = [2.0]*9
enr = 150.0

if wrapper_python.get_chimes_3b_order() > 0:
	wrapper_python.chimes_compute_3b_props(rij, vec, typ, frc, prs, enr)




rij = [2.26, 2.26, 3.25, 2.26, 2.26, 3.25]
vec = [[3.0, 3.0, 3.0],[3.0, 3.0, 3.0],[3.0, 3.0, 3.0],[3.0, 3.0, 3.0],[3.0, 3.0, 3.0],[3.0, 3.0, 3.0]]
typ = ["N", "U", "U", "N"]
frc = [None]*4
frc   [0] = [1.0,2.0,3.0]
frc   [1] = [4.0,5.0,6.0]
frc   [2] = [7.0,8.0,9.0]
frc   [3] = [10.0,11.0,12.0]
prs = [2.0]*9
enr = 150.0

if wrapper_python.get_chimes_4b_order() > 0:
	wrapper_python.chimes_compute_4b_props(rij, vec, typ, frc, prs, enr)


