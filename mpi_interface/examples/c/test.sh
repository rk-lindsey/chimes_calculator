alias gcc=gcc-11 # For Becky's osx situation

make clean-all
make all

ff=../../../serial_interface/tests/force_fields/published_params.HN3.2+3+4b.Tersoff.special.offsets.txt
cf=../../../serial_interface/tests/configurations/HN3.2gcc_3000K.OUTCAR_#000.xyz


mpirun -n 6 ./chimescalc-test_mpi-C $ff $cf 


