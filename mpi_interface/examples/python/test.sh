alias gcc=gcc-11 # For Becky's osx situation

# Compile libchimescalc-mpi_dl.so

make clean-all
make all

cp chimescalc_mpi_dl.so libchimescalc-mpi_dl.so

ff=../../../serial_interface/tests/force_fields/published_params.HN3.2+3+4b.Tersoff.special.offsets.txt
cf=../../../serial_interface/tests/configurations/HN3.2gcc_3000K.OUTCAR_#000.xyz
op=0

mpirun -np 4 python main.py $ff $cf $op ../../api 1 > /dev/null
