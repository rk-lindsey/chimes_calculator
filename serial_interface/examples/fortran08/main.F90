! ChIMES Calculator
! Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
! Contributing Author:  Nir Goldman (2020)
! Contributing Author: Bálint Aradi (2021)

program test_F08_api
  use chimes_serial08, only : ChimesCalc, ChimesCalc_init
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: GPa = 6.9479 ! convert kcal/mol.A^3 to GPa

  type(ChimesCalc) :: chimes
  integer :: io_num, fd
  character(1024) :: coord_file, param_file
  integer :: iatom, natom
  real(dp) :: latvecs(3, 3), stress(3, 3)
  real(dp) :: energy
  real(dp), allocatable :: coords(:,:), forces(:,:)
  character(len=10), allocatable :: atom_types(:)


  io_num = command_argument_count()
  if (io_num < 2) then
    print *, "To run: ./test_F08.x <parameter file> <xyz config. file>"
    print *, "Exiting code."
    error stop
  end if

  call get_command_argument(1, param_file)
  call get_command_argument(2, coord_file)

  open (newunit=fd, status='old', file=coord_file)
  read(fd, *) natom
  read(fd, *) latvecs(:, 1), latvecs(:, 2), latvecs(:, 3)

  allocate(forces(3, natom))
  allocate(atom_types(natom))
  allocate(coords(3, natom))

  read(fd, *) (atom_types(iatom), coords(:, iatom), iatom = 1, natom)
  close(fd)

  energy = 0.0_dp
  forces(:,:) = 0.0_dp
  stress(:,:) = 0.0_dp

  ! Initialize ChiMES calculator
  call ChimesCalc_init(chimes, trim(param_file), 0)

  ! Set that atom types. You need this call only once (unless atom types change)
  call chimes%set_atom_types(atom_types)

  ! Get contributions. This call can be done several times with various coords/latvecs.
  call chimes%calculate(coords, latvecs, energy, forces, stress)

  print *, "Success!"
  print '(A,1X, F0.6)', "Energy (kcal/mol)", energy
  print *, "Stress tensors (GPa)"
  print '(A,1X, F15.6)', "s_xx: ", stress(1, 1) * GPa
  print '(A,1X, F15.6)', "s_yy: ", stress(2, 2) * GPa
  print '(A,1X, F15.6)', "s_zz: ", stress(3, 3) * GPa
  print '(A,1X, F15.6)', "s_xy: ", stress(1, 2) * GPa
  print '(A,1X, F15.6)', "s_xz: ", stress(1, 3) * GPa
  print '(A,1X, F15.6)', "s_yz: ", stress(2, 3) * GPa
  print *, "Forces (kcal/mol)"
  print '(F15.6)', forces

#if DEBUG==1

  open (newunit=fd, status='replace', file='debug.dat')
  write(fd, '(F15.6)') energy
  write(fd, '(F15.6)') stress(1, 1) * GPa
  write(fd, '(F15.6)') stress(2, 2) * GPa
  write(fd, '(F15.6)') stress(3, 3) * GPa
  write(fd, '(F15.6)') stress(1, 2) * GPa
  write(fd, '(F15.6)') stress(1, 3) * GPa
  write(fd, '(F15.6)') stress(2, 3) * GPa
  write(fd, '(E15.6)') forces
  close(fd)

#endif

end program
