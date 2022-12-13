! ChIMES Calculator
! Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
! Contributing Author:  Nir Goldman (2020) 

      program test_F_api
      use chimescalc_serial
      use, intrinsic :: ISO_C_binding
      implicit none
      integer io_num, stat, small
      double precision, parameter :: GPa = 6.9479 ! convert kcal/mol.A^3 to GPa
      character(C_char), dimension(1025) :: c_file
      character(1024) :: coord_file, param_file
      CHARACTER ( len = 1024 ) :: wq_char
      integer :: i, j, k, l, natom, ns
      real(C_double) ::   lx, ly, lz
      real(C_double) :: stress(9)
      real(C_double) :: energy
      real(C_double) :: ca(3), cb(3), cc(3)
      real(C_double), allocatable :: xc(:), yc(:), zc(:), fx(:), fy(:), fz(:)
      character(len=10), allocatable, target :: atom_type(:), c_atom(:)
      character(len=2) :: atom
      TYPE(C_PTR), allocatable, dimension(10) :: stringPtr(:)
      integer lenstr

      small = 1

      io_num = command_argument_count()
      if ((io_num .ne. 2) .and. (io_num .ne. 3)) then
        print*,"To run: ./test_F.x <parameter file> <xyz config. file>"
        print*,"or"
        print*,"/test.x <parameter file> <xyz config. file> <allow_replicates(0/1)>"
        print*,"Exiting code.\n"
        STOP
      endif
      
      if (io_num .eq. 3) then
          call GET_COMMAND_ARGUMENT(3, wq_char)
          read(wq_char,*,iostat=stat)  small
      endif
      
      call GET_COMMAND_ARGUMENT(1, wq_char)
      param_file = trim(wq_char)
      call GET_COMMAND_ARGUMENT(2, wq_char)
      coord_file = trim(wq_char)  
      
      print*,"Read args:"      
      do i = 1, io_num
          call GET_COMMAND_ARGUMENT(i, wq_char)
          print*,i,trim(wq_char)
      enddo
      
      open (unit=10, status='old', file=coord_file)
      read(10,*)natom
      read(10,*)ca(1),ca(2),ca(3),cb(1),cb(2),cb(3),cc(1),cc(2),cc(3)

      allocate(fx(natom))
      allocate(fy(natom))
      allocate(fz(natom))
      allocate(atom_type(natom))
      allocate(c_atom(natom))
      allocate(xc(natom))
      allocate(yc(natom))
      allocate(zc(natom))
      allocate(stringPtr(natom))
      do i = 1, natom
        read(10,*)atom_type(i),xc(i),yc(i),zc(i)
        ! trim to remove whitespace
        ! add null char at end to make C-readable
        c_atom(i) = trim(atom_type(i))//char(0)
      enddo
      close(10)
      ! initialize force arrays
      fx(:) = 0d0
      fy(:) = 0d0
      fz(:) = 0d0
      ! initialize system energy
      energy = 0d0

      call f_set_chimes(small)
      
      print*,"fcheck-1"

      call f_init_chimes(trim(param_file) // c_null_char,  0) ! last '0' is the rank of the process
      
      print*,"fcheck-2"
      
      stress(:) = 0d0
      do ns = 1, natom
        stringPtr(ns) = c_loc(c_atom(ns))
      enddo
      
      print*,"fcheck-3"

      call f_calculate_chimes (natom, xc, yc, zc, stringPtr, ca,  &
      &      cb, cc, energy, fx, fy, fz, stress)

      print *, "Success!"
      print '(A,1X, F0.6)', "Energy (kcal/mol):",energy
      print *, "Stress tensors (GPa):"
      print '(A,1X, F15.6)', "s_xx: ",stress(1)*GPa
      print '(A,1X, F15.6)', "s_yy: ",stress(5)*GPa
      print '(A,1X, F15.6)', "s_zz: ",stress(9)*GPa
      print '(A,1X, F15.6)', "s_xy: ",stress(2)*GPa
      print '(A,1X, F15.6)', "s_xz: ",stress(3)*GPa
      print '(A,1X, F15.6)', "s_yz: ",stress(6)*GPa
      print *, "Forces (kcal/mol/A):"
      do i = 1, natom
         print '(F15.6)',fx(i)
         print '(F15.6)',fy(i)
         print '(F15.6)',fz(i)
      enddo
      
#if DEBUG==1

      open (unit = 20, status = 'replace', file='debug.dat')
      write(20,'(F15.6)') energy
      write(20,'(F15.6)') stress(1)*GPa
      write(20,'(F15.6)') stress(5)*GPa
      write(20,'(F15.6)') stress(9)*GPa
      write(20,'(F15.6)') stress(2)*GPa
      write(20,'(F15.6)') stress(3)*GPa
      write(20,'(F15.6)') stress(6)*GPa

      ! Changed format of forces to E15.7 to output the same number of
      ! digits as C (LEF) 08/02/21
      do i = 1, natom
         write(20,'(F15.6)') fx(i)
         write(20,'(F15.6)') fy(i)
         write(20,'(F15.6)') fz(i)
      enddo
      close(20)

#endif 

      end program
