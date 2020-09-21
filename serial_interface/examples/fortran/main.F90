      program test_F_api
      use wrapper
      use, intrinsic :: ISO_C_binding
      implicit none
      integer io_num
      double precision, parameter :: GPa = 6.9479 ! convert kcal/mol.A^3 to GPa
      character(C_char), dimension(80) :: c_file
      character(C_char), dimension(80) :: dummy_var
      character(80) :: coord_file, param_file
      CHARACTER ( len = 100 ) :: wq_char
      integer :: i, j, k, l, natom, ns
      integer(C_int) :: nlayer
      real(C_double) ::   lx, ly, lz
      real(C_double) :: stress(9)
      real(C_double) :: energy
      real(C_double) :: ca(3), cb(3), cc(3)
      real(C_double), allocatable :: xc(:), yc(:), zc(:), fx(:), fy(:), fz(:)
      character(C_char), dimension(80), allocatable, target :: atom_type(:) 
      character(len=80) :: atom
      TYPE(C_PTR), allocatable :: stringPtr(:)

      io_num = IARGC()
      if (io_num .lt. 2) then
        print*,"To run: ./test_F.x <parameter file> <xyz config. file>"
        print*,"Exiting code.\n"
        STOP
      endif
      call getarg(1, wq_char)
      read(wq_char,*)param_file
      call getarg(2, wq_char)
      read(wq_char,*)coord_file
      open (unit=10, status='old', file=coord_file)
      read(10,*)natom
      read(10,*)ca(1),ca(2),ca(3),cb(1),cb(2),cb(3),cc(1),cc(2),cc(3)

      allocate(fx(natom))
      allocate(fy(natom))
      allocate(fz(natom))
      allocate(atom_type(natom))
      allocate(xc(natom))
      allocate(yc(natom))
      allocate(zc(natom))
      allocate(stringPtr(natom))
      do i = 1, natom
        read(10,*)atom,xc(i),yc(i),zc(i)
        atom_type(i) = atom//C_NULL_CHAR
      enddo
      close(10)
      ! initialize force arrays
      fx(:) = 0d0
      fy(:) = 0d0
      fz(:) = 0d0
      ! initialize system energy
      energy = 0d0
      call f_set_chimes()
      nlayer = 1 
      c_file = string2Cstring(param_file)
      call f_init_chimes(c_file, nlayer)
      stress(:) = 0d0
      do ns = 1, natom
        stringPtr(ns) = c_loc(atom_type(ns))
      enddo
      call f_calculate_chimes (natom, xc, yc, zc, stringPtr, ca,  &
      &      cb, cc, energy, fx, fy, fz, stress)
      print*,'Total energy= ',energy
      open (unit = 20, status = 'replace', file='output_libf_serial.xyz')
      write(20,*)natom
      stress(:) = GPa*stress(:)
      write(20,*)lx, ly, lz, stress(1:9), energy
      do i = 1, natom
        write(20,*)atom_type(i),xc(i),yc(i),zc(i),fx(i), fy(i), fz(i)
      enddo
      close(20)
      end program
