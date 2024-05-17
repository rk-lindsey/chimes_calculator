! ChIMES Calculator
! Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
! Contributing Author:  Nir Goldman (2020)

      module chimescalc_serial
      use, intrinsic :: ISO_C_binding
      implicit none

      interface

        function f_chimes_open_instance() & 
      &          bind (C, name='chimes_open_instance')
          use, intrinsic :: ISO_C_binding, only :  C_ptr
          type(C_ptr) :: f_chimes_open_instance
        end function f_chimes_open_instance

        subroutine f_chimes_close_instance(handle) &
      &            bind (C, name='chimes_close_instance')
          use, intrinsic :: ISO_C_binding, only : C_ptr
          type(C_ptr), value :: handle
        end subroutine f_chimes_close_instance

        subroutine f_calculate_chimes_instance(handle, natom, xc, yc, zc, cptr, ca,  &
      &            cb, cc, energy, fx, fy, fz, stress)  &
      &            bind (C, name='calculate_chimes_instance')
          use, intrinsic :: ISO_C_binding, only : C_char, C_ptr, C_int, C_double
          implicit none
          integer(C_int), value :: natom
          type(c_ptr), dimension(*)  :: cptr(natom)
          type(C_ptr), value :: handle
          real(C_double) :: xc(natom), yc(natom), zc(natom)
          real(C_double) :: fx(natom), fy(natom), fz(natom)
          character(C_char), dimension(80) :: c_atom(natom)
          real(C_double) :: ca(3), cb(3), cc(3)
          real(C_double) :: energy
          real(C_double) :: stress(9)
        end subroutine f_calculate_chimes_instance

        subroutine f_calculate_chimes (natom, xc, yc, zc, cptr, ca,  &
      &            cb, cc, energy, fx, fy, fz, stress)  &
      &            bind (C, name='calculate_chimes')
          use, intrinsic :: ISO_C_binding, only : C_char, C_ptr, C_int, C_double
          implicit none
          integer(C_int), value :: natom
          type(c_ptr), dimension(*)  :: cptr(natom)
          real(C_double) :: xc(natom), yc(natom), zc(natom)
          real(C_double) :: fx(natom), fy(natom), fz(natom)
          character(C_char), dimension(80) :: c_atom(natom)
          real(C_double) :: ca(3), cb(3), cc(3)
          real(C_double) :: energy
          real(C_double) :: stress(9)
        end subroutine f_calculate_chimes

        subroutine f_set_chimes_instance(handle, small) bind &
      &            (C, name='set_chimes_serial_instance')
          use, intrinsic :: ISO_C_binding, only : C_ptr, C_int
          implicit none
          integer(C_int), value :: small
          type(C_ptr), value :: handle
        end subroutine f_set_chimes_instance

        subroutine f_set_chimes(small) bind &
      &            (C, name='set_chimes_serial')
          use, intrinsic :: ISO_C_binding, only : C_ptr, C_int
          implicit none
          integer(C_int), value :: small
        end subroutine f_set_chimes

        subroutine f_init_chimes_instance(handle, param_file, rank) &
      &            bind (C, name='init_chimes_serial_instance')
          use, intrinsic :: ISO_C_binding, only : C_ptr, C_int, C_char
          integer(C_int), value :: rank
          character (kind=C_char), dimension(*) :: param_file
          type(C_ptr), value :: handle
        end subroutine f_init_chimes_instance

        subroutine f_init_chimes(param_file, rank) &
      &            bind (C, name='init_chimes_serial')
          import C_int, C_char
          integer(C_int), intent(in) :: rank
          character (kind=C_char), dimension(*) :: param_file
        end subroutine f_init_chimes

      end interface

      contains

      function string2Cstring (string) result (C_string)
        use, intrinsic :: ISO_C_binding, only : C_char, C_NULL_CHAR
        character (len=*), intent(in) :: string
        character (len=1, kind=C_char) :: C_string (len_trim(string)+1)
        integer :: i, n
        n = len_trim (string)
        do i = 1, n
           C_string(i) = string(i:i)
        end do
        C_string(n+1) = C_NULL_CHAR
      end function string2Cstring

      end module chimescalc_serial
