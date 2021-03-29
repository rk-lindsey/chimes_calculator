! ChIMES Calculator
! Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
! Contributing Author:  Nir Goldman (2020)

module chimes_serial
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_double, c_char, c_null_char
  implicit none

  private
  public :: f_init_chimes, f_set_chimes, f_calculate_chimes, string_to_cstring

  interface

    subroutine f_init_chimes(param_file, rank) bind (C, name='init_chimes')
      import :: c_int, c_char
      character (kind=c_char), dimension(*) :: param_file
      integer(c_int), intent(in) :: rank
    end subroutine f_init_chimes

    subroutine f_set_chimes() bind(C, name='set_chimes')
    end subroutine f_set_chimes

    subroutine f_calculate_chimes (natom, xc, yc, zc, atom_types, ca, cb, cc, energy,&
        & fx, fy, fz, stress) bind(C, name='calculate_chimes')
      import :: c_char, c_ptr, c_int, c_double
      integer(c_int), value :: natom
      real(c_double), intent(in) :: xc(*), yc(*), zc(*)
      type(c_ptr), intent(in) :: atom_types(*)
      real(c_double), intent(in) :: ca(*), cb(*), cc(*)
      real(c_double), intent(inout) :: energy
      real(c_double), intent(inout) :: fx(*), fy(*), fz(*)
      real(c_double), intent(inout) :: stress(*)
    end subroutine f_calculate_chimes

  end interface

contains

  ! Converts a Fortran string to a C-string
  pure function string_to_cstring(string) result(cstring)
    character(len=*), intent(in) :: string
    character(len=(len_trim(string) + len(c_null_char)), kind=c_char) :: cstring

    cstring = trim(string) // c_null_char

  end function string_to_cstring

end module chimes_serial
