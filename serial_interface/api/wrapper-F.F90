      module wrapper
      use, intrinsic :: ISO_C_binding
      implicit none
      interface 
        subroutine f_calculate_chimes (natom, xc, yc, zc, cptr, ca,  &
      &            cb, cc, energy, fx, fy, fz, stress)  &
      &            bind (C, name='calculate_chimes')
          use, intrinsic :: ISO_C_binding, only : C_char, C_ptr, C_int, C_double
          type(c_ptr), dimension(*)  :: cptr(natom)
          integer(C_int) :: natom  
          real(C_double) :: xc(natom), yc(natom), zc(natom)
          real(C_double) :: fx(natom), fy(natom), fz(natom)
          character(C_char), dimension(80) :: c_atom(natom)
          real(C_double) :: ca(3), cb(3), cc(3) 
          real(C_double) :: energy
          real(C_double) :: stress(9)
        end subroutine f_calculate_chimes
        subroutine f_set_chimes() bind (C, name='set_chimes')
        end subroutine f_set_chimes
        subroutine f_init_chimes(param_file, nlayer) & 
      &            bind (C, name='init_chimes') 
          import C_int, C_char
          integer(C_int), intent(in) :: nlayer
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
        forall (i = 1:n)
           C_string(i) = string(i:i)
        end forall
        C_string(n+1) = C_NULL_CHAR
      end function string2Cstring
      end module
