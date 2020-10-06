      module wrapper
      use, intrinsic :: ISO_C_binding
      implicit none
      interface 

        subroutine f_chimes_compute_4b_props_fromf90(dr_4b, dist_4b, type1,  & 
       & type2, type3, type4, f4b, stress, sys_ener) & 
       & bind (C, name='chimes_compute_4b_props_fromf90')
          import C_double, C_char, C_ptr
          real(C_double), intent(in) :: dr_4b(6), dist_4b(3,6)
          real(C_double) :: stress(9)
          real(C_double) :: f4b(3,4)
          real(C_double) :: sys_ener
          character(C_char), dimension(80), intent(in) :: type1
          character(C_char), dimension(80), intent(in) :: type2
          character(C_char), dimension(80), intent(in) :: type3
          character(C_char), dimension(80), intent(in) :: type4
        end subroutine f_chimes_compute_4b_props_fromf90
	
        subroutine f_chimes_compute_3b_props_fromf90(dr_3b, dist_3b, type1,  & 
       & type2, type3, f3b, stress, sys_ener) & 
       & bind (C, name='chimes_compute_3b_props_fromf90')
          import C_double, C_char, C_ptr
          real(C_double), intent(in) :: dr_3b(3), dist_3b(3,3)
          real(C_double) :: stress(9)
          real(C_double) :: f3b(3,3)
          real(C_double) :: sys_ener
          character(C_char), dimension(80), intent(in) :: type1
          character(C_char), dimension(80), intent(in) :: type2
          character(C_char), dimension(80), intent(in) :: type3
        end subroutine f_chimes_compute_3b_props_fromf90
	
        subroutine f_chimes_compute_2b_props_fromf90(rij, dr, type1,  & 
       & type2, f2b, stress, sys_ener) & 
       & bind (C, name='chimes_compute_2b_props_fromf90')
          import C_double, C_char, C_ptr
          real(C_double), intent(in) :: rij, dr(3)
          real(C_double) :: stress(9)
          real(C_double) :: f2b(2,3)
          real(C_double) :: sys_ener
          character(C_char), dimension(80), intent(in) :: type1
          character(C_char), dimension(80), intent(in) :: type2
        end subroutine f_chimes_compute_2b_props_fromf90
	
        subroutine f_get_chimes_epot(sys_ener) & 
      &   bind (C, name='get_chimes_epot')
          import C_double
          real(C_double) :: sys_ener
        end subroutine f_get_chimes_epot
	
        subroutine f_set_chimes_epot(sys_ener) & 
      &   bind (C, name='set_chimes_epot')
          import C_double
          real(C_double) :: sys_ener
        end subroutine f_set_chimes_epot
	
        subroutine f_set_chimes() bind (C, name='set_chimes')
        end subroutine f_set_chimes
	
        subroutine f_init_chimes(rank) bind (C, name='init_chimes') 
          import C_int
          integer(C_int), intent(in) :: rank
        end subroutine f_init_chimes
	
        subroutine f_chimes_read_params(param_file) & 
      &   bind (C, name='chimes_read_params')
          import C_char
          character (kind=C_char), dimension(*) :: param_file
        end subroutine f_chimes_read_params
	
        function f_get_chimes_2b_order () result (order2b) &
          bind (C, name='get_chimes_2b_order')
          import :: C_int
          integer (C_int) :: order2b
        end function f_get_chimes_2b_order
	
        function f_get_chimes_3b_order () result (order3b) &
          bind (C, name='get_chimes_3b_order')
          import :: C_int
          integer (C_int) :: order3b
        end function f_get_chimes_3b_order
	
        function f_get_chimes_4b_order () result (order4b) &
          bind (C, name='get_chimes_4b_order')
          import :: C_int
          integer (C_int) :: order4b
        end function f_get_chimes_4b_order
	
        function f_get_chimes_max_2b_cutoff () result (rcut_2b) &
          bind (C, name='get_chimes_max_2b_cutoff')
          import :: C_double
          real(C_double) :: rcut_2b
        end function f_get_chimes_max_2b_cutoff
	
        function f_get_chimes_max_3b_cutoff () result (rcut_3b) &
          bind (C, name='get_chimes_max_3b_cutoff')
          import :: C_double
          real(C_double) :: rcut_3b
        end function f_get_chimes_max_3b_cutoff
	
        function f_get_chimes_max_4b_cutoff () result (rcut_4b) &
          bind (C, name='get_chimes_max_4b_cutoff')
          import :: C_double
          real(C_double) :: rcut_4b
        end function f_get_chimes_max_4b_cutoff
	
      end interface
      contains
      pure function string2Cstring (string) result (C_string)
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
