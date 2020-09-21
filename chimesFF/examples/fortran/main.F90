! Code author: Nir Goldman (2020)
      program test_F_api
      use wrapper
      use, intrinsic :: ISO_C_binding
      implicit none
      integer io_num
      double precision, parameter :: GPa = 6.9479 ! convert kcal/mol.A^3 to GPa
      character(C_char), dimension(80) :: c_file 
      character(C_char), dimension(80) :: dummy_var
      character(2000) :: coord_file, param_file
      CHARACTER ( len = 2000 ) :: wq_char
      integer :: i, j, k, l, natom
      integer(C_int) :: rank 
      real(C_double) ::   lx, ly, lz, ldummy
      double precision :: vol
      real(C_double) :: stress(9)
      double precision :: epot
      real(C_double) :: sys_ener
      double precision, allocatable :: ftot(:,:)
      real(C_double), allocatable :: xc(:), yc(:), zc(:)
      character(C_char), dimension(80), allocatable :: atom_type(:)
      character(C_char), dimension(80) :: type1, type2, type3, type4
      character(80) :: atype2b(2)
      real(C_double) :: f2b(3,2), f3b(3,3), f4b(3,4)
      integer(C_int) :: order2b, order3b, order4b
      real(C_double) :: rcut_2b, rcut_3b, rcut_4b
      real(C_double) :: xij, yij, zij, dr(3)
      real(C_double) :: xik, yik, zik
      real(C_double) :: xil, yil, zil
      real(C_double) :: xjk, yjk, zjk
      real(C_double) :: xjl, yjl, zjl
      real(C_double) :: xkl, ykl, zkl
      real(C_double) :: dr_3b(3), dist_3b(3,3)
      real(C_double) :: dr_4b(6), dist_4b(3,6)
      real(C_double) :: rij
      type(C_ptr) :: c_rij 
      
      io_num = command_argument_count()
      if (io_num .lt. 2) then
        print*,"To run: ./test_F.x <parameter file> <xyz config. file>"
        print*,"Exiting code.\n"
        STOP
      endif
      call GET_COMMAND_ARGUMENT(1, wq_char)
      param_file = trim(wq_char)
      call GET_COMMAND_ARGUMENT(2, wq_char)
      coord_file = trim(wq_char)

      open (unit=10, status='old', file=coord_file)
      read(10,*)natom
      read(10,*)lx,ldummy,ldummy,ldummy,ly,ldummy,ldummy,ldummy,lz
      vol = lx*ly*lz

      allocate(ftot(3,natom))
      allocate(atom_type(natom))
      allocate(xc(natom))
      allocate(yc(natom))
      allocate(zc(natom))
      do i = 1, natom
        read(10,*)atom_type(i),xc(i),yc(i),zc(i)
      enddo
      ! initialize ftot array
      do j = 1, natom
        do i = 1, 3
          ftot(i,j) = 0.0
        enddo
      enddo
      ! initialize f2b array, index order reversed automatically
      do j = 1, 2
        do i = 1, 3
          f2b(i,j) = 0.0
        enddo
      enddo
      close(10)
      call f_set_chimes()
      rank = 0
      call f_init_chimes(rank)
      c_file = string2Cstring(param_file)
      call  f_chimes_read_params(c_file)
      order2b=f_get_chimes_2b_order();
      order3b=f_get_chimes_3b_order();
      order4b=f_get_chimes_4b_order();
      rcut_2b = f_get_chimes_max_2b_cutoff();
      sys_ener = 0d0
      stress(:) = 0d0
      !call f_set_chimes_epot(sys_ener)
      do i = 1, natom-1
        do j = i+1, natom
          xij = xc(i) - xc(j)
          xij = xij - lx*nint(xij/lx)
          dr(1) = xij   
          yij = yc(i) - yc(j)
          yij = yij - lx*nint(yij/ly)
          dr(2) = yij   
          zij = zc(i) - zc(j)
          zij = zij - lx*nint(zij/lz)
          dr(3) = zij   
          rij = sqrt(xij*xij + yij*yij + zij*zij)
          f2b(:,1) = ftot(:,i)
          f2b(:,2) = ftot(:,j)
          
          type1 = string2Cstring(atom_type(i))
          type2 = string2Cstring(atom_type(j))
          if (rij .le. rcut_2b) then
          ! f2b, stress tensor, epot are all cumulative
             call f_chimes_compute_2b_props_fromf90(rij, dr, type1, & 
      &           type2, f2b, stress, sys_ener)
          endif
          !save results back in ftot
          ftot(:,i) = f2b(:,1)
          ftot(:,j) = f2b(:,2)
        enddo
      enddo
      print*,'2B energy= ',sys_ener
! 3-body calculation
      if (order3b .gt. 0) then
        rcut_3b = f_get_chimes_max_3b_cutoff();
        ! initialize f3b array, index order reversed automatically
        do j = 1, 3
          do i = 1, 3
            f3b(i,j) = 0.0
          enddo
        enddo
        do i = 1, natom-2
          do j = i+1, natom-1
            do k = j+1, natom
       ! compute relative coordinates and apply minimum image PBC
       ! order in chimesFF is ij, ik, jk
       ! ij pairs 
              xij = (xc(i) - xc(j));
              xij = xij - lx*nint(xij/lx);
              yij = (yc(i) - yc(j));
              yij = yij - ly*nint(yij/ly);
              zij = (zc(i) - zc(j));
              zij = zij - lz*nint(zij/lz);
              dist_3b(1,1) = xij;
              dist_3b(2,1) = yij;
              dist_3b(3,1) = zij;
              dr_3b(1) = sqrt(xij*xij + yij*yij + zij*zij);
       ! ik pairs
              xik = (xc(i) - xc(k));
              xik = xik - lx*nint(xik/lx);
              yik = (yc(i) - yc(k));
              yik = yik - ly*nint(yik/ly);
              zik = (zc(i) - zc(k));
              zik = zik - lz*nint(zik/lz);
              dist_3b(1,2) = xik;
              dist_3b(2,2) = yik;
              dist_3b(3,2) = zik;
              dr_3b(2) = sqrt(xik*xik + yik*yik + zik*zik);
       ! jk pairs 
              xjk = (xc(j) - xc(k));
              xjk = xjk - lx*nint(xjk/lx);
              yjk = (yc(j) - yc(k));
              yjk = yjk - ly*nint(yjk/ly);
              zjk = (zc(j) - zc(k));
              zjk = zjk - lz*nint(zjk/lz);
              dist_3b(1,3) = xjk;
              dist_3b(2,3) = yjk;
              dist_3b(3,3) = zjk;
              dr_3b(3) = sqrt(xjk*xjk + yjk*yjk + zjk*zjk);
              f3b(:,1) = ftot(:,i)
              f3b(:,2) = ftot(:,j)
              f3b(:,3) = ftot(:,k)
          
              type1 = string2Cstring(atom_type(i))
              type2 = string2Cstring(atom_type(j))
              type3 = string2Cstring(atom_type(k))
              if (dr_3b(1) .le. rcut_3b) then
                if (dr_3b(2) .le. rcut_3b) then
                  if (dr_3b(3) .le. rcut_3b) then
                    ! f2b, stress tensor, epot are all cumulative
                    call f_chimes_compute_3b_props_fromf90(dr_3b, dist_3b, type1, & 
      &                  type2, type3, f3b, stress, sys_ener)
                  endif
                endif
              endif
            !save results back in ftot
              ftot(:,i) = f3b(:,1)
              ftot(:,j) = f3b(:,2)
              ftot(:,k) = f3b(:,3)
            enddo
          enddo
        enddo
        print*,'3B energy= ',sys_ener
      endif
! 4-body calculation
      if (order4b .gt. 0) then
        rcut_4b = f_get_chimes_max_3b_cutoff();
        ! initialize f3b array, index order reversed automatically
        do j = 1, 4
          do i = 1, 3
            f4b(i,j) = 0.0
          enddo
        enddo
        do i = 1, natom-3
          do j = i+1, natom-3
            do k = j+1, natom-1
              do l = k+1, natom
                ! compute relative coordinates and apply minimum image PBC
                ! order in chimesFF is: ij, ik, il, jk, jl, kl
                ! ij pairs 
                xij = (xc(i) - xc(j));
                xij = xij - lx*nint(xij/lx);
                yij = (yc(i) - yc(j));
                yij = yij - ly*nint(yij/ly);
                zij = (zc(i) - zc(j));
                zij = zij - lz*nint(zij/lz);
                dist_4b(1,1) = xij;
                dist_4b(2,1) = yij;
                dist_4b(3,1) = zij;
                dr_4b(1) = sqrt(xij*xij + yij*yij + zij*zij);
                ! ik pairs 
                xik = (xc(i) - xc(k));
                xik = xik - lx*nint(xik/lx);
                yik = (yc(i) - yc(k));
                yik = yik - ly*nint(yik/ly);
                zik = (zc(i) - zc(k));
                zik = zik - lz*nint(zik/lz);
                dist_4b(1,2) = xik;
                dist_4b(2,2) = yik;
                dist_4b(3,2) = zik;
                dr_4b(2) = sqrt(xik*xik + yik*yik + zik*zik);
                ! il pairs 
                xil = (xc(i) - xc(l));
                xil = xil - lx*nint(xil/lx);
                yil = (yc(i) - yc(l));
                yil = yil - ly*nint(yil/ly);
                zil = (zc(i) - zc(l));
                zil = zil - lz*nint(zil/lz);
                dist_4b(1,3) = xil;
                dist_4b(2,3) = yil;
                dist_4b(3,3) = zil;
                dr_4b(3) = sqrt(xil*xil + yil*yil + zil*zil);
                ! jk pairs 
                xjk = (xc(j) - xc(k));
                xjk = xjk - lx*nint(xjk/lx);
                yjk = (yc(j) - yc(k));
                yjk = yjk - ly*nint(yjk/ly);
                zjk = (zc(j) - zc(k));
                zjk = zjk - lz*nint(zjk/lz);
                dist_4b(1,4) = xjk;
                dist_4b(2,4) = yjk;
                dist_4b(3,4) = zjk;
                dr_4b(4) = sqrt(xjk*xjk + yjk*yjk + zjk*zjk);
                ! jl pairs 
                xjl = (xc(j) - xc(l));
                xjl = xjl - lx*nint(xjl/lx);
                yjl = (yc(j) - yc(l));
                yjl = yjl - ly*nint(yjl/ly);
                zjl = (zc(j) - zc(l));
                zjl = zjl - lz*nint(zjl/lz);
                dist_4b(1,5) = xjl;
                dist_4b(2,5) = yjl;
                dist_4b(3,5) = zjl;
                dr_4b(5) = sqrt(xjl*xjl + yjl*yjl + zjl*zjl);
                ! kl pairs 
                xkl = (xc(k) - xc(l));
                xkl = xkl - lx*nint(xkl/lx);
                ykl = (yc(k) - yc(l));
                ykl = ykl - ly*nint(ykl/ly);
                zkl = (zc(k) - zc(l));
                zkl = zkl - lz*nint(zkl/lz);
                dist_4b(1,6) = xkl;
                dist_4b(2,6) = ykl;
                dist_4b(3,6) = zkl;
                dr_4b(6) = sqrt(xkl*xkl + ykl*ykl + zkl*zkl);
                ! set up 4b force array, f4b
                f4b(:,1) = ftot(:,i)
                f4b(:,2) = ftot(:,j)
                f4b(:,3) = ftot(:,k)
                f4b(:,4) = ftot(:,l)
                type1 = string2Cstring(atom_type(i))
                type2 = string2Cstring(atom_type(j))
                type3 = string2Cstring(atom_type(k))
                type4 = string2Cstring(atom_type(l))
                if (dr_4b(1) .le. rcut_4b) then
                  if (dr_4b(2) .le. rcut_4b) then
                    if (dr_4b(3) .le. rcut_4b) then
                      if (dr_4b(4) .le. rcut_4b) then
                        if (dr_4b(5) .le. rcut_4b) then
                          if (dr_4b(6) .le. rcut_4b) then
                          ! f2b, stress tensor, epot are all cumulative
                            call f_chimes_compute_4b_props_fromf90(dr_4b, dist_4b, type1, & 
      &                          type2, type3, type4, f4b, stress, sys_ener)
                          endif
                        endif
                      endif
                    endif
                  endif
                endif
              !save results back in ftot
                ftot(:,i) = f4b(:,1)
                ftot(:,j) = f4b(:,2)
                ftot(:,k) = f4b(:,3)
                ftot(:,l) = f4b(:,4)
              enddo
            enddo
          enddo
        enddo
        print*,'4B energy= ',sys_ener
      endif
      open (unit = 20, status = 'replace', file='output_libf.xyz')
      write(20,*)natom
      stress(:) = GPa*(stress(:))/vol
      write(20,*)lx, ly, lz, stress(1:9), sys_ener
      do i = 1, natom
        write(20,*)atom_type(i),xc(i),yc(i),zc(i),ftot(:,i)
      enddo
      close(20)
      end program
