! ChIMES Calculator
! Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
! Contributing Author:  Nir Goldman (2020)
! Modern Fortran (Fortran 2003) interface: Bálint Aradi (2021)

!> Provides simplified (a.k.a. serial) ChIMES interface using modern Fortran.
!>
!> Note: Currently the simplified ChIMES interface only allows the creation of one ChIMES instance
!> without the possibility of destroying that instance. Therefore, the Fortran interface only
!> allows the creation of one ChIMES-calculator.
!>
module chimescalc_serial08
  use, intrinsic :: iso_c_binding, only : c_ptr, c_double, c_char, c_loc, c_null_char
  use chimescalc_serial, only : f_init_chimes, f_set_chimes, f_calculate_chimes, string2Cstring
  implicit none

  private
  public :: ChimesCalc, ChimesCalc_init

  !> Default precision for reals
  integer, parameter :: dp = c_double

  !> Internal variable counting nr. of active calculator instances (only 1 allowed)
  integer :: active_instances_ = 0

  ! Workaround: gfortran 7, gfortran 8
  ! Deferred length character arrays are not handled correclty, leading to trashed memory
  ! Therefore a fixed character length array is used for the species names, instead
  integer, parameter :: species_name_length = 20


  !> Chimes calculator
  type :: ChimesCalc
    private

    !> Whether the calculator is active
    logical :: is_active_ = .false.

    !> C representation of the atom type names.
    character(len=species_name_length+1), pointer :: c_atom_types_(:) => null()

    !> C pointers referencing the atom type names
    type(c_ptr), allocatable :: c_atom_type_ptrs_(:)

  contains

    procedure :: set_atom_types => ChimesCalc_set_atom_types
    procedure :: calculate => ChimesCalc_calculate
    final :: ChimesCalc_final

  end type ChimesCalc


contains

  !> Initializes ChIMES calculator
  subroutine ChimesCalc_init(this, paramfile, nrank, small)

    !> Initialized instance on exit.
    type(ChimesCalc), intent(inout) :: this

    !> Name of the parameter file to use for the initialization
    character(*), intent(in) :: paramfile

    !> Rank parameter (?)
    integer, intent(in) :: nrank
    integer, intent(in) :: small

    if (this%is_active_) then
      error stop "This ChimesCalc instance has already been initialized"
    end if
    if (active_instances_ /= 0) then
      error stop "Only one chimes_calc instance is allowed"
    end if
    this%is_active_ = .true.
    active_instances_ = 1

    call f_set_chimes(small)
    call f_init_chimes(string2Cstring(paramfile), nrank)

  end subroutine ChimesCalc_init


  !> Destroys ChIMES calculator
  subroutine ChimesCalc_final(this)
    type(ChimesCalc), intent(inout) :: this

    this%is_active_ = .false.
    ! As the C-interface does not offer finalization, we do not allow for a new ChIMES instance
    ! even if the exising one had been destroyed.
    !active_instances_ = 0
    if (associated(this%c_atom_types_)) then
      deallocate(this%c_atom_types_)
    end if

  end subroutine ChimesCalc_final


  !> Sets the atom types (needs the type name of each atom)
  !>
  !> This has to be called in order to set up the atom types, before calculate() can be invoked.
  !> It must be called only once, unless the atom types change (e.g. atoms are reordered).
  !>
  subroutine ChimesCalc_set_atom_types(this, atom_types)

    !> Instance.
    class(ChimesCalc), intent(inout) :: this

    !> Type name of each atom. Shape: [nr_of_atoms]
    character(*), intent(in) :: atom_types(:)

    integer :: natom, iat

    if (associated(this%c_atom_types_)) then
      deallocate(this%c_atom_types_)
    end if
    if (allocated(this%c_atom_type_ptrs_)) then
      deallocate(this%c_atom_type_ptrs_)
    end if

    natom = size(atom_types)
    !allocate(character(len=(len(atom_types) + 1)) :: this%c_atom_types_(natom))
    allocate(this%c_atom_types_(natom))
    allocate(this%c_atom_type_ptrs_(natom))
    do iat = 1, natom
      !this%c_atom_types_(iat) = trim(atom_types(iat)) // c_null_char
      this%c_atom_types_(iat) = &
          & atom_types(iat)(1 : min(len_trim(atom_types(iat)), species_name_length)) // c_null_char
      this%c_atom_type_ptrs_(iat) = c_loc(this%c_atom_types_(iat))
    end do

  end subroutine ChimesCalc_set_atom_types


  !> Calculates (adds) various ChIMES contributions
  !>
  !> Note: You have to call set_atom_types() once before you can call this method.
  !>
  subroutine ChimesCalc_calculate(this, coords, latvecs, energy, forces, stress)

    !> Instance.
    class(ChimesCalc), intent(inout) :: this

    !> Coordinates of the atoms. Shape: [3, nr_of_atoms]
    real(dp), intent(in) :: coords(:,:)

    !> Lattice vectors. Shape: [3, 3], first index runs over x,y,z, second over lattice vectors.
    real(dp), intent(in) :: latvecs(:,:)

    !> Variable which should be increased by the ChIMES energy.
    real(dp), intent(inout) :: energy

    !> Forces, which ChIMES contribution should be added to. Shape: [3, nr_of_atoms].
    real(dp), intent(inout) :: forces(:,:)

    !> Stress tensor, which the ChIMES contribution should be added to. Shape: [3, 3].
    real(dp), intent(inout), contiguous, target :: stress(:,:)

    real(dp), pointer :: stress_1d(:)

    stress_1d(1 : size(stress)) => stress
    call f_calculate_chimes(size(this%c_atom_types_), coords(1, :), coords(2, :), coords(3, :),&
        & this%c_atom_type_ptrs_, latvecs(:, 1), latvecs(:, 2), latvecs(:, 3), energy,&
        & forces(1, :), forces(2, :), forces(3, :), stress_1d)

  end subroutine ChimesCalc_calculate


end module chimescalc_serial08
