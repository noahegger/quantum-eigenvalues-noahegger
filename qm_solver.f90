!-----------------------------------------------------------------------
!Module: qm_solver
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This module addresses the relevant physics to solve the time-independent
!! schrodinger equation for the three distinct potentials. While some of 
!! the bare mathematics is contained in eigen_solver, this module
!! calls the relevant functions and subroutines needed to build the 
!! hamiltonian (where some auxiliary subroutines and functions exist in
!! the hamiltonian module) and solves the TISE for each distinct potential. 
!! While in addition to solving the TISE for the three distinct potentials, 
!! this module also constructs the sampling points for which all potentials 
!! are functions of, i.e. it constructs the x position array.
!! 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! sample_box
!! solve_infinite_well
!! analytic_infinite_well
!! solve_harmonic_oscillator
!! analytic_harmonic_oscillator
!! solve_woods_saxon
!! vary_r_saxon
!!----------------------------------------------------------------------
module qm_solver
use types
use hamiltonian, only : calculate_harmonic_v, kinetic_energy, normalize, & 
& woods_saxon_v 
implicit none

! Natural units
real(dp), parameter :: hbar = 197.3
real(dp), parameter :: mass = 939

private
public solve_infinite_well, sample_box, analytic_infinite_well, solve_harmonic_oscillator, &
    analytic_harmonic_oscillator, solve_woods_saxon, vary_r_saxon

contains

!-----------------------------------------------------------------------
!! Subroutine: sample_box
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine takes in user input length of potential well as well 
!! as the number of sampling points. Constructs the array containing 
!! the x-positions that discretizes the wave function.
!!----------------------------------------------------------------------
!! Input:
!!
!! length		real 		half length of potential well
!! n_points		integer 	number of sampling points
!!
!-----------------------------------------------------------------------
!! Output:
!!
!! x_vector(:)	real 		array containing sampling points between -L and L
!!
!-----------------------------------------------------------------------
subroutine sample_box(length, n_points, x_vector)
    implicit none
    integer, intent(in) :: n_points
    real(dp), intent(in) :: length
    real(dp), intent(out), allocatable :: x_vector(:)
    real(dp) :: dx
    integer :: i
    if(allocated(x_vector)) deallocate(x_vector)
    allocate(x_vector(1:n_points))

    ! Step size
    dx = 2._dp*length/(n_points-1)

    ! Build x array, start at -L until reaching L
    i = 1
    do while(-length + dx*(i-1) <= length)
        x_vector(i) = -length + dx*(i-1)
        i = i + 1
    enddo

    ! Now we have an array of equally spaced points from -L to L

end subroutine sample_box

!-----------------------------------------------------------------------
!! Subroutine: solve_infinite_well
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine constructs the kinetic energy tri-diagonal matrix
!! as well as the potential diagonal matrix by calling other subroutines
!! such as those contained in the hamiltonian module. The subroutine
!! retrieves 3 1-D arrays, specifically the kinetic energy diagonal
!! array, the kinetic energy off diagonal array, and the potential
!! energy diagonal array. It then calls dstev() (a routine in lapack)
!! to find the eigenfunctions and eigenenergies associated with the hamiltonian. 
!! Finally, the subroutine uses the normalize routine in the hamiltonian 
!! module to normalize the first three wave functions associated with the
!! 3 lowest energies. The lowest 3 energies and the corresponding normalized
!! wave functions are sent back to main.f90.
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector     real        array containing points from -L to L
!! n_points     integer     number of sampling points
!! length       real        half length of potential well
!! n_energies   integer     number of eigenvalues and corresponding eigenfunction
!!
!-----------------------------------------------------------------------
!! Output:
!!
!! energies(:)              real        array containing lowest 3 energies 
!! normal_functions(:,:)    real        array with first 3 normalized eigenfunctions
!!
!-----------------------------------------------------------------------
subroutine solve_infinite_well(x_vector, n_points, length, n_energies, energies, normal_functions)
    implicit none
    integer, intent(in) :: n_points, n_energies
    real(dp), intent(in) :: x_vector(:)
    real(dp), intent(in) :: length
    real(dp), intent(out), allocatable :: energies(:), normal_functions(:,:)
    real(dp), allocatable :: ke_diag(:), ke_offdiag(:)
    real(dp), allocatable :: dstev_work(:), wave_functions(:,:)
    integer :: info, i

    if(allocated(energies)) deallocate(energies)

    allocate(ke_diag(1:n_points))
    allocate(ke_offdiag(1:n_points-1))
    allocate(energies(1:n_energies))
    allocate(wave_functions(1:n_points,1:n_points))
    allocate(dstev_work(1:2*n_points - 2))
    allocate(normal_functions(1:n_points,1:n_energies))

    ! Hamiltonian is just the kinetic energy as V = 0.
    ! Call kinetic_energy which builds diagonal and off diagonal term arrays
    call kinetic_energy(n_points, length, ke_diag, ke_offdiag)

    ! Send the hamiltonian diagonal and off diagonal into the lapack
    ! routine in order to return the array of eigenenergies as well as an n_points by n_points matrix containing the wave 
    ! functions partially normalized along the columns. 
    call dstev('V', n_points, ke_diag, ke_offdiag, wave_functions, n_points, dstev_work, info)

    ! Normalize wave functions
    call normalize(wave_functions, n_energies, length, n_points, normal_functions)

    ! Make array of size n_energies (3 in this case) to hold eigenenergies
    ! computed by dstev() 
    do i = 1, n_energies
        energies(i) = ke_diag(i)
    enddo
end subroutine solve_infinite_well

!-----------------------------------------------------------------------
!! Subroutine: analytic_infinite_well
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine computes the analytic eigenenergies for the infinite 
!! potential well for use in comparison to the numerical method.
!!----------------------------------------------------------------------
!! Input:
!!
!! length       real        half length of potential well
!!
!-----------------------------------------------------------------------
!! Output:
!!
!! analytic(:)  real        array contianing 3 lowest analytic eigenenergies
!!
!-----------------------------------------------------------------------
subroutine analytic_infinite_well(length, n_energies, analytic)
    implicit none
    real(dp), intent(in) :: length
    integer, intent(in) :: n_energies
    real(dp), intent(out), allocatable :: analytic(:)
    integer :: i
    allocate(analytic(1:n_energies))

    ! Calculate the ground state, first excited state, and second excited
    ! state using analytic formula. 
    do i = 1, n_energies
        analytic(i) = (i**2)*(hbar**2)*(pi**2)/(8._dp*mass*(length**2))
    enddo
end subroutine analytic_infinite_well

!-----------------------------------------------------------------------
!! Subroutine: solve_harmonic_oscillator
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine does the same thing as solve_infinite_well, except
!! for the harmonic oscillator potential.
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector     real        array of equally spaced sampling points from -L to L
!! n_points     integer     number of sampling points
!! n_energies   integer     number of eigenvalues/eigenfunctions
!! length       real        half length of potential well
!!
!-----------------------------------------------------------------------
!! Output:
!!
!!
!-----------------------------------------------------------------------
subroutine solve_harmonic_oscillator(x_vector, n_points, length, n_energies, energies, normal_functions)
    implicit none
    integer, intent(in) :: n_points, n_energies
    real(dp), intent(in) :: x_vector(:)
    real(dp), intent(in) :: length
    real(dp), intent(out), allocatable :: energies(:), normal_functions(:,:)
    real(dp), allocatable :: potential(:), ke_diag(:), ke_offdiag(:)
    real(dp), allocatable :: dstev_work(:), wave_functions(:,:), h_diag(:)
    integer :: info, i, j

    if(allocated(energies)) deallocate(energies)
    ! Need to allocate arrays which construct the hamiltonian 
    ! as well as arrays need to find eigenvalues and eigenvectors
    ! of the hamiltonian. 
    allocate(ke_diag(1:n_points))
    allocate(ke_offdiag(1:n_points-1))
    allocate(potential(1:n_points))
    allocate(energies(1:n_energies))
    allocate(wave_functions(1:n_points,1:n_points))
    allocate(dstev_work(1:2*n_points - 2))
    allocate(normal_functions(1:n_points,1:n_energies))
    allocate(h_diag(1:n_points))

    ! Diagonal terms of the potential
    call calculate_harmonic_v(x_vector, n_points, potential)

    ! Diagonal and off diagonal terms of kinetic energy
    call kinetic_energy(n_points, length, ke_diag, ke_offdiag)

    ! Diagonal element of hamiltonian
    do j = 1, n_points
        h_diag(j) = ke_diag(j) + potential(j)
    enddo

    ! Use lapack routine to return the array of eigenvalues (replaces h_diag) 
    ! as well as an n_points by n_points matrix containing the wave functions
    ! partially normalized along the columns. 
    call dstev('V', n_points, h_diag, ke_offdiag, wave_functions, n_points, dstev_work, info)

    ! Normalize the eigenfunctions
    call normalize(wave_functions, n_energies, length, n_points, normal_functions)

    ! Construct array of eigenvalues returned from dstev().
    do i = 1, n_energies
        energies(i) = h_diag(i)
    enddo
end subroutine solve_harmonic_oscillator

!-----------------------------------------------------------------------
!! Subroutine: analytic_harmonic_oscillator
!-----------------------------------------------------------------------
!! By: Noah Egger
!! 
!! Computes the analytic eigenenergies for the harmonic oscillator for
!! use in comparison to the numerical method.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! length       real        half length of potential well
!!
!-----------------------------------------------------------------------
!! Output:
!!
!! analytic(:)  real        array containing analytic eigenenergies
!!
!-----------------------------------------------------------------------
subroutine analytic_harmonic_oscillator(length, n_energies, analytic)
    implicit none
    real(dp), intent(in) :: length
    integer, intent(in) :: n_energies
    real(dp), intent(out), allocatable :: analytic(:)
    integer :: i
    allocate(analytic(1:n_energies))

    ! Calculate eigenenergies with analytic formula
    do i = 1, n_energies
        analytic(i) = (i-0.5_dp)*hbar**2/mass
    enddo
end subroutine analytic_harmonic_oscillator

!-----------------------------------------------------------------------
!! Subroutine: solve_woods_saxon
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine does the same thing as solve_infinite_well, except
!! for the woods saxon potential.
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector         real        array of equally spaced points from -L to L
!! n_points         integer     number of sampling points
!! n_energies       integer     number of eigenvalues/eigenvectors (3 lowest energy states)
!! length           real        half length of potential well
!! radius           real        radius of woods saxon potential
!!
!-----------------------------------------------------------------------
!! Output:
!!
!! energies(:)              real        array containing eigenenergies
!! normal_functions(:,:)    real        array of normalized eigenfunctions
!!
!-----------------------------------------------------------------------
subroutine solve_woods_saxon(x_vector, n_points, length, radius, n_energies, energies, normal_functions)
    implicit none
    integer, intent(in) :: n_points, n_energies
    real(dp), intent(in) :: x_vector(:)
    real(dp), intent(in) :: length, radius
    real(dp), intent(out), allocatable :: energies(:), normal_functions(:,:)
    real(dp), allocatable :: potential(:), ke_diag(:), ke_offdiag(:)
    real(dp), allocatable :: dstev_work(:), wave_functions(:,:), h_diag(:)
    integer :: info, i, j

    if(allocated(energies)) deallocate(energies)

    ! Allocate arrays to hold information used for constructing hamiltonian
    ! and thus finding eigenenergies/eigenvectors
    allocate(ke_diag(1:n_points))
    allocate(ke_offdiag(1:n_points-1))
    allocate(potential(1:n_points))
    allocate(energies(1:n_energies))
    allocate(wave_functions(1:n_points,1:n_points))
    allocate(dstev_work(1:2*n_points - 2))
    allocate(normal_functions(1:n_points,1:n_energies))
    allocate(h_diag(1:n_points))

    ! Calculate diagonal terms of potential
    call woods_saxon_v(x_vector, n_points, radius, potential)

    ! Build KE arrays holding diagonal and off diagonal terms
    call kinetic_energy(n_points, length, ke_diag, ke_offdiag)

    ! Construct diagonal elements of hamiltonian
    do j = 1, n_points
        h_diag(j) = ke_diag(j) + potential (j)
    enddo

    ! Use lapack routine to return array of eigenvalues (diagonalized hamiltonian)
    ! and a 2D array containing wave functions
    call dstev('V', n_points, h_diag, ke_offdiag, wave_functions, n_points, dstev_work, info)

    ! Normalize wave functions
    ! Return eigenvalues in a 1D array
    call normalize(wave_functions, n_energies, length, n_points, normal_functions)
    do i = 1, n_energies
        energies(i) = h_diag(i)
    enddo

end subroutine solve_woods_saxon

!-----------------------------------------------------------------------
!! Subroutine: vary_r_saxon
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine finds first 3 eigenenergies corresponding to a 
!! particular radius. Results are later plotted to see how each
!! energy changes as a funtion of radius.
!!----------------------------------------------------------------------
!! Input:
!!
!! length           real        half length of potential well
!! radius           real        radius of the woods saxon potential
!! x_vector         real        array containing equally space sampling points from -L and +L
!! n_points         integer     number of sampling points
!!
!-----------------------------------------------------------------------
!! Output:
!!
!! energies(:)      real        eigenenergies array
!!
!-----------------------------------------------------------------------
subroutine vary_r_saxon(x_vector, length, n_points, n_energies, radius, energies)
    implicit none
    real(dp), intent(in) :: length, radius, x_vector(:)
    integer, intent(in) :: n_points, n_energies
    real(dp), intent(out), allocatable :: energies(:)
    real(dp), allocatable :: potential(:), ke_diag(:), ke_offdiag(:)
    real(dp), allocatable :: h_diag(:), z(:,:), d(:)
    integer :: i, info, j
    allocate(h_diag(1:n_points))
    allocate(z(1:n_points, 1:n_points))
    allocate(d(1:2*n_points-2))
    allocate(potential(1:n_points))
    allocate(ke_diag(1:n_points))
    allocate(ke_offdiag(1:n_points-1))
    allocate(energies(1:n_energies))

    ! Construct diagonal terms of potential
    call woods_saxon_v(x_vector, n_points, radius, potential)

    ! Construct diagonal and off diagonal arrays for KE
    call kinetic_energy(n_points, length, ke_diag, ke_offdiag)

    ! Construct diagonal terms of hamiltonian
    do j = 1, n_points
        h_diag(j) = ke_diag(j) + potential(j)
    enddo

    ! Use lapack routine to return array of eigenvalues (diagonalized hamiltonian)
    ! and a 2D array containing wave functions
    call dstev('N', n_points, h_diag, ke_offdiag, z, n_points, d, info)

    ! Return eigenvalues in 1D array
    do i = 1, n_energies
        energies(i) = h_diag(i)
    enddo
    ! Read_write module will loop through various radius values
    ! and routine energies array containing eigenenergies for given radius
end subroutine vary_r_saxon

end module qm_solver