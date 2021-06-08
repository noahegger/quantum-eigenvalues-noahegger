!-----------------------------------------------------------------------
!Module: hamiltonian
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This module houses the subroutines necessary for constructing the kinetic
!! energy matrix as well as the potential energy matrix. One must construct
!! both the diagonal and off-diagonal terms for both matrices, as well
!! as for the distinct potentials, i.e. the infinite well, harmonic oscillator,
!! and the woods saxon. More specifically, the potential energy will be a
!! diagonal matrix and the kinetic energy will be a tridiagonal matrix, thus
!! making the hamiltonian a tridiagonal matrix as well. Additionally, this 
!! module holds the subroutine that normalizes the eigenfunctions received.
!!
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! harmonic_v
!! kinetic_energy
!! normalize
!! woods_saxon_v
!!----------------------------------------------------------------------
module hamiltonian
use types
implicit none

! Since more than one subroutine in this module will use the value for
! hbar and mass, it would be a good idea to define them here as 
! parameters

real(dp), parameter :: hbar = 197.3
real(dp), parameter :: mass = 939
real(dp), parameter :: v_0 = 50
real(dp), parameter :: a = 0.2

private
public :: calculate_harmonic_v, kinetic_energy, normalize, &
& woods_saxon_v

contains

!-----------------------------------------------------------------------
!! Subroutine: calculate_harmonic_v
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine constructs the diagonal elements for the harmonic 
!! oscillator potential defined as V = ((h/2m)^2)*x^2 where h is hbar and 
!! x is the position we are sampling for the potential. This 1-D array is 
!! then used to construct the full matrix. 
!!----------------------------------------------------------------------
!! Input:
!! 
!! x_vector(:)      real        array containing all sampling points between -L and +L
!! n_points         integer     number of sampling points
!!
!-----------------------------------------------------------------------
!! Output:
!!
!! potential(:)     real        array containing the potential evaluated at corresponding x
!!
!-----------------------------------------------------------------------

subroutine calculate_harmonic_v(x_vector, n_points, potential)
    implicit none
    real(dp), intent(out), allocatable :: potential(:)
    real(dp), intent(in) :: x_vector(:)
    integer, intent(in) :: n_points
    integer :: i, size_v, size_x
	if(allocated(potential)) deallocate(potential)
    allocate(potential(1:n_points))

    ! Fill the potential array which contains the diagonal values
    do i = 1, n_points
        potential(i) = ((hbar**2)/(2._dp*mass))*x_vector(i)**2
    enddo
end subroutine calculate_harmonic_v

!-----------------------------------------------------------------------
!! Subroutine: kinetic_energy
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine constructs and fills the kinetic energy diagonal and
!! off diagonal arrays for use in the qm_solver.f90 module, in which the 
!! hamiltonians are constructed and eigenvectors and eigenvalues are found. 
!!----------------------------------------------------------------------
!! Input:
!!
!! length           real        half length of the potential well
!! n_points         integer     number of sampling points
!!
!-----------------------------------------------------------------------
!! Output:
!!
!! ke_diag          real        array containing kinetic energy diagonal terms
!! ke_offdiag       real        array containing kinetic energy off diagonal terms
!!
!-----------------------------------------------------------------------

subroutine kinetic_energy(n_points, length, ke_diag, ke_offdiag)
    implicit none
    real(dp), intent(in) :: length
    integer, intent(in) :: n_points
    real(dp), allocatable, intent(out) :: ke_diag(:), ke_offdiag(:)
    real(dp) :: dx
    integer :: i, j
    if(allocated(ke_diag)) deallocate(ke_diag)
    if(allocated(ke_offdiag)) deallocate(ke_offdiag)
    allocate(ke_diag(1:n_points))
    allocate(ke_offdiag(1:n_points-1))

    ! Step size
    dx = 2._dp*length/(n_points-1._dp)

    ! Construct kinetic energy diagonal
    do i = 1, n_points
        ke_diag(i) = (hbar**2)/mass/(dx**2)
    enddo

    ! Construct kinetic energy off diagonal
    do j = 1, n_points-1
        ke_offdiag(j) = -0.5_dp*(hbar**2)/mass/(dx**2)
    enddo

    ! Arrays will called by another subroutine to solve TISE and find
    ! eigenvalues and eigenvectors 

end subroutine kinetic_energy

!-----------------------------------------------------------------------
!! Subroutine: normalize
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine receives the un-normalized wave functions, wave_functions,
!! which are generated from the hamiltonian using the dstev() routine in lapack.
!! This subroutine takes a 2-D array of wave functions and normalizes the 
!! first three eigenfunctions (first three columns). The normalization constant 
!! is put into an array. The eigenfunctions are then each multipled by their 
!! corresponding normalization constant and a new 2-D array is constructed containg 
!! the normalized eigenfunctions as the columns of the matrix.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! wave_functions(:,:)      real        2-D array containing un-normalized eigenfunctions
!! length                   real        half length of the potential well
!! n_points                 integer     number of sampling points
!! n_energies               integer     number of lowest eigenvalues and corresponding eigenfunctions 
!!
!-----------------------------------------------------------------------
!! Output:
!!
!! normal_function(:,:)     real        array containing the first three normalized eigenfunctions
!!
!-----------------------------------------------------------------------

subroutine normalize(wave_functions, n_energies, length, n_points, normal_function)
    implicit none
    real(dp), intent(in) :: wave_functions(:,:), length
    integer, intent(in) :: n_points, n_energies
    real(dp), allocatable :: squared(:), normalization(:)
    real(dp), intent(out), allocatable :: normal_function(:,:)
    real(dp) :: sum, dx
    integer :: i, j, k, l

    allocate(squared(1:n_energies))
    allocate(normalization(1:n_energies))
    allocate(normal_function(1:n_points, 1:n_energies))

    ! Step size
    dx = 2._dp*length/(n_points-1._dp)

    ! Square the column vector element wise and add up the individiual terms
    ! Construct new array where each element is the sum of squared elements for 
    ! a particular column.

    do i = 1, n_energies
        sum = 0._dp
        do j = 1, n_points
            sum = sum + wave_functions(j,i)**2
        enddo
        ! Store sum of squared column elements in array called squared
        squared(i) = sum
    enddo

    ! Construct array that holds normalization constant for each column vector
    do j = 1, n_energies
        normalization(j) = 1/sqrt(squared(j)*dx)
    enddo

    ! Multiplying each column element by the normalization constant will
    ! result in a normalized wave function array.

    do k = 1, n_energies
        do l = 1, n_points
            normal_function(l,k) = wave_functions(l,k)*normalization(k)
        enddo
    enddo
end subroutine normalize

!-----------------------------------------------------------------------
!! Subroutine: woods_saxon_v
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine constructs the diagonal elements of the woods saxon potential
!! described as V = -v_0/(1+exp((x-R)/a)) where h is hbar, x is the position
!! at which we are sampling the potential, and R is the radius of the potential.
!! This 1-D array is then sent back to the subroutine that calls it to then
!! construct the full matrix describing the potential. 
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector(:)      real        array containing all sampling points between -L and +L
!! radius           real        radius of the woods saxon potential
!! n_points         integer     number of sampling points
!!
!-----------------------------------------------------------------------
!! Output:
!!
!! potential(:)     real        array containing the potential evaluated at corresponding x
!!
!-----------------------------------------------------------------------

subroutine woods_saxon_v(x_vector, n_points, radius, potential)
    implicit none
    real(dp), intent(out), allocatable :: potential(:)
    real(dp), intent(in) :: x_vector(:)
    real(dp), intent(in) :: radius
    integer, intent(in) :: n_points
    integer :: i, size_x
    if(allocated(potential)) deallocate(potential)
    allocate(potential(1:n_points))

    size_x = size(x_vector)

    ! Check size of the x_vector array matches the number of sampling points. 
    if (size_x /= n_points) then
        print*, 'x_vector is not proper length'
        stop
    endif

    ! Construct potential array which contains the diagonal values for woods
    ! saxon
    do i = 1, n_points
        potential(i) = -v_0/(1+exp((abs(x_vector(i))-radius)/a))
    enddo
end subroutine woods_saxon_v

end module hamiltonian