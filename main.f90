!-----------------------------------------------------------------------------
!! Program: schrodinger_solution
!! By: Noah Egger
!!
!! This modules calls all defined subroutines needed to retrieve the 3 lowest 
!! eigenvalues, their corresponding eigenfunctions, prints those energies, and 
!! writes to a file that contains the normalized wave functions along with the 
!! x-positions where the points were sampled. We implement this process for 3 
!! distinct potentials; the infinite well, the harmonic oscillator, and the 
!! woods saxon.
!-----------------------------------------------------------------------------
program schrodinger_solution 

use types
use read_write, only: read_input, write_probability_density, print_energies, write_woods_saxon_energies
use qm_solver, only: sample_box, solve_infinite_well, analytic_infinite_well, solve_harmonic_oscillator,&
    analytic_harmonic_oscillator, solve_woods_saxon
implicit none

integer :: n_points
real(dp) :: length, radius

integer, parameter :: n_energies = 3
!real(dp) :: energies(1:n_energies), analytic_energies(1:n_energies)
real(dp), allocatable :: wave_functions(:,:), energies(:), analytic_energies(:)
real(dp), allocatable :: x_vector(:)
real(dp), parameter :: r_min = 2._dp, r_max = 10._dp
allocate(x_vector(1:n_points))

! Prompts user for input: length of potential, radius of woods saxon potential, and numer of sampling points
call read_input(length, n_points, radius)

! Constructs 1D array of equally spaced sampling points
call sample_box(length, n_points, x_vector)

! Solving particle in a box for infinite well
call solve_infinite_well(x_vector, n_points, length, n_energies, energies, wave_functions)

! Call routine to calculate eigenenergies analytically
call analytic_infinite_well(length, n_energies, analytic_energies)

! Call routine to print analytic and numerical energies
call print_energies('Infinite Well', energies, analytic_energies, n_energies)

! Writes to a file the normalized wave functions and the corresponding x locations
call write_probability_density('infinite_well_wf.dat', x_vector, wave_functions)

! Solving harmonic oscillator
call solve_harmonic_oscillator(x_vector, n_points, length, n_energies, energies, wave_functions)

! Calculate first 3 energies analytically
call analytic_harmonic_oscillator(length, n_energies, analytic_energies)

! Print numerical and analytic energies
call print_energies('Harmonic oscillator', energies, analytic_energies, n_energies)

! Writes normalized squared wave functions along with sampling locations
call write_probability_density('harmonic_oscillator_wf.dat', x_vector, wave_functions)

! Solving Woods Saxon, finding eigenfunctions and eigenenergies
call solve_woods_saxon(x_vector, n_points, length, radius, n_energies, energies, wave_functions)

! Writes normalized squared wave functions along with sampling locations
call write_probability_density('woods_saxon_wf.dat', x_vector, wave_functions)

! Woods Saxon lowest 3 energies as a function of radius from 2fm to 8fm
call write_woods_saxon_energies('woods_saxon_ener.dat', x_vector, n_points, n_energies, length, r_min, r_max)

end program schrodinger_solution