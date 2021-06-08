!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Rodrigo Navarro Perez/Noah Egger
!!
!! Contains subroutines and functions related to reading input from the
!! user and  writing output into a text file
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! read_advanced_input
!!----------------------------------------------------------------------
!! Included functions:
!!
!! read_real
!! read_integer
!-----------------------------------------------------------------------
module read_write
use types
use qm_solver, only : vary_r_saxon
implicit none

private
public :: read_input, read_real, read_integer, write_probability_density, &
print_energies, write_woods_saxon_energies

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! Rodrigo Navarro Perez/Noah Egger
!!
!! Displays a message describing what the program does and the expected
!! input. After that it uses the `read_real` and `read_integer`
!! functions to assign values to the different parameters.
!!----------------------------------------------------------------------
!! Output:
!!
!! n_points     integer     number of grid points the discretized wave function
!! length       real        length of the box
!! radius       real        radius of the Woods-Saxon potential
!-----------------------------------------------------------------------
subroutine read_input(length, n_points, radius)
    implicit none
    real(dp), intent(out) :: length, radius
    integer, intent(out) ::  n_points

    print*, 'This program finds numerical solutions to the 1-D '
    print*, 'Schrodinger equation for three distinct potentials.'
    print*, 'Please enter values for the length of the potential'
    print*, 'well, in addition to the number of sampling points'
    print*, '"N", and the radius of the Woods-Saxon Potential "r".'

    ! Function are called below to ask for this user input
    ! Each input is verified to be of correct type in a do loop
    length = read_real('length L of the potential well:')
    n_points = read_integer('number of points N on the lattice:')
    radius = read_real('radius of the Woods Saxon potential:')    

end subroutine read_input

!-----------------------------------------------------------------------
!! Function: read_integer
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Verifies the input given by the user is of the correct type.
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Output:
!!
!! x        integer     A positive non negative number given by the user
!-----------------------------------------------------------------------
integer function read_integer(name) result(x)
    implicit none
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print *, 'Provide a nonzero positive value for the '//trim(name)//':'
  
     do
        read(*,'(a)',iostat=ierror) string
        ! If input is nonempty, proceed
        if(string.ne.'') then
            read(string,*,iostat=ierror) x
            ! If input can be made into a number, proceed
            if (ierror == 0) then
                ! if number is positive, we can exit the loop.
                if (x > 0) exit       
                    print *, "'"//trim(string)//"'"// 'cannot be negative or zero, please provide a positive number' 
            else
                print *, "'"//trim(string)//"'"//' is not a number, please provide a number'
            endif          
        else
            print *, 'that was an empty input, please provide a positive, non-zero number'
        
        endif
    enddo
end function read_integer
!-----------------------------------------------------------------------
!! Function: read_real
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Verifies the input given by the user is of the correct type. 
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Output:
!!
!! x        real        A positive non negative number given by the user
!-----------------------------------------------------------------------
real(dp) function read_real(name) result(x)
    implicit none
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror
 
    print *, 'Provide a nonzero positive value for the '//trim(name)//':'

    do
        read(*,'(a)',iostat=ierror) string
        ! If input is not empty, proceed
        if(string.ne.'') then
            read(string,*,iostat=ierror) x
            ! If input can be made into a number, proceed
            if (ierror == 0) then
                ! if number is positive, we can exit the loop.
                if (x > 0) exit       
                    print *, "'"//trim(string)//"'"// 'cannot be negative or zero, please provide a positive number' 
            else
                print *, "'"//trim(string)//"'"//' is not a number, please provide a number'
            endif          
        else
            print *, 'that was an empty input, please provide a positive, non-zero, number'
        
        endif
    enddo
end function read_real

!-----------------------------------------------------------------------
!! Function: write_probability_density
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine writes the probabilty density of the normalized
!! eigenfunctions that correspond to the eigenenergies. This routine is
!! called for the infinite well, harmonic oscillator, and woods saxon potentials. 
!! The probability densities are written along with their corresponding x
!! value from x_vector array.
!!----------------------------------------------------------------------
!! Input:
!!
!! file_name             character   A string with the name of the file to be created and written to 
!! x_vector              real        array containing equally spaced points from -L to L
!! wave_functions(:,:)   real        2-D array containing normalized eigenfunctions  
!!
!!----------------------------------------------------------------------
!! Output:
!!
!! 
!-----------------------------------------------------------------------


subroutine write_probability_density(file_name, x_vector, wave_functions)
    implicit none
    character(len=*), intent(in) :: file_name
    real(dp), intent(in) :: x_vector(:), wave_functions(:,:)
    integer :: unit, i, size_x
    size_x = size(x_vector)

    ! Open file, write column headers, write data
    open(newunit=unit,file=trim(file_name))
    write(unit,'(4a28)') 'x', 'ground state', '1st excited', '2nd excited'

    ! Record x value and corresponding probability density at that location
    do i = 1, size_x
    ! may need to indent
    write(unit,*) x_vector(i), wave_functions(i,1)**2, wave_functions(i,2)**2, wave_functions(i,3)**2
    enddo
    close(unit)
    
end subroutine write_probability_density

!-----------------------------------------------------------------------
!! Function: print_energies
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine prints the numerical and analytic eigenenergies.
!!----------------------------------------------------------------------
!! Input:
!!
!! name          character   A string with a title of the model the energies correspond to     
!!----------------------------------------------------------------------
!! Output:
!!
!! numerical(:)  array       Contains numerical energies
!! analytic(:)   array       Contains analytic energies
!!
!-----------------------------------------------------------------------


subroutine print_energies(name, numerical, analytic, n_energies)
    implicit none
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: numerical(:), analytic(:)
    integer, intent(in) :: n_energies
    integer :: i

    ! Check analytic and numerical energy arrays are same size
    if (size(numerical) /= size(analytic)) then
        print*, "arrays size don't match in print_energies"
        stop
    endif

    print*, 'Comparing numerical and analytic solutions in'
    print*, trim(name) 
    print'(a9,2a15)', 'number', 'numerical', 'analytic'
    
    do i = 1, n_energies
        print*, i, numerical(i), analytic(i)
    enddo
end subroutine print_energies

!-----------------------------------------------------------------------
!! Function: write_woods_saxon_energies
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine writes to a file the eigenenergies for the woods saxon 
!! potential as a function of changing radius. This subroutine interates 
!! through radii from 2 fm to 10 fm with a step of "dr". For each radius, 
!! a subroutine vary_saxon_r is called, returning the eigenenergies.
!! For each specific radius, the corresponding eigenenergies and radius 
!! are written to file. 
!!----------------------------------------------------------------------
!! Input:
!!
!! file_name     character   A string with a brief description of the value being asked for
!! n_points      integer     number of sampling points      
!! length        real        half length of potential well
!! r_min         real        minimum radius
!! r_max         real        maximum radius
!! x_vector(:)   real        array containing equall spaced sampling points between -L and +L
!!----------------------------------------------------------------------
!! Output:
!!
!-----------------------------------------------------------------------

subroutine write_woods_saxon_energies(file_name, x_vector, n_points, n_energies, length, r_min, r_max)
    implicit none
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: n_points, n_energies
    real(dp), intent(in) :: length, r_min, r_max, x_vector(:)
    real(dp), allocatable :: energies(:) 
    real(dp) :: dr, radius 
    integer :: unit2
    allocate(energies(1:n_energies))

    ! Step size
    dr = (r_max-r_min)/(n_points-1._dp)

    ! Open unit, create file, write headers for columns. 
    open(newunit=unit2, file=trim(file_name))
    write(unit2,'(4a28)') 'Radius', 'Ground', 'First Excited', 'Second Excited'

    ! Radius starts at r_min 

    radius = r_min
    do 
        if(radius>r_max) exit
    ! For the first radius we call vary_r_saxon and get first three energies back, then
    ! iterate through all radii until r_max is reached. 
        call vary_r_saxon(x_vector, length, n_points, n_energies, radius, energies)

    ! Write energies corresponding to radius to file
        write(unit2, *) radius, energies(1), energies(2), energies(3)
    ! Increase radius to radius + dr, then run through loop again. 
        radius = radius + dr
    enddo


end subroutine write_woods_saxon_energies


end module read_write
