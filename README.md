# Solving Schrödinger Equation with Eigenvalues: Noah Egger

# Program Goal

This program solves the time-independent Schrödinger equation (TISE henceforth) for three distinct potentials: the infinite well, harmonic oscillator, and Woods Saxon potential. For the infinite well, 'V = 0'. For the harmonic oscillator, 'V = (hbar^2/2m)x^2'. For the Woods Saxon, 'V = -V_0/(1+exp((x-R)/a))'. Upon compiling the program, the user is asked to provide the half length of the potential well 'L' as well as the number of sampling points 'N' and the radius of the Woods Saxon potential. 

The length of the well is discretized into equally space points from '-L' to 'L', thus allowing us to discretize the wave function. The spacing in the points 'dx' is equal to 2L/(n_points-1). Using the potential and kinetic energy functions allows us to construct the hamiltonian. The kinetic energy is a tri-diagonal matrix, with '(h_bar^2)/(mass dx^2)' on the diagonals and '-(h_bar^2)/(2 mass dx^2)' along the off-diagonals. The potential energy matrix, which is distinct for each model, has values along the diagonal diagonal. The Hamiltonian is the sum of these matrices. We use the routine 'dstev()' within 'lapack' to solve the TISE and return the eigenenergies and corresponding eigenfunctions of any particular hamiltonian. This program returns the first 3 lowest eigenenergies for each potential model. Additionally, the program calculates the eigenenergies analytically for usage as comparison to the numerical calculations. For each potential model, the corresponding 3 eigenfunctions are normalized and written to files along with the x sampling positions from -L to L. We also write to these same files the ground, 1st, and 2nd excited states' energy densities. These files are `infinite_well_wf.dat`, `harmonic_oscillator_wf.dat`, and `woods_saxon_wf.dat` which are all plotted in the `.ipynb` file.

Finally, the 3 lowest energies of the Woods Saxon potential are calculated as a function of varying radius (from 2fm to 10fm). The step size 'dr' is given by '8/(n-1)'. These energies are written to a file `woods_saxon_ener.dat` and plotted as a function of radius in the `.ipynb` file.

# Directions for Usage

Navigate to your directory and ensure all relevant files are contained. Type "make" into the terminal and press "enter". The files will compile and the executable 'woods_saxon' will be created. Type "./woods_saxon" into the terminal and press "enter". The command line will prompt you to enter the length of the potential well "L", the number of sampling points "N, and the radius of the woods saxon potential "radius". After each prompt, enter a number and press "enter". After this, the terminal will print the numerical and analytic energies for the lowest 3 energy states. In addition, the program will create the four files `infinite_well_wf.dat`, `harmonic_oscillator_wf.dat`, `woods_saxon_wf.dat`, and `woods_saxon_ener.dat`. Open `plots_analysis.ipynb` with jupyter notebook and run the program.

# Files Contained

`read_write.f90` reads command line input, writes the probability densities to files, and writes the energies as a function of radius for woods saxon to a file.  
`hamiltonian.f90` constructs the kinetic energy diagonals and off diagonals, potential diagonals , and normalization of the eigenfunctions for each model.   

`qm_solver.f90` calls  subroutines to construct the hamiltonian for each model, then retreives the eigenvalues and eigenfunctions for each model. Also constructs sampling points array and finds the three lowest energies via the analytic method. Addtionally, finds the three lowest eigenergies as a function of varying radius for the woods saxon potential.  

`types.f90` contains the argument types, integers and reals for calculation.   

`plots_analysis_assignment04.ipynb` illustrates plotting results. 

`makefile` handles compilation order. Typing "make" into terminal compiles the .f90 files and creates executable.   

`main.f90` houses the main calls to run the program, including the write routines. 

