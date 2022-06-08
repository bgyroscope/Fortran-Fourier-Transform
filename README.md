# Fortran-Fourier-Transform
Use RK4 solver to find Fourier Transform

Osc_main.f90 and Osc_solver.f90 are the two Fourtran files that once compiled will produce executable that preforms the actual calculations. 

The main outputs are: 
* r_w.dat &rarr; the omega grid used. 
* fft_final.dat &rarr; the final Fourier transform
* fft_exact.dat &rarr; the exact result 
* simpleRK4.dat &rarr; the output of the RK4 differential equation solver. 

This example has been done with the function $f(t) $ 
