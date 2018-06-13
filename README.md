# Simulating-3D-Lattice-Gauge-Model
A research project on simulating the properties of three-dimensional z2 lattice gauge model
1) The equilibrium properties of 3D Lattice Gauge Model (also called Classic Toric Code model)

   Here we apply Metropolis algorithm to simulate the equilibrium properties of the 3D Lattice Gauge Model at different temperature. 
   
   to.h: header file with defination of arrays, functions etc.
   
   to_equ.cxx : functions are defined within this file.
   
   to.cxx: main program, calculating the order parameter, energy for different temperatures with Monte Carlo.
   
   to_binaverage.cxx : calculating the averages and error bars of the original measurements after running to.cxx.
  
 2) The dynamical properties of 3D Lattice Gauge Model (also called Classic Toric Code model)
 
   In this part, we employed Simulated Annealing to calculate the dynamical properties of the correlations within the 3D Lattice Gauge Model.
   
   corr.h: header file with dininitions of arrays, functions etc.
   
   corr_quench.cxx: functions are defined within this file.
   
   corr.cxx: main program, calculating each measurement with simulated annealing.
   
   corr_binaverage.cxx: calculating the averages and error bars of the original measurements after running corr_quench.cxx.
