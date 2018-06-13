{\rtf1\ansi\ansicpg1252\cocoartf1504\cocoasubrtf830
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 A research project on simulating the properties of three-dimensional z2 lattice gauge model\
\
1) The equilibrium properties of 3D Lattice Gauge Model (also called Classic Toric Code model)\
\
   Here we apply Metropolis algorithm to simulate the equilibrium properties of the 3D Lattice Gauge Model at different temperature. \
 The process for simulation is:\
   a) Build the lattice for the 3D Lattice Gauge Model\
   b) Initialize the system with random spins\
   c) Performing Monte Carlo steps to equilibrium the system at the initial temperature Ti (Ti comes from input)\
   d) Take measurements: energy, magnetization etc.\
   e) Decrease the temperature by dT\
   f) Repeat c)-e) till the final temperature Tf is reached\
   g) Write results to files\
   \
   The folder contains the following code:\
   \
   to.h: header file with declarations of arrays, functions etc.\
   \
   to_equ.cxx : functions are defined within this file, such as initialization, Monte Carlo update, taking measurement etc.\
   \
   to.cxx: main program, calculating the order parameter, energy for different temperatures with Monte Carlo.\
   \
   to_binaverage.cxx : calculating the averages and error bars of the original measurements after running to.cxx.\
   \
   To compile:\
   g++ to_equ.cxx to.cxx\
   \
   To Run:\
   ./a.out\
   \
   The input file \'93input.txt\'94 has the following format:\
   ll dd\
   init istp mstp nbin\
   filename confname\
   ti tf\
   \
    ll: system size\
   dd: system dimension, here dd=3\
   init: an indicator about whether of not to initialize the system (0 yes, otherwise no)\
   istp: number of Monte Carlo steps to equilibrium the system\
   mstp: number of repetitions of taking measurement within one bin\
   nbin: number of bins\
   filename: the output filename\
   confname: file name for storing the configuration\
   ti: initial temperature\
   tf: final temperature\
  \
   \
   \
   \
  \
 2) The dynamical properties of 3D Lattice Gauge Model (also called Classic Toric Code model)\
 \
   In this part, we employed Simulated Annealing to calculate the dynamical properties of the correlations within the 3D Lattice Gauge Model.\
  \
 The process for simulation is:\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0    a) Build the lattice for the 3D Lattice Gauge Model\
   b) Initialize the system with random spins\
   c) Performing Monte Carlo steps to equilibrium the system at the initial temperature Ti (Ti comes from input)\
   d) Take measurements: energy, magnetization etc.\
   e) Decrease the temperature by dT\
   f)  Perform one Monte Carlo Step\
   e) Repeat d) to f) till final temperature Tf is reached \
   f) Write results to files\
   \
   The folder contains the following code:\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0   \
   corr.h: header file with declarations of arrays, functions etc.\
 \
   corr_quench.cxx: functions are defined within this file, such as Monte Carlo update, decreasing temperature, taking measurements etc.\
  \
   corr.cxx: main program, taking measurement at each step with simulated annealing.\
    \
   corr_binaverage.cxx: calculating the averages and error bars of the original measurements after running corr_quench.cxx.\
\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0    To compile:\
   g++ to_equ.cxx to.cxx\
   \
   To Run:\
   ./a.out\
\
   The input file \'93input.txt\'94 has the following format:\
   ll dd\
   init istp mstp nbin\
   qtime qpnt\
   filename confname\
   ti tf\
\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0 \
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0    ll: system size\
   dd: system dimension, here dd=3\
   init: an indicator about whether of not to initialize the system (0 yes, otherwise no)\
   istp: number of Monte Carlo steps to equilibrium the system at Ti\
   mstp: number of repetitions of taking measurement within one bin\
   nbin: number of bins\
   qtime: total annealing time\
   qpnt:  number of data points to save and write\
   filename: the output filename\
   confname: file name for storing the configuration\
   ti: initial temperature\
   tf: final temperature\
  \
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0 \
}