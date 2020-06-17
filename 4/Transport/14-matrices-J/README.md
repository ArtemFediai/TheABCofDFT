# 14 matrices

Code to calculate transport properties of 1D-quantum system connected to two periodic electrodes.

**How to use**: 

1.  14_matrices_main.py is the main program which manages the software stored in modules/.
2.  The system has 3 options for executing: -c, -b, -p they stand for: config, build and plot
*  Use -b and the path to the yaml file to build your empirical tight binding 14 matrices. They are saved in the path specified in the <settings>.yml file. 
*  Use -c and the yaml settings file to compute transmission (T) and Density of States (DOS) for your system.
*  Use -p to plot with "mathplotlib" both quantities.
3.  The location of the 14 input matrices **(all numpy arrays stored as *.dat)** has to be specified in the <settings>.yml file.
4.  The <settings>.yml is stored in the settings folder. 
5.  The folder /tests contains examples where the software has been tested (single tight-binding chain in different configurations, 4-zigzag tube).
  