Matlab code for topology optimization with binary design variables and sequential integer linear programming.

Authors:
Raghavendra Sivapuram (University of California), rsivapur@eng.ucsd.edu;
Renato Picelli (University of São Paulo), rpicelli@usp.br;

Implementation of the Topology Optimization of Binary Structures (TOBS) method.

Numerical features:
1. Problem linearization;
2. Move limits (constraints relaxation);
3. Sensitivity filtering;
4. Integer programming*.

*This code uses the mixed-integer linear programming solver "intlinprog" from Matlab. For better performance, we recommend the use of the CPLEX library, free to download at IBM website. After installing CPLEX, the installation path and "options.Optimizer" at TopOpt/TOBS.m must be edited.

Publications:
[1] R. Sivapuram, R. Picelli, Topology optimization of binary structures using Integer
Linear Programming, Finite Elements in Analysis and Design 139 (2018) 49–61
[2] R. Sivapuram, R. Picelli, Y.M. Xie, Topology optimization of binary microstructures involving various non-volume constraints, Computational Materials Science (2018) 154, 405-425

Disclaimer:
The authors reserve all rights but do not guaranty that the code is free from errors. Furthermore, they shall not be liable in any event caused by the use of the program.
