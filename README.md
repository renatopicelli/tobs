# Fluid Topology Optimization using TOBS approach

This repository uses the more open source libraries possible to optimize fluid flow using the TOBS method.

## Author of this repository:

Bruno Caldas de Souza (University of Sao Paulo) , bruno.caldas@usp.br;

## Author of the TOBS method:

Raghavendra Sivapuram (University of California), rsivapur@eng.ucsd.edu;

Renato Picelli (University of São Paulo), rpicelli@usp.br;

## Requirements (conda is recommended):

1. Python 3.5;
2. Fenics 2018.1.0;
3. Dolfin Adjoint 2018.1.0;
4. Octave to read the .m files;
5. CPlex* for Integer optimization.

*Unfortunately this one is proprietary to IBM, however there is an free student version available.

## How to use:

First, try to install the conda environment by:
```
conda env create -f environment.yml -n fluid_tobs
```

Then install Fenics 2018.1.0 and Dolfin-Adjoint 2018.1.0 from source:
```
git clone --branch 2018.1.0 git@bitbucket.org:fenics-project/dolfin.git
```

```
git clone --branch 2018.1.0 git@bitbucket.org:dolfin-adjoint/pyadjoint.git
```

Then, activate the environment:
```
conda activate fluid_tobs
```

In order to run the basic simulation type in your conda environment:
```
python stokes_tobs.py
```

## The results:

The results can be seen in the output folder that is created and should look like this:
![Alt Text](output/example.gif)

## Publications:

[1] R. Sivapuram, R. Picelli, Topology optimization of binary structures using Integer
Linear Programming, Finite Elements in Analysis and Design 139 (2018) 49–61

[2] R. Sivapuram, R. Picelli, Y.M. Xie, Topology optimization of binary microstructures involving various non-volume constraints, Computational Materials Science (2018) 154, 405-425
