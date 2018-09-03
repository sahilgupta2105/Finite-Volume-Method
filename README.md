# Lid Driven Cavity Solver based on FVM

Implementation of lid driven cavity solver based on SIMPLE algorithm. The method used to solve system of linear equations is Gauss-Seidel. Moreover, the 3D version of the code implements central differencing scheme for discretization of velocity terms and comes in two variants: 3D and 2D arrays. This refers to the storage to the fluid properties such as pressure, velocity etc. The 2D arrays is more efficient as compared to the 3D arrays storage scheme.
