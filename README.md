
## symmFast
Reference implementation of a fast bilinear algorithm for the symmetric matrix times vector kernel in distributed memory.

### Theory
The information content of a symmetric matrix is `n(n+1)/2`, but even standard high performance algorithms for symmetric matrix times vector multiplication have bilinear cost `n^2` since they ignore symmetry altogether.
However, it is possible to construct algorithms for structured matrices that have an information-content-optimal bilinear cost.

For example, [Solomonik and Demmel's 2021 paper](https://www.degruyter.com/document/doi/10.1515/cmam-2019-0075/html) proposes an algebraic algorithm with bilinear complexity `n(n+1)/2` for the symmetric matrix times vector product, though it requires an increased number of additions.
Especially in settings where multiplication is more expensive than addition, like for complex numbers, such algorithms may offer appreciable speedups over canonical structure-unaware algorithms.

Such speedups are particularly valuable since the symmetric matrix times vector multiply is one of the most important kernels in numerical linear algebra, being the basis for iterative linear and eigenpair solvers on symmetric matrices.

### Implementation and contents
The project was originally going to use a larger PETSc framework, which remains in the repository.
However, in the interest of development speed, the project has been contained (with its own makefile) to `bench/ex2.cpp`.

This folder also contains the slurm script used to perform strong scaling benchmarking (`strong_scaling.sh`), the output files of said runs for real and complex arrays (`output_<type>.sh`), and a Python processing script.

### Installation
Requires: working C++ compiler and PETSc.

Use the makefile in (`bench`).

