
## symmFast
Reference implementation of a fast bilinear algorithm for the symmetric matrix times vector kernel.

### Theory
The information content of a symmetric matrix is `n(n+1)/2`, but even standard high performance algorithms for symmetric matrix times vector multiplication have bilinear cost `n^2` since they ignore symmetry altogether.
However, it is possible to construct algorithms for structured matrices that have an information-content-optimal bilinear cost.

For example, [Solomonik and Demmel's 2021 paper](https://www.degruyter.com/document/doi/10.1515/cmam-2019-0075/html) propose an algebraic algorithm with bilinear complexity `n(n+1)/2` for the symmetric matrix times vector product.
Especially in settings where multiplication is more expensive than addition, like for complex numbers, such algorithms may offer appreciable speedups over canonical structure-unaware algorithms.

Such speedups are particularly valuable since the symmetric matrix times vector multiply is one of the most important kernels in numerical linear algebra, being the basis for iterative linear and eigenpair solvers on symmetric matrices.

### This project
Existing linear algebra libraries like [LAPACK](https://github.com/Reference-LAPACK/lapack) and various BLAS implementations have the benefit of decades of improvements and optimizations, meaning it would be unfair to compare a relatively unoptimized reference version of the fast bilinear algorithm to them.
This repository serves not only as reference implementation for the fast bilinear algorithms and schemes for their parallelization, but also as a fair testbed between algorithms, since the canonical and fast bilinear algorithms are implemented with the same level of optimizations.
