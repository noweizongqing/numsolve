# Linear solving: relation to mizchi/numbt

## The overlap is real

`mizchi/numbt` is a NumPy-style numerical library for MoonBit. Its
`LapackMat` layer exposes LAPACK-backed `fmat_solve`, `fmat_inv`,
`fmat_cholesky`, `fmat_qr`, `fmat_svd`, `fmat_det`, and `fmat_lstsq` through a
native C FFI. On macOS it uses Apple Accelerate; on Linux it expects OpenBLAS
and LAPACK.

That means a dense, native `Float` solve is already available in the MoonBit
ecosystem. `numsolve` must not pretend otherwise. The two libraries overlap on
the basic question "solve `A · x = b`"; this document explains where
`numsolve` extends that capability instead of merely repeating it.

## Different axes

| Axis | `mizchi/numbt` | `numsolve` |
| --- | --- | --- |
| Implementation | C FFI to BLAS/LAPACK | Pure MoonBit |
| State type | `Float` views or `LapackMat` byte storage | `Array[Array[Double]]` |
| Target | Native, with platform BLAS dependencies | Portable MoonBit, designed for Wasm and native |
| Solve granularity | One-shot in-place `fmat_solve` | Reusable `LU` / `Cholesky` factor objects |
| Multiple right-hand sides | Caller loops over a solve | `solve_many` keeps the factor and answers a batch |
| Transposed system | Caller must arrange storage or refactor | `LU::solve_transposed` reuses the same factors |
| Residual diagnostics | Not exposed at the pure MoonBit layer | `residual_norm` uses a compensated inner product |
| Iterative refinement | Not exposed | `LU::solve_refined` reuses the factor and a compensated residual |
| Iterative solvers | Direct LAPACK methods | `jacobi`, `gauss_seidel`, and conjugate-gradient methods |
| Structured tridiagonal systems | Caller uses a dense solve | `solve_tridiagonal` is O(n) and is reused by the spline module |
| Error reporting | LAPACK status translated to `Bool` / option values | Explicit `LinError` variants with the failing index, size, tolerance, or residual |

The boundary is therefore not "numbt versus numsolve" but "accelerated native
dense linear algebra versus a portable pure-MoonBit numerical toolkit." A
project that can use Accelerate or OpenBLAS should prefer `numbt` for large
dense native matrices. A Wasm deployment, a project without a native
toolchain, or a model built around reusable factors, tridiagonal structure,
iterative methods, or explicit numerical diagnostics needs the capabilities
that `numsolve` adds.

## What `numsolve` adds

### Reusable factorizations

`LU::decompose` and `Cholesky::decompose` return factor objects. A caller can
solve `A · x = b`, `A · X = B`, `Aᵀ · x = b`, or form the inverse from one
factorization. The LAPACK layer is optimized for a single in-place call; the
pure-MoonBit layer makes the factorization itself a value that can be reused
and inspected through its dimensions and determinant.

### Compensated residual and refinement

The residual `b - A · x` is computed with an error-free product and a
compensated sum. Without that, an ill-conditioned system can have a residual
that is numerically zero even when the forward error is not. `solve_refined`
turns that residual into a correction equation and solves it with the
existing factor. It reduces the residual without refactoring, while the
documentation still states the limit: forward error on a very ill-conditioned
matrix is bounded by the condition number and by double precision.

### Iterative solvers

`jacobi` and `gauss_seidel` add a stationary-iteration family. They need no
factorization, use only matrix-vector products, and stop when the residual
meets a tolerance. They are appropriate when the matrix is diagonally
dominant or large and the caller wants an iterative method; they raise
`ZeroDiagonal` or `NoConvergence` instead of returning a plausible-looking
last iterate.

`conjugate_gradient` adds a Krylov method for symmetric positive definite
systems, and `preconditioned_conjugate_gradient` uses the diagonal as a
Jacobi preconditioner. These methods choose search directions from the
residual history rather than updating one component at a time. A non-positive
curvature is reported as `Breakdown`, which distinguishes an indefinite or
non-symmetric matrix from a run that merely needs more iterations.

### Tridiagonal structure

`solve_tridiagonal` implements the Thomas algorithm in O(n), and the cubic
spline module uses it directly. A dense LAPACK solve can solve the same
system, but it does not preserve the structure or the linear complexity, and
it cannot be used by the Wasm build of the spline module without the native
dependency.

## Interoperation, not replacement

The libraries can coexist in one workspace:

- use `numbt` when the data is already `Float`/`LapackMat`, the target is
  native, and the hot path is a large dense BLAS/LAPACK operation;
- use `numsolve` when the target must be portable, when the state is
  `Double`, when the system has tridiagonal or iteratively solved structure,
  or when the caller needs reusable factors, explicit errors, and residual
  diagnostics;
- convert at the boundary deliberately. `numbt` stores LAPACK matrices in
  byte-backed storage for zero-copy native calls; `numsolve` stores rows as
  `Array[Double]`. The conversion is an explicit data-layout decision, not a
  hidden dependency.

`numsolve` does not vendor, wrap, or reimplement LAPACK. It provides a
smaller, auditable pure-MoonBit kernel with a different set of guarantees and
a documented extension path.

## Verification

The linear tests cover:

- partial pivoting against an elimination without row exchanges;
- reusable `solve_many`, transpose solve, determinant, and inverse;
- Cholesky, tridiagonal, residual, and refinement behavior;
- Jacobi and Gauss-Seidel on a diagonally dominant system;
- conjugate gradient and Jacobi-preconditioned conjugate gradient on an SPD system;
- zero diagonal, tolerance, iteration-limit, and no-convergence errors;
- Hilbert 8 as a deliberately ill-conditioned system;
- empty and one-row systems as boundary cases.
