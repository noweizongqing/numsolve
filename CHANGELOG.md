# Changelog

## 0.6.0

- Add `VectorTrajectory::truncate_at` for terminal event trajectories.
- Preserve the derivative at the interpolated truncation point.
- Add terminal workflow tests, installation instructions, contribution guide,
  changelog, and license documentation.

## 0.5.1

- Sync generated interfaces with current MoonBit `Int32` aliases.
- Restore GitHub Actions compatibility with the current formatter.

## 0.5.0

- Add cubic Hermite dense output (`interpolate_at`, `sample_at`).
- Add first-event detection with direction filtering (`first_event`).

## 0.4.1

- Migrate source formatting to the current MoonBit formatter.
- Restore GitHub Actions on the latest toolchain.

## 0.4.0

- Add conjugate gradient and Jacobi-preconditioned conjugate gradient.
- Add general banded elimination.
- Add transpose and multiple-right-hand-side solves.

## 0.3.0

- Add compensated residuals and iterative refinement.
- Add Jacobi and Gauss-Seidel iterations.
- Document the relationship to `mizchi/numbt`.

## 0.2.0

- Add vector ODE systems and the `solve_ivp` entry point.

## 0.1.0

- Add root finding, scalar ODEs, linear systems, interpolation, CLI, tests, and
  benchmarks.
