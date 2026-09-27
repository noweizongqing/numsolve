///|
/// Numerical solvers in pure MoonBit: root finding, ordinary differential
/// equations, linear systems, and interpolation.
///
/// * root finding (`roots.mbt`) - where a function of one variable is zero:
///   `bisection`, `false_position`, `secant`, `newton` and `brent`, each
///   asked for a `tolerance` and a `max_iterations`, and each reporting what
///   went wrong — a bracket with no sign change, a slope of zero, a negative
///   tolerance, a search that ran out of iterations — as a `RootError`
///   carrying the values involved.
/// * ordinary differential equations (`ode.mbt`) - `euler`, `midpoint` and
///   `rk4` over a fixed number of steps, and `rk45`, the Dormand-Prince 5(4)
///   stepper that chooses its own step size from the error it estimates.
/// * linear systems (`linear.mbt`) - an `LU` factorization with partial
///   pivoting, a `Cholesky` one for symmetric positive definite matrices, and
///   the Thomas recurrence for a tridiagonal system.
/// * interpolation (`interp.mbt`) - a line and a cubic spline through a table
///   of points, the spline under a natural, clamped or not-a-knot boundary.
///
/// Every module is written against a published algorithm and checked against
/// a reference implementation of it, and the tests measure the convergence
/// the algorithm claims rather than only its final answer: bisection's error
/// against `(b - a) / 2^k`, Newton's against its own square, the secant's
/// against the golden ratio's power.
///
/// ```moonbit nocheck
/// ///|
/// let root = @numsolve.brent(fn(x) { x * x - 2.0 }, 1.0, 2.0)
/// ```
