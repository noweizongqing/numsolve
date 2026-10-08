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
/// * ordinary differential equations (`ode.mbt` and `ode_system.mbt`) -
///   scalar `euler`, `midpoint`, `rk4` and adaptive `rk45`, plus vector-state
///   counterparts and a common `solve_ivp` entry point. The system solvers
///   use absolute and relative component tolerances.
/// * linear systems (`linear.mbt`) - an `LU` factorization with partial
///   pivoting, a `Cholesky` one for symmetric positive definite matrices, and
///   the Thomas recurrence for a tridiagonal system.
/// * interpolation (`interp.mbt`) - a line and a cubic spline through a table
///   of points, the spline under a natural, clamped or not-a-knot boundary.
///
/// The algorithms are standard published methods. Tests compare selected
/// results with analytic solutions or SciPy and measure the convergence the
/// algorithms claim rather than only their final answer: bisection's error
/// against `(b - a) / 2^k`, Newton's against its own square, the secant's
/// against the golden ratio's power, and the ODE steppers against their
/// advertised orders.
///
/// ```moonbit nocheck
/// ///|
/// let root = @numsolve.brent(fn(x) { x * x - 2.0 }, 1.0, 2.0)
/// ```
