///|
/// Stationary iterative solvers for diagonally dominant linear systems:
/// Jacobi and Gauss-Seidel.
///
/// These methods complement the direct factorizations in `linear.mbt`. They
/// need no factorization and use only matrix-vector products, but they
/// converge only for suitable matrices and need an iteration budget. They are
/// deliberately separate from the dense direct solve so that a caller has to
/// choose the method and its tolerance explicitly.

///|
/// A vector of `size` zeros.
fn iterative_zeros(size : Int) -> Array[Double] {
  let values : Array[Double] = Array(capacity=size)
  let mut i = 0
  while i < size {
    values.push(0.0)
    i = i + 1
  }
  values
}

///|
/// Shared Jacobi / Gauss-Seidel loop, returning the solution and the number of
/// iterations taken. `use_latest` is true for Gauss-Seidel and false for
/// Jacobi.
fn stationary_solve(
  a : Array[Array[Double]],
  b : Array[Double],
  tolerance : Double,
  max_iterations : Int,
  use_latest : Bool,
) -> (Array[Double], Int) raise LinError {
  if tolerance <= 0.0 {
    raise LinError::NonPositiveTolerance(tolerance)
  }
  if max_iterations <= 0 {
    raise LinError::NonPositiveMaxIterations(max_iterations)
  }
  let size = checked_square(a)
  let _ = checked_size(size, b)
  if size == 0 {
    return ([], 0)
  }
  let mut i = 0
  while i < size {
    if a[i][i] == 0.0 {
      raise LinError::ZeroDiagonal(i)
    }
    i = i + 1
  }
  let x = iterative_zeros(size)
  if residual_norm(a, b, x) <= tolerance {
    return (x, 0)
  }
  let mut iteration = 0
  let mut converged = false
  while iteration < max_iterations && !converged {
    if use_latest {
      let mut row = 0
      while row < size {
        let mut acc = b[row]
        let mut column = 0
        while column < size {
          if column != row {
            acc = acc - a[row][column] * x[column]
          }
          column = column + 1
        }
        x[row] = acc / a[row][row]
        row = row + 1
      }
    } else {
      let previous = x.copy()
      let mut row = 0
      while row < size {
        let mut acc = b[row]
        let mut column = 0
        while column < size {
          if column != row {
            acc = acc - a[row][column] * previous[column]
          }
          column = column + 1
        }
        x[row] = acc / a[row][row]
        row = row + 1
      }
    }
    iteration = iteration + 1
    if residual_norm(a, b, x) <= tolerance {
      converged = true
    }
  }
  if !converged {
    raise LinError::NoConvergence(iteration, residual_norm(a, b, x))
  }
  (x, iteration)
}

///|
/// Solve `A · x = b` by the Jacobi iteration: at each pass every component is
/// updated from the previous iterate, so the updates can be computed in any
/// order.
///
/// The method converges when the diagonal is non-zero and the spectral radius
/// of the iteration matrix is below one; a diagonally dominant matrix is the
/// common sufficient case. A run that exhausts `max_iterations` raises
/// `NoConvergence` rather than returning the last iterate as if it were a
/// solution.
pub fn jacobi(
  a : Array[Array[Double]],
  b : Array[Double],
  tolerance? : Double = 1.0e-10,
  max_iterations? : Int = 1000,
) -> Array[Double] raise LinError {
  stationary_solve(a, b, tolerance, max_iterations, false).0
}

///|
/// Solve `A · x = b` by the Gauss-Seidel iteration: each component uses the
/// values already updated in the current pass, so the method usually needs
/// fewer iterations than Jacobi on the same matrix.
///
/// It has the same diagonal, tolerance and iteration-budget contract as
/// `jacobi`, and it raises `NoConvergence` when the budget runs out.
pub fn gauss_seidel(
  a : Array[Array[Double]],
  b : Array[Double],
  tolerance? : Double = 1.0e-10,
  max_iterations? : Int = 1000,
) -> Array[Double] raise LinError {
  stationary_solve(a, b, tolerance, max_iterations, true).0
}
