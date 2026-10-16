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
fn iterative_zeros(size : Int32) -> Array[Double] {
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
  max_iterations : Int32,
  use_latest : Bool,
) -> (Array[Double], Int32) raise LinError {
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
  tolerance~ : Double = 1.0e-10,
  max_iterations~ : Int32 = 1000,
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
  tolerance~ : Double = 1.0e-10,
  max_iterations~ : Int32 = 1000,
) -> Array[Double] raise LinError {
  stationary_solve(a, b, tolerance, max_iterations, true).0
}

///|
/// The matrix-vector product used by the iterative methods.
fn iterative_matrix_vector(
  a : Array[Array[Double]],
  x : Array[Double],
) -> Array[Double] raise LinError {
  let size = checked_square(a)
  let _ = checked_size(size, x)
  let product : Array[Double] = Array(capacity=size)
  let mut i = 0
  while i < size {
    let (value, error) = compensated_dot(a[i], x)
    product.push(value + error)
    i = i + 1
  }
  product
}

///|
/// The Euclidean norm of a vector.
fn iterative_norm(values : Array[Double]) -> Double {
  compressed_dot_self(values).sqrt()
}

///|
/// A small wrapper around the compensated dot product for a vector with
/// itself, kept separate so the norm and the CG loops cannot drift apart.
fn compressed_dot_self(values : Array[Double]) -> Double {
  let (sum, error) = compensated_dot(values, values)
  sum + error
}

///|
/// Add `scale * direction` to `target` in place.
fn iterative_add_scaled(
  target : Array[Double],
  scale : Double,
  direction : Array[Double],
) -> Unit {
  let mut i = 0
  while i < target.length() {
    target[i] = target[i] + scale * direction[i]
    i = i + 1
  }
}

///|
/// Solve a symmetric positive definite system with the conjugate-gradient
/// method.
///
/// Each pass uses one matrix-vector product and two inner products, and the
/// residual decreases in the A-norm. A non-positive `pᵀ·A·p` means the
/// matrix is not positive definite or the iteration has broken down, and is
/// raised as `Breakdown` rather than silently continued.
pub fn conjugate_gradient(
  a : Array[Array[Double]],
  b : Array[Double],
  tolerance~ : Double = 1.0e-10,
  max_iterations~ : Int32 = 1000,
) -> Array[Double] raise LinError {
  if tolerance <= 0.0 {
    raise LinError::NonPositiveTolerance(tolerance)
  }
  if max_iterations <= 0 {
    raise LinError::NonPositiveMaxIterations(max_iterations)
  }
  let size = checked_square(a)
  let _ = checked_size(size, b)
  if size == 0 {
    return []
  }
  let x = iterative_zeros(size)
  let residual = b.copy()
  if iterative_norm(residual) <= tolerance {
    return x
  }
  let mut direction = residual.copy()
  let mut residual_dot = compressed_dot_self(residual)
  let mut iteration = 0
  while iteration < max_iterations {
    let product = iterative_matrix_vector(a, direction)
    let curvature = {
      let (sum, error) = compensated_dot(direction, product)
      sum + error
    }
    if curvature <= 0.0 {
      raise LinError::Breakdown(iteration)
    }
    let alpha = residual_dot / curvature
    iterative_add_scaled(x, alpha, direction)
    iterative_add_scaled(residual, -alpha, product)
    let next_dot = compressed_dot_self(residual)
    iteration = iteration + 1
    if next_dot.sqrt() <= tolerance {
      return x
    }
    let beta = next_dot / residual_dot
    let next_direction = Array(capacity=size)
    let mut i = 0
    while i < size {
      next_direction.push(residual[i] + beta * direction[i])
      i = i + 1
    }
    direction = next_direction
    residual_dot = next_dot
  }
  raise LinError::NoConvergence(iteration, residual_norm(a, b, x))
}

///|
/// Solve a symmetric positive definite system with Jacobi-preconditioned
/// conjugate gradient. The preconditioner is the diagonal of `A`, which
/// costs one division per component and often reduces the iteration count of
/// systems whose diagonal varies by orders of magnitude.
pub fn preconditioned_conjugate_gradient(
  a : Array[Array[Double]],
  b : Array[Double],
  tolerance~ : Double = 1.0e-10,
  max_iterations~ : Int32 = 1000,
) -> Array[Double] raise LinError {
  if tolerance <= 0.0 {
    raise LinError::NonPositiveTolerance(tolerance)
  }
  if max_iterations <= 0 {
    raise LinError::NonPositiveMaxIterations(max_iterations)
  }
  let size = checked_square(a)
  let _ = checked_size(size, b)
  if size == 0 {
    return []
  }
  let diagonal : Array[Double] = Array(capacity=size)
  let mut i = 0
  while i < size {
    if a[i][i] == 0.0 {
      raise LinError::ZeroDiagonal(i)
    }
    diagonal.push(a[i][i])
    i = i + 1
  }
  let x = iterative_zeros(size)
  let residual = b.copy()
  let mut z : Array[Double] = Array(capacity=size)
  i = 0
  while i < size {
    z.push(residual[i] / diagonal[i])
    i = i + 1
  }
  let mut direction = z.copy()
  let mut rho = {
    let (sum, error) = compensated_dot(residual, z)
    sum + error
  }
  if rho.sqrt() <= tolerance {
    return x
  }
  let mut iteration = 0
  while iteration < max_iterations {
    let product = iterative_matrix_vector(a, direction)
    let curvature = {
      let (sum, error) = compensated_dot(direction, product)
      sum + error
    }
    if curvature <= 0.0 {
      raise LinError::Breakdown(iteration)
    }
    let alpha = rho / curvature
    iterative_add_scaled(x, alpha, direction)
    iterative_add_scaled(residual, -alpha, product)
    z = Array(capacity=size)
    i = 0
    while i < size {
      z.push(residual[i] / diagonal[i])
      i = i + 1
    }
    let next_rho = {
      let (sum, error) = compensated_dot(residual, z)
      sum + error
    }
    iteration = iteration + 1
    if iterative_norm(residual) <= tolerance {
      return x
    }
    let beta = next_rho / rho
    let next_direction = Array(capacity=size)
    i = 0
    while i < size {
      next_direction.push(z[i] + beta * direction[i])
      i = i + 1
    }
    direction = next_direction
    rho = next_rho
  }
  raise LinError::NoConvergence(iteration, residual_norm(a, b, x))
}
