///|
/// The diagonally dominant system used by the iterative tests. Its exact
/// solution is [1, 2, 3].
fn iterative_sample_matrix() -> Array[Array[Double]] {
  [[4.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 2.0]]
}

///|
/// The right-hand side of the sample matrix for the exact solution [1, 2, 3].
fn iterative_sample_rhs() -> Array[Double] {
  [6.0, 10.0, 8.0]
}

///|
test "Jacobi and Gauss-Seidel solve a diagonally dominant system" {
  let a = iterative_sample_matrix()
  let b = iterative_sample_rhs()
  let expected : Array[Double] = [1.0, 2.0, 3.0]
  let jacobi_solution = jacobi(a, b, tolerance=1.0e-12, max_iterations=200)
  let gauss_seidel_solution = gauss_seidel(
    a,
    b,
    tolerance=1.0e-12,
    max_iterations=200,
  )
  assert_true(linear_close_vector(jacobi_solution, expected, 1.0e-9))
  assert_true(linear_close_vector(gauss_seidel_solution, expected, 1.0e-9))
}

///|
test "Gauss-Seidel takes no more iterations than Jacobi on the same system" {
  let a = iterative_sample_matrix()
  let b = iterative_sample_rhs()
  let (jacobi_solution, jacobi_iterations) = stationary_solve(
    a, b, 1.0e-12, 200, false,
  )
  let (gauss_seidel_solution, gauss_seidel_iterations) = stationary_solve(
    a, b, 1.0e-12, 200, true,
  )
  assert_true(gauss_seidel_iterations <= jacobi_iterations)
  assert_true(residual_norm(a, b, jacobi_solution) <= 1.0e-12)
  assert_true(residual_norm(a, b, gauss_seidel_solution) <= 1.0e-12)
}

///|
test "an iterative solve with a zero diagonal is refused" {
  let failed = try {
    let _ = jacobi([[0.0, 1.0], [1.0, 1.0]], [1.0, 1.0], max_iterations=10)
    false
  } catch {
    LinError::ZeroDiagonal(row) => row == 0
    _ => false
  }
  assert_true(failed)
}

///|
test "an iterative solve that cannot converge reports its budget" {
  let failed = try {
    let _ = gauss_seidel(
      [[1.0, 2.0], [2.0, 1.0]],
      [1.0, 1.0],
      tolerance=1.0e-12,
      max_iterations=20,
    )
    false
  } catch {
    LinError::NoConvergence(iterations, residual) =>
      iterations == 20 && residual > 1.0e-12
    _ => false
  }
  assert_true(failed)
}

///|
test "iterative solvers validate tolerance and iteration limits" {
  let a = iterative_sample_matrix()
  let b = iterative_sample_rhs()
  let failed = try {
    let _ = jacobi(a, b, tolerance=0.0)
    false
  } catch {
    LinError::NonPositiveTolerance(value) => value == 0.0
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = gauss_seidel(a, b, max_iterations=0)
    false
  } catch {
    LinError::NonPositiveMaxIterations(limit) => limit == 0
    _ => false
  }
  assert_true(failed)
}

///|
test "an empty iterative system is already solved" {
  assert_eq(jacobi([], []).length(), 0)
  assert_eq(gauss_seidel([], []).length(), 0)
}
