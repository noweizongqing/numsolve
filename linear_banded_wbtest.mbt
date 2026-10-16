///|
test "a banded solve agrees with the tridiagonal solver" {
  let a : Array[Array[Double]] = [
    [4.0, 1.0, 0.0],
    [2.0, 3.0, 1.0],
    [0.0, 1.0, 2.0],
  ]
  let rhs : Array[Double] = [6.0, 11.0, 8.0]
  let expected : Array[Double] = [1.0, 2.0, 3.0]
  assert_true(
    linear_close_vector(solve_banded(1, 1, a, rhs), expected, 1.0e-12),
  )
  assert_true(
    linear_close_vector(
      solve_banded(1, 1, a, rhs),
      solve_tridiagonal([2.0, 1.0], [4.0, 3.0, 2.0], [1.0, 1.0], rhs),
      1.0e-12,
    ),
  )
}

///|
test "a general banded system keeps its declared structure" {
  let a : Array[Array[Double]] = [
    [4.0, 1.0, 2.0, 0.0],
    [1.0, 3.0, 1.0, 2.0],
    [0.0, 1.0, 3.0, 1.0],
    [0.0, 0.0, 1.0, 2.0],
  ]
  let rhs : Array[Double] = [7.0, 7.0, 5.0, 3.0]
  assert_true(
    linear_close_vector(
      solve_banded(1, 2, a, rhs),
      [1.0, 1.0, 1.0, 1.0],
      1.0e-12,
    ),
  )
}

///|
test "a matrix outside the declared band is refused" {
  let a : Array[Array[Double]] = [
    [4.0, 1.0, 0.0],
    [2.0, 3.0, 1.0],
    [1.0, 0.0, 2.0],
  ]
  let failed = try {
    let _ = solve_banded(1, 1, a, [1.0, 2.0, 3.0])
    false
  } catch {
    LinError::NotBanded(row, column) => row == 2 && column == 0
    _ => false
  }
  assert_true(failed)
}

///|
test "a banded solve validates bandwidths and dimensions" {
  let failed = try {
    let _ = solve_banded(-1, 0, [[1.0]], [1.0])
    false
  } catch {
    LinError::InvalidBandwidth(width) => width == -1
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = solve_banded(1, 1, [[1.0, 0.0], [0.0, 1.0]], [1.0])
    false
  } catch {
    LinError::SizeMismatch(needed, given) => needed == 2 && given == 1
    _ => false
  }
  assert_true(failed)
}

///|
test "an empty banded system is already solved" {
  assert_eq(solve_banded(1, 1, [], []).length(), 0)
}
