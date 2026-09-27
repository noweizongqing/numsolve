///|
/// Whether two values agree to within `tolerance`.
fn linear_close(actual : Double, expected : Double, tolerance : Double) -> Bool {
  (actual - expected).abs() <= tolerance
}

///|
/// Whether two vectors agree entry by entry.
fn linear_close_vector(
  actual : Array[Double],
  expected : Array[Double],
  tolerance : Double,
) -> Bool {
  if actual.length() != expected.length() {
    return false
  }
  let mut i = 0
  while i < actual.length() {
    if !linear_close(actual[i], expected[i], tolerance) {
      return false
    }
    i = i + 1
  }
  true
}

///|
/// Gaussian elimination with no row exchanges, for the test that measures what
/// pivoting is worth: the same arithmetic as `LU::decompose` with the search
/// for a pivot taken out, and the same failure when the diagonal entry it is
/// left with is zero.
fn linear_no_pivot_solve(
  a : Array[Array[Double]],
  b : Array[Double],
) -> Array[Double] raise LinError {
  let size = a.length()
  let f : Array[Array[Double]] = Array(capacity=size)
  let mut i = 0
  while i < size {
    let row = a[i].copy()
    row.push(b[i])
    f.push(row)
    i = i + 1
  }
  let mut column = 0
  while column < size {
    if f[column][column] == 0.0 {
      raise LinError::Singular(column)
    }
    let mut row = column + 1
    while row < size {
      let multiplier = f[row][column] / f[column][column]
      let mut j = column
      while j <= size {
        f[row][j] = f[row][j] - multiplier * f[column][j]
        j = j + 1
      }
      row = row + 1
    }
    column = column + 1
  }
  let x : Array[Double] = Array(capacity=size)
  let mut k = 0
  while k < size {
    x.push(0.0)
    k = k + 1
  }
  let mut row = size - 1
  while row >= 0 {
    let mut acc = f[row][size]
    let mut j = row + 1
    while j < size {
      acc = acc - f[row][j] * x[j]
      j = j + 1
    }
    x[row] = acc / f[row][row]
    row = row - 1
  }
  x
}

///|
/// The matrix the factorizations are shown on, and its right-hand side. Its
/// exact solution is [0.5, 0, 0] and its determinant is 4, so both the answer
/// and the determinant can be compared against numbers rather than against
/// another implementation.
fn linear_sample_matrix() -> Array[Array[Double]] {
  [[2.0, 1.0, 1.0], [4.0, 3.0, 3.0], [8.0, 7.0, 9.0]]
}

///|
/// The LU factors of the sample matrix, against values computed from the same
/// elimination in NumPy: the largest element of the first column is exchanged
/// into the pivot, which is why the factors are those of the permuted matrix
/// and why the permutation is reported alongside them.
test "an LU decomposition is the elimination it names" {
  let a = linear_sample_matrix()
  let decomposition = LU::decompose(a)
  assert_eq(decomposition.size(), 3)
  assert_eq(decomposition.swaps, 2)
  assert_eq(decomposition.pivots[0], 2)
  assert_eq(decomposition.pivots[1], 0)
  assert_eq(decomposition.pivots[2], 1)
  let factors = decomposition.factors
  let expected : Array[Array[Double]] = [
    [8.0, 7.0, 9.0],
    [0.25, -0.75, -1.25],
    [0.5, 0.6666666666666666, -0.6666666666666667],
  ]
  let mut i = 0
  while i < 3 {
    assert_true(linear_close_vector(factors[i], expected[i], 1.0e-15))
    i = i + 1
  }
}

///|
/// The matrix that was factored is not written over: the factors are built in
/// a copy, so a caller can keep the matrix it handed in.
test "a factorization leaves its argument alone" {
  let a = linear_sample_matrix()
  let _ = LU::decompose(a)
  assert_true(linear_close(a[0][0], 2.0, 0.0))
  assert_true(linear_close(a[2][2], 9.0, 0.0))
  let _ = Cholesky::decompose([[4.0, 2.0], [2.0, 3.0]])
  assert_true(linear_close(a[1][1], 3.0, 0.0))
}

///|
/// The solution, the determinant, and the inverse of the same system. The
/// inverse is checked twice: against the values NumPy produced for it, and by
/// multiplying it back into the matrix, which is the property that has to
/// hold — an inverse that is itself inverted, rather than one that is right.
test "an LU solve, its determinant and its inverse" {
  let decomposition = LU::decompose(linear_sample_matrix())
  assert_true(
    linear_close_vector(
      decomposition.solve([1.0, 2.0, 4.0]),
      [0.5, 0.0, 0.0],
      1.0e-15,
    ),
  )
  assert_true(linear_close(decomposition.determinant(), 4.0, 1.0e-15))
  let inverse = decomposition.inverse()
  let expected : Array[Array[Double]] = [
    [1.5000000000000002, -0.4999999999999998, 0.0],
    [-3.0, 2.4999999999999996, -0.5],
    [0.9999999999999999, -1.4999999999999998, 0.5],
  ]
  let mut i = 0
  while i < 3 {
    assert_true(linear_close_vector(inverse[i], expected[i], 1.0e-12))
    i = i + 1
  }
  let a = linear_sample_matrix()
  let mut row = 0
  while row < 3 {
    let mut column = 0
    while column < 3 {
      let mut acc = 0.0
      let mut k = 0
      while k < 3 {
        acc = acc + a[row][k] * inverse[k][column]
        k = k + 1
      }
      assert_true(
        linear_close(acc, if row == column { 1.0 } else { 0.0 }, 1.0e-12),
      )
      column = column + 1
    }
    row = row + 1
  }
}

///|
/// `solve_many` is the same answer as one `solve` per right-hand side, in one
/// call: the factorization is paid for once.
test "a solve with many right-hand sides" {
  let decomposition = LU::decompose(linear_sample_matrix())
  let left : Array[Array[Double]] = [
    [1.0, 2.0, 4.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
  ]
  let solutions = decomposition.solve_many(left)
  assert_eq(solutions.length(), 3)
  let mut i = 0
  while i < 3 {
    assert_true(
      linear_close_vector(solutions[i], decomposition.solve(left[i]), 0.0),
    )
    i = i + 1
  }
  assert_true(
    linear_close_vector(
      solutions[1],
      [1.5000000000000002, -3.0, 0.9999999999999999],
      1.0e-15,
    ),
  )
}

///|
/// The Cholesky factorization of a symmetric positive definite matrix against
/// SciPy's, and the property that defines it: `L · Lᵀ` is the matrix it was
/// built from, to within the rounding of the multiplications.
test "a Cholesky factorization is the matrix times its transpose" {
  let m : Array[Array[Double]] = [
    [4.0, 2.0, 1.0],
    [2.0, 3.0, 1.0],
    [1.0, 1.0, 2.0],
  ]
  let decomposition = Cholesky::decompose(m)
  assert_eq(decomposition.size(), 3)
  let lower = decomposition.lower
  let expected : Array[Array[Double]] = [
    [2.0, 0.0, 0.0],
    [1.0, 1.4142135623730951, 0.0],
    [0.5, 0.35355339059327373, 1.2747548783981961],
  ]
  let mut i = 0
  while i < 3 {
    assert_true(linear_close_vector(lower[i], expected[i], 1.0e-15))
    i = i + 1
  }
  let mut row = 0
  while row < 3 {
    let mut column = 0
    while column < 3 {
      let mut acc = 0.0
      let mut k = 0
      while k < 3 {
        acc = acc + lower[row][k] * lower[column][k]
        k = k + 1
      }
      assert_true(linear_close(acc, m[row][column], 1.0e-12))
      column = column + 1
    }
    row = row + 1
  }
}

///|
/// A Cholesky solve and the determinant, which comes out as the square of the
/// product of the diagonal: 13 for this matrix, exactly.
test "a Cholesky solve and its determinant" {
  let m : Array[Array[Double]] = [
    [4.0, 2.0, 1.0],
    [2.0, 3.0, 1.0],
    [1.0, 1.0, 2.0],
  ]
  let decomposition = Cholesky::decompose(m)
  assert_true(
    linear_close_vector(
      decomposition.solve([1.0, 4.0, 9.0]),
      [-1.2307692307692308, 0.5384615384615384, 4.846153846153847],
      1.0e-12,
    ),
  )
  assert_true(linear_close(decomposition.determinant(), 13.0, 1.0e-12))
}

///|
/// A matrix that is symmetric but not positive definite has a pivot that is
/// not positive, and Cholesky has no square root to take of it: the second
/// pivot of [[1, 2], [2, 1]] is −3, and the error says which one it was.
test "a matrix that is not positive definite is refused" {
  let failed = try {
    let _ = Cholesky::decompose([[1.0, 2.0], [2.0, 1.0]])
    false
  } catch {
    LinError::NotPositiveDefinite(index) => index == 1
    _ => false
  }
  assert_true(failed)
}

///|
/// The tridiagonal solver on the system whose three diagonals are
/// [1, 2], [3, 4, 5] and [1, 2], against the values NumPy gives for the same
/// matrix written out in full.
test "a tridiagonal system is solved in one pass" {
  let solution = solve_tridiagonal([1.0, 2.0], [3.0, 4.0, 5.0], [1.0, 2.0], [
    1.0, 2.0, 3.0,
  ])
  assert_true(
    linear_close_vector(
      solution,
      [0.27906976744186046, 0.16279069767441862, 0.5348837209302325],
      1.0e-15,
    ),
  )
}

///|
/// The sizes a tridiagonal system has to have: three diagonals of `n`, `n − 1`
/// and `n − 1`, and a right-hand side of `n`. Anything else is refused with
/// the length that was needed and the length that was given.
test "the tridiagonal solver checks the shape of its bands" {
  let failed = try {
    let _ = solve_tridiagonal([1.0], [3.0, 4.0, 5.0], [1.0, 2.0], [
      1.0, 2.0, 3.0,
    ])
    false
  } catch {
    LinError::SizeMismatch(needed, given) => needed == 2 && given == 1
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = solve_tridiagonal([1.0, 2.0], [3.0, 4.0, 5.0], [1.0, 2.0], [
      1.0, 2.0,
    ])
    false
  } catch {
    LinError::SizeMismatch(needed, given) => needed == 3 && given == 2
    _ => false
  }
  assert_true(failed)
}

///|
/// A zero pivot inside the band leaves nothing to exchange it with: the
/// tridiagonal solver refuses rather than dividing by it, and says which
/// diagonal entry it was.
test "the tridiagonal solver refuses a zero pivot" {
  let failed = try {
    let _ = solve_tridiagonal([1.0], [0.0, 0.0], [1.0], [1.0, 2.0])
    false
  } catch {
    LinError::Singular(index) => index == 0
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = solve_tridiagonal([1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0], [
      1.0, 2.0, 3.0,
    ])
    false
  } catch {
    LinError::Singular(index) => index == 1
    _ => false
  }
  assert_true(failed)
}

///|
/// A matrix that is not square, and a right-hand side of the wrong length, are
/// both refused with the two numbers that disagree.
test "the shape of a system is checked before it is solved" {
  let failed = try {
    let _ = LU::decompose([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    false
  } catch {
    LinError::NotSquare(rows, columns) => rows == 2 && columns == 3
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = LU::decompose([[1.0, 2.0], [3.0]])
    false
  } catch {
    LinError::NotSquare(rows, columns) => rows == 2 && columns == 1
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = solve_linear_system(linear_sample_matrix(), [1.0, 2.0])
    false
  } catch {
    LinError::SizeMismatch(needed, given) => needed == 3 && given == 2
    _ => false
  }
  assert_true(failed)
}

///|
/// A singular matrix has no pivot to eliminate with, whatever row is exchanged
/// into place, and the elimination stops at the column where it found that
/// out rather than dividing by zero.
test "a singular matrix is refused" {
  let failed = try {
    let _ = LU::decompose([[1.0, 2.0], [2.0, 4.0]])
    false
  } catch {
    LinError::Singular(index) => index == 1
    _ => false
  }
  assert_true(failed)
}

///|
/// Pivoting decides the answer, not just its accuracy. The first matrix has a
/// zero where the elimination would divide, so without row exchanges there is
/// no answer at all; the second has a pivot so small that the multipliers
/// below it are 1e18, and the answer that comes back without exchanging rows
/// is [0, 1] for an exact solution of [1, 1].
test "pivoting is what makes the answer right" {
  let zero_pivot : Array[Array[Double]] = [[0.0, 2.0], [3.0, 4.0]]
  assert_true(
    linear_close_vector(
      solve_linear_system(zero_pivot, [2.0, 7.0]),
      [1.0, 1.0],
      1.0e-15,
    ),
  )
  let failed = try {
    let _ = linear_no_pivot_solve(zero_pivot, [2.0, 7.0])
    false
  } catch {
    LinError::Singular(index) => index == 0
    _ => false
  }
  assert_true(failed)
  let tiny_pivot : Array[Array[Double]] = [[1.0e-18, 1.0], [1.0, 1.0]]
  assert_true(
    linear_close_vector(
      solve_linear_system(tiny_pivot, [1.0, 2.0]),
      [1.0, 1.0],
      1.0e-12,
    ),
  )
  assert_true(
    linear_close_vector(
      linear_no_pivot_solve(tiny_pivot, [1.0, 2.0]),
      [0.0, 1.0],
      1.0e-12,
    ),
  )
}

///|
/// The Hilbert matrix of order 8 — `H[i][j] = 1 / (i + j + 1)`, the textbook
/// ill-conditioned matrix, with a condition number of 1.5e10 — solved for the
/// right-hand side made of its row sums, whose exact solution is all ones.
///
/// The error is 1.4e-07, and the point of the test is that it is that and not
/// worse: at this condition number about ten of the sixteen digits a double
/// has are left, and an elimination that did not pivot would spend the rest
/// of them. The bound asserted is 1e-5, wide enough that a different order of
/// summation inside the same algorithm — LAPACK reaches 2.5e-07 here, NumPy
/// 1.3e-07 — cannot fail it, and narrow enough that losing a digit or two to a
/// change in the pivoting would.
test "a Hilbert system is solved to the accuracy it has left" {
  let size = 8
  let hilbert : Array[Array[Double]] = Array(capacity=size)
  let rhs : Array[Double] = Array(capacity=size)
  let mut i = 0
  while i < size {
    let row : Array[Double] = Array(capacity=size)
    let mut j = 0
    while j < size {
      row.push(1.0 / (i + j + 1).to_double())
      j = j + 1
    }
    hilbert.push(row)
    rhs.push(row.fold(init=0.0, fn(acc, value) { acc + value }))
    i = i + 1
  }
  let solution = solve_linear_system(hilbert, rhs)
  let mut worst = 0.0
  let mut k = 0
  while k < size {
    let error = (solution[k] - 1.0).abs()
    if error > worst {
      worst = error
    }
    k = k + 1
  }
  assert_true(worst <= 1.0e-5)
}

///|
/// A system of no rows is square, consistent and already solved: there is
/// nothing to factor, the determinant of an empty product is one, and solving
/// it returns nothing rather than failing.
test "a system of no rows is solved by doing nothing" {
  let decomposition = LU::decompose([])
  assert_eq(decomposition.size(), 0)
  assert_true(linear_close(decomposition.determinant(), 1.0, 0.0))
  assert_eq(decomposition.solve([]).length(), 0)
  let cholesky = Cholesky::decompose([])
  assert_eq(cholesky.size(), 0)
  assert_true(linear_close(cholesky.determinant(), 1.0, 0.0))
  assert_eq(solve_tridiagonal([], [], [], []).length(), 0)
}

///|
/// A one-by-one system goes through both factorizations: the LU of [[4]] is
/// itself, the Cholesky of it is [[2]], and a tridiagonal system of one row is
/// a division.
test "a system of one row is a division" {
  let lu = LU::decompose([[4.0]])
  assert_true(linear_close(lu.determinant(), 4.0, 0.0))
  assert_true(linear_close_vector(lu.solve([8.0]), [2.0], 0.0))
  let cholesky = Cholesky::decompose([[4.0]])
  assert_true(linear_close(cholesky.lower[0][0], 2.0, 0.0))
  assert_true(linear_close(cholesky.determinant(), 4.0, 0.0))
  assert_true(linear_close_vector(cholesky.solve([8.0]), [2.0], 0.0))
  assert_true(
    linear_close_vector(solve_tridiagonal([], [4.0], [], [8.0]), [2.0], 0.0),
  )
}
