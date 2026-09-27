///|
/// Linear systems: `A · x = b` for a dense square matrix, and the one banded
/// case worth its own code.
///
/// A matrix is an `Array[Array[Double]]` of rows, and a right-hand side is an
/// `Array[Double]` of the same length as the matrix has rows.
///
/// `LU` is the factorization `P · A = L · U` — a lower triangular matrix with
/// ones on its diagonal, an upper triangular one, and the row exchanges that
/// partial pivoting asked for. It is built by Gaussian elimination with the
/// largest remaining element in the column chosen as the pivot at every step,
/// which is what keeps a small leading element from being divided through and
/// destroying the answer. `LU::solve` then does one forward and one back
/// substitution, and asking again with a different `b` costs only those two
/// passes.
///
///     | 2  1  1 |          | 8   7    9 |
///     | 4  3  3 |   ->     | 1/4 -3/4 -5/4 |   with rows permuted, two swaps
///     | 8  7  9 |          | 1/2  2/3 -2/3 |
///
/// `Cholesky` is the factorization for a symmetric positive definite matrix,
/// `A = L · Lᵀ`, at half the work of an LU and with no pivoting: the pivots of
/// such a matrix are positive, so the square roots are real. It is the one to
/// use when the matrix is known to be of that kind, and it says so with an
/// error when it is not.
///
/// `solve_tridiagonal` is the Thomas algorithm for a matrix that is zero
/// everywhere but three central diagonals, which it takes as three arrays: an
/// `O(n)` solve for the systems that come out of a cubic spline or a
/// finite-difference grid, and the reason `interp.mbt` can build a spline in
/// one pass.

///|
/// A system that cannot be set up, or that has no single solution.
pub suberror LinError {
  /// The matrix is not square, or its rows are not all the same length. The
  /// payloads are the number of rows and the length of the first row. A
  /// matrix of no rows is square, and is solved by doing nothing.
  NotSquare(Int, Int)
  /// An array is not as long as the system needs it to be. The payloads are
  /// the length that was needed and the length that was given: the number of
  /// rows of the matrix against the length of a right-hand side, or `n - 1`
  /// against the length of a sub- or super-diagonal.
  SizeMismatch(Int, Int)
  /// A pivot came out exactly zero, and no row below it could be exchanged
  /// for it either: the matrix is singular and the system has no unique
  /// answer. The payload is the index of the column the elimination was in
  /// when it found this, which is also the last row that could have had a
  /// pivot.
  Singular(Int)
  /// Cholesky was given a matrix that is not positive definite. The payload
  /// is the index of the diagonal entry whose pivot came out zero or
  /// negative — the first one that did.
  NotPositiveDefinite(Int)
}

///|
/// Whether every row of `a` is as long as the first one, and the matrix is
/// therefore rectangular. Returns the number of columns, or raises with the
/// number of rows and the length of the offending row.
fn checked_columns(a : Array[Array[Double]]) -> Int raise LinError {
  let rows = a.length()
  if rows == 0 {
    return 0
  }
  let columns = a[0].length()
  let mut i = 1
  while i < rows {
    if a[i].length() != columns {
      raise LinError::NotSquare(rows, a[i].length())
    }
    i = i + 1
  }
  columns
}

///|
/// Whether `a` is square, and how wide it is.
fn checked_square(a : Array[Array[Double]]) -> Int raise LinError {
  let rows = a.length()
  let columns = checked_columns(a)
  if rows != columns {
    raise LinError::NotSquare(rows, columns)
  }
  columns
}

///|
/// A `rows` by `columns` matrix of zeros.
fn zero_matrix(rows : Int, columns : Int) -> Array[Array[Double]] {
  let matrix : Array[Array[Double]] = Array(capacity=rows)
  let mut i = 0
  while i < rows {
    let row : Array[Double] = Array(capacity=columns)
    let mut j = 0
    while j < columns {
      row.push(0.0)
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  matrix
}

///|
/// The rows of `a`, copied: a factorization writes over its working matrix,
/// and the matrix the caller handed in is not its to write over.
fn copy_rows(a : Array[Array[Double]]) -> Array[Array[Double]] {
  let rows = a.length()
  let copy : Array[Array[Double]] = Array(capacity=rows)
  let mut i = 0
  while i < rows {
    copy.push(a[i].copy())
    i = i + 1
  }
  copy
}

///|
/// Whether `b` has one entry per row of the system, and the size it is.
fn checked_size(rows : Int, b : Array[Double]) -> Int raise LinError {
  if b.length() != rows {
    raise LinError::SizeMismatch(rows, b.length())
  }
  rows
}

///|
/// The LU factorization `P · A = L · U` of a square matrix, with partial
/// pivoting.
///
/// The two factors are held in one array: below the diagonal is `L` without
/// its unit diagonal, which is not stored because it is ones, and on and above
/// the diagonal is `U`. `pivots` holds the row permutation — row `i` of the
/// permuted matrix is row `pivots[i]` of `A` — and `swaps` counts the
/// exchanges, which is what the sign of the determinant turns on.
///
/// Nothing here is exposed: a decomposition is used through `solve`,
/// `determinant` and `inverse`, which are the questions it answers.
pub struct LU {
  priv factors : Array[Array[Double]]
  priv pivots : Array[Int]
  priv swaps : Int
}

///|
/// Factor `a` in place — in a copy of it, that is — by Gaussian elimination
/// with partial pivoting.
///
/// At each column the largest remaining element of that column is exchanged
/// into the diagonal before the column is eliminated with it. Without that,
/// a pivot that is small makes every multiplier below it large, and the
/// rounding of the elimination is multiplied by as much: `[[1e-18, 1], [1, 1]]`
/// has the exact solution `[1, 1]`, and elimination in the order the rows come
/// in answers `[0, 1]` — the pivot is not zero, so nothing complains.
pub fn LU::decompose(a : Array[Array[Double]]) -> LU raise LinError {
  let size = checked_square(a)
  let factors = copy_rows(a)
  let pivots : Array[Int] = Array(capacity=size)
  let mut i = 0
  while i < size {
    pivots.push(i)
    i = i + 1
  }
  let mut swaps = 0
  let mut column = 0
  while column < size {
    let mut pivot = column
    let mut row = column + 1
    while row < size {
      if factors[row][column].abs() > factors[pivot][column].abs() {
        pivot = row
      }
      row = row + 1
    }
    if factors[pivot][column] == 0.0 {
      raise LinError::Singular(column)
    }
    if pivot != column {
      let held = factors[pivot]
      factors[pivot] = factors[column]
      factors[column] = held
      let held_pivot = pivots[pivot]
      pivots[pivot] = pivots[column]
      pivots[column] = held_pivot
      swaps = swaps + 1
    }
    let diagonal = factors[column][column]
    row = column + 1
    while row < size {
      let multiplier = factors[row][column] / diagonal
      factors[row][column] = multiplier
      let mut j = column + 1
      while j < size {
        factors[row][j] = factors[row][j] - multiplier * factors[column][j]
        j = j + 1
      }
      row = row + 1
    }
    column = column + 1
  }
  { factors, pivots, swaps, }
}

///|
/// How wide the factored system is.
pub fn LU::size(self : LU) -> Int {
  self.factors.length()
}

///|
/// Solve `A · x = b` with the factors: permute `b` the way the rows were
/// permuted, substitute forward through `L`, and back through `U`.
pub fn LU::solve(self : LU, b : Array[Double]) -> Array[Double] raise LinError {
  let size = checked_size(self.size(), b)
  let y : Array[Double] = Array(capacity=size)
  let mut i = 0
  while i < size {
    let mut acc = b[self.pivots[i]]
    let mut j = 0
    while j < i {
      acc = acc - self.factors[i][j] * y[j]
      j = j + 1
    }
    y.push(acc)
    i = i + 1
  }
  let x : Array[Double] = Array(capacity=size)
  let mut k = 0
  while k < size {
    x.push(0.0)
    k = k + 1
  }
  let mut row = size - 1
  while row >= 0 {
    let mut acc = y[row]
    let mut j = row + 1
    while j < size {
      acc = acc - self.factors[row][j] * x[j]
      j = j + 1
    }
    x[row] = acc / self.factors[row][row]
    row = row - 1
  }
  x
}

///|
/// Solve `A · x = b` once for each right-hand side, in one pass over the
/// factors per column: the cost of factoring is paid once and the cost of the
/// substitutions is paid per right-hand side.
pub fn LU::solve_many(
  self : LU,
  bs : Array[Array[Double]],
) -> Array[Array[Double]] raise LinError {
  let solutions : Array[Array[Double]] = Array(capacity=bs.length())
  for b in bs {
    solutions.push(self.solve(b))
  }
  solutions
}

///|
/// The determinant, as the product of the diagonal of `U` with the sign the
/// row exchanges gave it. A singular matrix never gets this far: its
/// decomposition raises.
pub fn LU::determinant(self : LU) -> Double {
  let size = self.size()
  let mut det = 1.0
  let mut i = 0
  while i < size {
    det = det * self.factors[i][i]
    i = i + 1
  }
  if self.swaps % 2 == 1 {
    -det
  } else {
    det
  }
}

///|
/// The inverse, by solving `A · X = I` one unit vector at a time. It is the
/// same work as `n` solves, and is worth it only when the inverse is what is
/// wanted: solving for a right-hand side is cheaper and better conditioned
/// than multiplying by an inverse that had to be formed first.
pub fn LU::inverse(self : LU) -> Array[Array[Double]] raise LinError {
  let size = self.size()
  let identity : Array[Array[Double]] = Array(capacity=size)
  let mut i = 0
  while i < size {
    let unit : Array[Double] = Array(capacity=size)
    let mut j = 0
    while j < size {
      unit.push(if i == j { 1.0 } else { 0.0 })
      j = j + 1
    }
    identity.push(unit)
    i = i + 1
  }
  let columns = self.solve_many(identity)
  let inverse = zero_matrix(size, size)
  let mut row = 0
  while row < size {
    let mut column = 0
    while column < size {
      inverse[row][column] = columns[column][row]
      column = column + 1
    }
    row = row + 1
  }
  inverse
}

///|
/// Solve `A · x = b` by factoring `A` and substituting. The one-call form of
/// `LU::decompose` followed by `LU::solve`.
pub fn solve_linear_system(
  a : Array[Array[Double]],
  b : Array[Double],
) -> Array[Double] raise LinError {
  LU::decompose(a).solve(b)
}

///|
/// The Cholesky factorization `A = L · Lᵀ` of a symmetric positive definite
/// matrix, as the lower factor.
///
/// The diagonal entries are built by taking a square root, and the argument of
/// that square root is positive exactly when the matrix is positive definite:
/// a matrix that is symmetric and not positive definite — or not symmetric,
/// which this does not check — has a pivot that is zero or negative, and is
/// refused with the index of that pivot.
///
/// The work is `n³/6` against the LU's `n³/3`, and there is no pivoting to pay
/// for: the property that makes the pivots positive is what makes the method
/// stable without row exchanges.
pub struct Cholesky {
  priv lower : Array[Array[Double]]
}

///|
/// Factor a symmetric positive definite matrix.
///
/// Only the lower triangle of the argument is read: an upper triangle that
/// disagrees with it — a matrix that is not symmetric — is ignored rather
/// than averaged in, and the factors are those of the symmetric matrix the
/// lower triangle describes.
pub fn Cholesky::decompose(a : Array[Array[Double]]) -> Cholesky raise LinError {
  let size = checked_square(a)
  let lower = zero_matrix(size, size)
  let mut row = 0
  while row < size {
    let mut column = 0
    while column <= row {
      let mut acc = a[row][column]
      let mut k = 0
      while k < column {
        acc = acc - lower[row][k] * lower[column][k]
        k = k + 1
      }
      if row == column {
        if acc <= 0.0 {
          raise LinError::NotPositiveDefinite(row)
        }
        lower[row][column] = acc.sqrt()
      } else {
        lower[row][column] = acc / lower[column][column]
      }
      column = column + 1
    }
    row = row + 1
  }
  { lower, }
}

///|
/// How wide the factored system is.
pub fn Cholesky::size(self : Cholesky) -> Int {
  self.lower.length()
}

///|
/// Solve `A · x = b` with the factors: substitute forward through `L`, then
/// back through `Lᵀ`, which is the same array read the other way.
pub fn Cholesky::solve(
  self : Cholesky,
  b : Array[Double],
) -> Array[Double] raise LinError {
  let size = checked_size(self.size(), b)
  let y : Array[Double] = Array(capacity=size)
  let mut i = 0
  while i < size {
    let mut acc = b[i]
    let mut j = 0
    while j < i {
      acc = acc - self.lower[i][j] * y[j]
      j = j + 1
    }
    y.push(acc / self.lower[i][i])
    i = i + 1
  }
  let x : Array[Double] = Array(capacity=size)
  let mut k = 0
  while k < size {
    x.push(0.0)
    k = k + 1
  }
  let mut row = size - 1
  while row >= 0 {
    let mut acc = y[row]
    let mut j = row + 1
    while j < size {
      acc = acc - self.lower[j][row] * x[j]
      j = j + 1
    }
    x[row] = acc / self.lower[row][row]
    row = row - 1
  }
  x
}

///|
/// The determinant, as the square of the product of the diagonal of `L`: the
/// two triangular factors have the same diagonal, so `det(L · Lᵀ)` is
/// `det(L)²`, and it is positive the way a positive definite matrix's
/// determinant has to be.
pub fn Cholesky::determinant(self : Cholesky) -> Double {
  let size = self.size()
  let mut det = 1.0
  let mut i = 0
  while i < size {
    det = det * self.lower[i][i] * self.lower[i][i]
    i = i + 1
  }
  det
}

///|
/// Solve a tridiagonal system by the Thomas algorithm: eliminate the
/// sub-diagonal forward, then substitute back.
///
/// The matrix is given as its three diagonals rather than as a matrix — the
/// sub-diagonal first, then the diagonal, then the super-diagonal, each with
/// one fewer entry than the diagonal — because a tridiagonal matrix stored as
/// rows is mostly zeros, and `O(n)` work on three arrays is what the shape is
/// worth.
///
///     | d₀ u₀  0 |      lower = [l₀, l₁]
///     | l₀ d₁ u₁ |      diagonal = [d₀, d₁, d₂]
///     |  0 l₁ d₂ |      upper = [u₀, u₁]
///
/// A zero on the diagonal after elimination is a zero pivot and nothing to
/// exchange it with within the band: the matrix is singular, and the index of
/// that entry is what the error carries.
pub fn solve_tridiagonal(
  lower : Array[Double],
  diagonal : Array[Double],
  upper : Array[Double],
  rhs : Array[Double],
) -> Array[Double] raise LinError {
  let size = checked_size(diagonal.length(), rhs)
  let bands = if size == 0 { 0 } else { size - 1 }
  if lower.length() != bands {
    raise LinError::SizeMismatch(bands, lower.length())
  }
  if upper.length() != bands {
    raise LinError::SizeMismatch(bands, upper.length())
  }
  let modified_upper : Array[Double] = Array(capacity=size)
  let modified_rhs : Array[Double] = Array(capacity=size)
  let mut i = 0
  while i < size {
    let mut pivot = diagonal[i]
    if i > 0 {
      pivot = pivot - lower[i - 1] * modified_upper[i - 1]
    }
    if pivot == 0.0 {
      raise LinError::Singular(i)
    }
    modified_upper.push(if i < bands { upper[i] / pivot } else { 0.0 })
    modified_rhs.push(
      if i == 0 {
        rhs[i] / pivot
      } else {
        (rhs[i] - lower[i - 1] * modified_rhs[i - 1]) / pivot
      },
    )
    i = i + 1
  }
  let x : Array[Double] = Array(capacity=size)
  let mut k = 0
  while k < size {
    x.push(0.0)
    k = k + 1
  }
  let mut row = size - 1
  while row >= 0 {
    x[row] = modified_rhs[row] -
      modified_upper[row] * (if row < bands { x[row + 1] } else { 0.0 })
    row = row - 1
  }
  x
}
