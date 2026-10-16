///|
/// Banded linear systems: a dense matrix whose non-zero entries are confined
/// to a band around the diagonal, solved without walking the full dense
/// elimination.
///
/// The API takes the matrix in dense row form and validates the declared
/// bandwidth. The elimination then searches for a pivot only within the lower
/// bandwidth and updates only the columns that the band and its fill-in can
/// reach, so a system with bandwidths `p` and `q` costs `O(n · (p+q)²)`
/// instead of the `O(n³)` of the dense factorization. This is the same
/// structural idea as `solve_tridiagonal`, with a general pair of
/// bandwidths.

///|
/// A vector of `size` zeros.
fn banded_zeros(size : Int32) -> Array[Double] {
  let values : Array[Double] = Array(capacity=size)
  let mut i = 0
  while i < size {
    values.push(0.0)
    i = i + 1
  }
  values
}

///|
/// Whether `a` has no non-zero entry outside the band described by the two
/// bandwidths.
fn checked_banded(
  a : Array[Array[Double]],
  lower_bandwidth : Int32,
  upper_bandwidth : Int32,
) -> Int32 raise LinError {
  let size = checked_square(a)
  let mut row = 0
  while row < size {
    let lower_edge = if row > lower_bandwidth {
      row - lower_bandwidth
    } else {
      0
    }
    let upper_edge = if row + upper_bandwidth < size - 1 {
      row + upper_bandwidth
    } else {
      size - 1
    }
    let mut column = 0
    while column < size {
      if (column < lower_edge || column > upper_edge) && a[row][column] != 0.0 {
        raise LinError::NotBanded(row, column)
      }
      column = column + 1
    }
    row = row + 1
  }
  size
}

///|
/// Solve a banded system `A · x = b` by elimination with partial pivoting
/// restricted to the lower bandwidth.
///
/// `lower_bandwidth` is the number of diagonals below the main diagonal that
/// may be non-zero, and `upper_bandwidth` is the number above it. Entries
/// outside that band are refused with `NotBanded`, because a dense matrix
/// that violates the declared structure would otherwise be solved with the
/// wrong arithmetic assumptions.
pub fn solve_banded(
  lower_bandwidth : Int32,
  upper_bandwidth : Int32,
  a : Array[Array[Double]],
  b : Array[Double],
) -> Array[Double] raise LinError {
  if lower_bandwidth < 0 {
    raise LinError::InvalidBandwidth(lower_bandwidth)
  }
  if upper_bandwidth < 0 {
    raise LinError::InvalidBandwidth(upper_bandwidth)
  }
  let size = checked_banded(a, lower_bandwidth, upper_bandwidth)
  let _ = checked_size(size, b)
  if size == 0 {
    return []
  }
  let factors = copy_rows(a)
  let rhs = b.copy()
  let mut column = 0
  while column < size {
    let last_row = if column + lower_bandwidth < size - 1 {
      column + lower_bandwidth
    } else {
      size - 1
    }
    let mut pivot = column
    let mut row = column + 1
    while row <= last_row {
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
      let held_rhs = rhs[pivot]
      rhs[pivot] = rhs[column]
      rhs[column] = held_rhs
    }
    let pivot_value = factors[column][column]
    row = column + 1
    while row <= last_row {
      let multiplier = factors[row][column] / pivot_value
      factors[row][column] = 0.0
      let last_column = if column + upper_bandwidth + lower_bandwidth < size - 1 {
        column + upper_bandwidth + lower_bandwidth
      } else {
        size - 1
      }
      let mut j = column + 1
      while j <= last_column {
        factors[row][j] = factors[row][j] - multiplier * factors[column][j]
        j = j + 1
      }
      rhs[row] = rhs[row] - multiplier * rhs[column]
      row = row + 1
    }
    column = column + 1
  }
  let x = banded_zeros(size)
  let mut row = size - 1
  while row >= 0 {
    let mut acc = rhs[row]
    let last_column = if row + upper_bandwidth + lower_bandwidth < size - 1 {
      row + upper_bandwidth + lower_bandwidth
    } else {
      size - 1
    }
    let mut j = row + 1
    while j <= last_column {
      acc = acc - factors[row][j] * x[j]
      j = j + 1
    }
    x[row] = acc / factors[row][row]
    row = row - 1
  }
  x
}
