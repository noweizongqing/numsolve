///|
/// Interpolation: a curve through a table of points, and the value it takes
/// between them.
///
/// A table is two arrays, `xs` and `ys`, of the same length, with the `xs` in
/// increasing order: `LinearInterp` draws a straight line between each pair of
/// points, and `CubicSpline` fits one cubic per interval so that the curve, its
/// slope and its second derivative are all continuous across every knot.
///
/// `LinearInterp` is the line through the two knots that surround the point
/// asked for, and nothing else — the cheapest interpolation there is, with a
/// kink at every knot.
///
/// `CubicSpline` is the piecewise cubic whose second derivatives `M[i]` at the
/// knots are given by `h[i-1] · M[i-1] + 2(h[i-1] + h[i]) · M[i] + h[i] · M[i+1]
/// = 6 · ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1])`, a tridiagonal
/// system with `n` unknowns and `n - 2` equations — which is why a spline needs
/// two more conditions, and why `SplineBoundary` exists:
///
/// | `SplineBoundary` | The two extra conditions | The end conditions it gives |
/// | --- | --- | --- |
/// | `Natural` | `M[0] = M[n-1] = 0` | a curve that is free to bend: the one to use when nothing is known about the ends |
/// | `Clamped(d0, dn)` | `S'(x0) = d0`, `S'(xn) = dn` | a curve with prescribed slopes at the two ends |
/// | `NotAKnot` | `S'''` is continuous at the second and the second-to-last knot | no assumption about the ends at all: the two outermost cubics are the same polynomial |
///
/// The system is solved by `solve_tridiagonal` from `linear.mbt`, so a spline
/// through `n` points costs `O(n)` rather than the `O(n³)` of a dense
/// elimination.
///
///     (0, 1) (1, 2) (2, 0) (3, 4), natural
///     S(0.5) = 1.95    S(1.5) = 0.775    S(2.5) = 1.325
///
/// A point outside the range of the table has no answer — a spline and a line
/// both describe the interval they were built from, not the line beyond it —
/// so `evaluate` returns `None` there instead of extrapolating.

///|
/// A table that cannot be interpolated.
pub suberror InterpError {
  /// Fewer than two points were given: there is no interval to interpolate
  /// over. The payload is the number of points.
  TooFewKnots(Int)
  /// `xs` and `ys` are not the same length. The payloads are the two lengths.
  UnevenLengths(Int, Int)
  /// The knots are not in increasing order: a line or a spline through a table
  /// that doubles back is not a function of `x`. The payload is the first knot
  /// that is out of order.
  NotSorted(Double)
  /// Two knots are the same point, which would ask for two values at it. The
  /// payload is the knot that repeats.
  DuplicateKnots(Double)
}

///|
/// Check the two arrays a table is given as, and return the number of points.
fn checked_table(
  xs : Array[Double],
  ys : Array[Double],
) -> Int raise InterpError {
  let size = xs.length()
  if size < 2 {
    raise InterpError::TooFewKnots(size)
  }
  if ys.length() != size {
    raise InterpError::UnevenLengths(size, ys.length())
  }
  let mut i = 1
  while i < size {
    if xs[i] < xs[i - 1] {
      raise InterpError::NotSorted(xs[i])
    }
    if xs[i] == xs[i - 1] {
      raise InterpError::DuplicateKnots(xs[i])
    }
    i = i + 1
  }
  size
}

///|
/// The index of the interval a point falls in: the largest `i` with
/// `xs[i] <= x`, capped at the last interval so that the right endpoint
/// belongs to the interval ending there.
fn locate(xs : Array[Double], x : Double) -> Int {
  let mut low = 0
  let mut high = xs.length() - 1
  while low < high {
    let middle = (low + high + 1) / 2
    if xs[middle] <= x {
      low = middle
    } else {
      high = middle - 1
    }
  }
  if low > xs.length() - 2 {
    xs.length() - 2
  } else {
    low
  }
}

///|
/// Whether `x` is within the range the table was built from.
fn inside_range(xs : Array[Double], x : Double) -> Bool {
  x >= xs[0] && x <= xs[xs.length() - 1]
}

///|
/// The straight line through a table of points.
pub struct LinearInterp {
  priv xs : Array[Double]
  priv ys : Array[Double]
}

///|
/// Build the line through `xs` and `ys`.
pub fn LinearInterp::new(
  xs : Array[Double],
  ys : Array[Double],
) -> LinearInterp raise InterpError {
  ignore(checked_table(xs, ys))
  { xs: xs.copy(), ys: ys.copy(), }
}

///|
/// The value of the line at `x`, or `None` if `x` is outside the table.
pub fn LinearInterp::evaluate(self : LinearInterp, x : Double) -> Double? {
  if !inside_range(self.xs, x) {
    return None
  }
  let i = locate(self.xs, x)
  let width = self.xs[i + 1] - self.xs[i]
  let weight = (x - self.xs[i]) / width
  Some(self.ys[i] * (1.0 - weight) + self.ys[i + 1] * weight)
}

///|
/// How the two ends of a cubic spline are pinned down: the two conditions that
/// turn `n - 2` equations in `n` unknowns into a system.
pub(all) enum SplineBoundary {
  /// The second derivative is zero at both ends — the curve a ruler would take
  /// through the points if it were free at the ends.
  Natural
  /// The first derivative is `d0` at the first knot and `dn` at the last.
  Clamped(Double, Double)
  /// No end condition: the third derivative is continuous at the second and
  /// the second-to-last knot, so the two outermost cubics are one polynomial
  /// and those two knots are not knots at all.
  NotAKnot
}

///|
/// A cubic spline: one cubic per interval, joined so that the value, the slope
/// and the second derivative agree at every knot.
pub struct CubicSpline {
  priv xs : Array[Double]
  priv ys : Array[Double]
  /// The second derivatives at the knots, the unknowns of the tridiagonal
  /// system the boundary conditions complete.
  priv moments : Array[Double]
}

///|
/// The widths of the intervals of a table: `h[i] = xs[i + 1] - xs[i]`.
fn interval_widths(xs : Array[Double]) -> Array[Double] {
  let widths : Array[Double] = Array(capacity=xs.length() - 1)
  let mut i = 0
  while i < xs.length() - 1 {
    widths.push(xs[i + 1] - xs[i])
    i = i + 1
  }
  widths
}

///|
/// Six times the change in slope across a knot — `6 · (s[i] - s[i-1])` for the
/// slopes of the two intervals meeting there — which is the right-hand side of
/// the equation that knot contributes.
fn knot_rhs(ys : Array[Double], widths : Array[Double], knot : Int) -> Double {
  6.0 *
  (
    (ys[knot + 1] - ys[knot]) / widths[knot] -
    (ys[knot] - ys[knot - 1]) / widths[knot - 1]
  )
}

///|
/// The tridiagonal system over the second derivatives under a boundary
/// condition that names one end curvature each: row 0 for the first knot and
/// row `n - 1` for the last, beside the interior rows for knots 1 to `n - 2`.
fn spline_system(
  xs : Array[Double],
  ys : Array[Double],
  boundary : SplineBoundary,
) -> (Array[Double], Array[Double], Array[Double], Array[Double]) {
  let size = xs.length()
  let widths = interval_widths(xs)
  let lower : Array[Double] = Array(capacity=size - 1)
  let diagonal : Array[Double] = Array(capacity=size)
  let upper : Array[Double] = Array(capacity=size - 1)
  let rhs : Array[Double] = Array(capacity=size)
  let mut i = 0
  while i < size {
    diagonal.push(0.0)
    rhs.push(0.0)
    if i < size - 1 {
      lower.push(widths[i])
      upper.push(widths[i])
    }
    i = i + 1
  }
  let mut knot = 1
  while knot < size - 1 {
    diagonal[knot] = 2.0 * (widths[knot - 1] + widths[knot])
    rhs[knot] = knot_rhs(ys, widths, knot)
    knot = knot + 1
  }
  match boundary {
    Natural => {
      diagonal[0] = 1.0
      upper[0] = 0.0
      rhs[0] = 0.0
      diagonal[size - 1] = 1.0
      lower[size - 2] = 0.0
      rhs[size - 1] = 0.0
    }
    Clamped(slope0, slopen) => {
      let first = widths[0]
      diagonal[0] = 2.0 * first
      upper[0] = first
      rhs[0] = 6.0 * ((ys[1] - ys[0]) / first - slope0)
      let last = widths[size - 2]
      lower[size - 2] = last
      diagonal[size - 1] = 2.0 * last
      rhs[size - 1] = 6.0 * (slopen - (ys[size - 1] - ys[size - 2]) / last)
    }
    NotAKnot => abort("the not-a-knot condition is not a row of this system")
  }
  let rows = (lower, diagonal, upper, rhs)
  rows
}

///|
/// The second derivatives of a not-a-knot spline.
///
/// The condition is not a row of the system the other two boundaries write
/// into: "the third derivative is continuous at the knot next to the end"
/// relates three unknowns — `M[0]`, `M[1]` and `M[2]` — where row 0 has room
/// for two. So `M[0]` is eliminated instead of being solved for, using
/// `(M[1] - M[0]) / h[0] = (M[2] - M[1]) / h[1]`, and the row of the first
/// interior knot takes its place: the unknowns are then `M[1]` to `M[n - 2]`,
/// and the two end curvatures are read back off the relation once they are
/// known. The same is done at the other end, where the row of knot `n - 2`
/// absorbs `M[n - 1]`.
///
/// A table of three points has one knot to eliminate and no third derivative
/// to make continuous that is not already: the cubic through three points with
/// no curvature at either end of it is a parabola, and its second derivative is
/// the same at all three — which is what SciPy answers with as well.
fn not_a_knot_moments(xs : Array[Double], ys : Array[Double]) -> Array[Double] {
  let size = xs.length()
  let widths = interval_widths(xs)
  if size == 3 {
    let flat = knot_rhs(ys, widths, 1) / (3.0 * (widths[0] + widths[1]))
    return [flat, flat, flat]
  }
  let reduced = size - 2
  let lower : Array[Double] = Array(capacity=reduced - 1)
  let diagonal : Array[Double] = Array(capacity=reduced)
  let upper : Array[Double] = Array(capacity=reduced - 1)
  let rhs : Array[Double] = Array(capacity=reduced)
  let mut index = 0
  while index < reduced {
    let knot = index + 1
    if knot == 1 {
      let first = widths[0]
      let next = widths[1]
      diagonal.push((first + next) * (first + 2.0 * next))
      upper.push(next * next - first * first)
      rhs.push(knot_rhs(ys, widths, knot) * next)
    } else if knot == size - 2 {
      let last = widths[size - 2]
      let before = widths[size - 3]
      lower.push(before * before - last * last)
      diagonal.push((before + last) * (2.0 * before + last))
      rhs.push(knot_rhs(ys, widths, knot) * before)
    } else {
      lower.push(widths[knot - 1])
      diagonal.push(2.0 * (widths[knot - 1] + widths[knot]))
      upper.push(widths[knot])
      rhs.push(knot_rhs(ys, widths, knot))
    }
    index = index + 1
  }
  let middle = solve_tridiagonal(lower, diagonal, upper, rhs) catch {
    _ => abort("a spline system is diagonally dominant")
  }
  let moments : Array[Double] = Array(capacity=size)
  moments.push(middle[0] - widths[0] / widths[1] * (middle[1] - middle[0]))
  for value in middle {
    moments.push(value)
  }
  let last = reduced - 1
  moments.push(
    middle[last] +
    widths[size - 2] / widths[size - 3] * (middle[last] - middle[last - 1]),
  )
  moments
}

///|
/// Fit a cubic spline through `xs` and `ys` under the boundary condition.
///
/// A `NotAKnot` spline needs three knots to have a knot that is not one: with
/// two points it is refused, and with the other two boundaries a table of two
/// points is a straight line, whose second derivatives are zero.
pub fn CubicSpline::new(
  xs : Array[Double],
  ys : Array[Double],
  boundary : SplineBoundary,
) -> CubicSpline raise InterpError {
  let size = checked_table(xs, ys)
  let moments = match boundary {
    NotAKnot =>
      if size < 3 {
        raise InterpError::TooFewKnots(size)
      } else {
        not_a_knot_moments(xs, ys)
      }
    _ =>
      if size == 2 {
        [0.0, 0.0]
      } else {
        let system = spline_system(xs, ys, boundary)
        // The system is diagonally dominant with strictly increasing knots,
        // and irreducible — every row has an off-diagonal entry — so it is
        // never singular and the elimination never divides by zero.
        solve_tridiagonal(system.0, system.1, system.2, system.3) catch {
          _ => abort("a spline system is diagonally dominant")
        }
      }
  }
  { xs: xs.copy(), ys: ys.copy(), moments, }
}

///|
/// The second derivatives at the knots, in the order the knots are given: the
/// `M` of the system above, and the same numbers the cubics are built from.
///
/// They are zero at both ends exactly when the spline is natural.
pub fn CubicSpline::moments(self : CubicSpline) -> Array[Double] {
  self.moments.copy()
}

///|
/// The coefficients of the cubic on each interval as a polynomial in the
/// global `x`, one row `[a, b, c, d]` per interval of `a · x³ + b · x² + c · x + d`.
///
/// Two neighbouring rows describe the same value, the same slope and the same
/// second derivative at the knot they share — which is what makes the spline
/// `C²`, and what a reader checking the fit rather than sampling it wants.
pub fn CubicSpline::coefficients(self : CubicSpline) -> Array[Array[Double]] {
  let rows : Array[Array[Double]] = Array(capacity=self.xs.length() - 1)
  let mut i = 0
  while i < self.xs.length() - 1 {
    let width = self.xs[i + 1] - self.xs[i]
    let slope = (self.ys[i + 1] - self.ys[i]) / width
    let left = self.moments[i]
    let right = self.moments[i + 1]
    let a = (right - left) / (6.0 * width)
    let b = left / 2.0
    let c = slope - width * (2.0 * left + right) / 6.0
    let d = self.ys[i]
    rows.push([
      a,
      b - 3.0 * a * self.xs[i],
      c - 2.0 * b * self.xs[i] + 3.0 * a * self.xs[i] * self.xs[i],
      d -
      c * self.xs[i] +
      b * self.xs[i] * self.xs[i] -
      a * self.xs[i] * self.xs[i] * self.xs[i],
    ])
    i = i + 1
  }
  rows
}

///|
/// The value of the spline at `x`, or `None` if `x` is outside the table.
pub fn CubicSpline::evaluate(self : CubicSpline, x : Double) -> Double? {
  if !inside_range(self.xs, x) {
    return None
  }
  let i = locate(self.xs, x)
  let width = self.xs[i + 1] - self.xs[i]
  let offset = x - self.xs[i]
  Some(
    self.ys[i] +
    offset *
    (
      (self.ys[i + 1] - self.ys[i]) / width -
      width * (2.0 * self.moments[i] + self.moments[i + 1]) / 6.0 +
      offset *
      (
        self.moments[i] / 2.0 +
        offset * (self.moments[i + 1] - self.moments[i]) / (6.0 * width)
      )
    ),
  )
}
