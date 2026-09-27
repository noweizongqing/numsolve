///|
/// Interpolation tests: a table of knots and the values a line and a spline
/// take between them.
///
/// Every golden below was computed first with SciPy 1.18's `CubicSpline` —
/// natural, clamped as `bc_type=((1, 0.0), (1, 0.0))`, and not-a-knot — and
/// then written down in full, so the numbers in the assertions are the
/// reference implementation's rather than this one's.

///|
/// The tolerance the spline goldens are held to. SciPy's elimination and this
/// one's Thomas recurrence differ in the order they add the same terms, so the
/// agreement is to the last few digits rather than exact.
const INTERP_EPS : Double = 1.0e-12

///|
/// Whether two numbers agree to `tolerance`.
fn interp_close(actual : Double, expected : Double, tolerance : Double) -> Bool {
  (actual - expected).abs() <= tolerance
}

///|
/// The knots the spline goldens are sampled on.
fn interp_xs() -> Array[Double] {
  [0.0, 1.0, 2.0, 3.0]
}

///|
/// The values at those knots.
fn interp_ys() -> Array[Double] {
  [1.0, 2.0, 0.0, 4.0]
}

///|
/// The same values at knots that are not evenly spaced.
fn interp_uneven_xs() -> Array[Double] {
  [0.0, 1.0, 3.0, 4.0]
}

///|
/// The values the exponential-spline test samples.
fn interp_exp_ys() -> Array[Double] {
  [@math.exp(0.0), @math.exp(0.25), @math.exp(0.5), @math.exp(1.0)]
}

///|
/// A spline the tests know is well formed.
fn interp_spline(
  xs : Array[Double],
  ys : Array[Double],
  boundary : SplineBoundary,
) -> CubicSpline {
  CubicSpline::new(xs, ys, boundary) catch {
    _ => abort("the golden table is a valid one")
  }
}

///|
/// Whether building a table raises `TooFewKnots` with `count` points.
fn interp_too_few(f : () -> Unit raise InterpError, count : Int) -> Bool {
  try {
    f()
    false
  } catch {
    InterpError::TooFewKnots(got) => got == count
    _ => false
  }
}

///|
/// Whether building a table raises `UnevenLengths` with the two lengths.
fn interp_uneven(
  f : () -> Unit raise InterpError,
  knots : Int,
  values : Int,
) -> Bool {
  try {
    f()
    false
  } catch {
    InterpError::UnevenLengths(got_knots, got_values) =>
      got_knots == knots && got_values == values
    _ => false
  }
}

///|
/// Whether building a table raises `NotSorted` at `knot`.
fn interp_unsorted(f : () -> Unit raise InterpError, knot : Double) -> Bool {
  try {
    f()
    false
  } catch {
    InterpError::NotSorted(got) => got == knot
    _ => false
  }
}

///|
/// Whether building a table raises `DuplicateKnots` at `knot`.
fn interp_duplicate(f : () -> Unit raise InterpError, knot : Double) -> Bool {
  try {
    f()
    false
  } catch {
    InterpError::DuplicateKnots(got) => got == knot
    _ => false
  }
}

///|
/// The value of one row of `coefficients` at `x`, with the row read as a
/// polynomial in the global `x`.
fn interp_polynomial(row : Array[Double], x : Double) -> Double {
  row[0] * x * x * x + row[1] * x * x + row[2] * x + row[3]
}

///|
/// The slope of that polynomial.
fn interp_slope(row : Array[Double], x : Double) -> Double {
  3.0 * row[0] * x * x + 2.0 * row[1] * x + row[2]
}

///|
/// The second derivative of that polynomial.
fn interp_curvature(row : Array[Double], x : Double) -> Double {
  6.0 * row[0] * x + 2.0 * row[1]
}

///|
/// The value a line takes between the knots of the golden table — the average
/// of the two ends, because the knots are one apart.
test "a line through the golden table" {
  let line = LinearInterp::new(interp_xs(), interp_ys()) catch {
    _ => abort("the golden table is a valid one")
  }
  assert_true(interp_close(line.evaluate(0.5).unwrap(), 1.5, INTERP_EPS))
  assert_true(interp_close(line.evaluate(1.5).unwrap(), 1.0, INTERP_EPS))
  assert_true(interp_close(line.evaluate(2.5).unwrap(), 2.0, INTERP_EPS))
  assert_true(interp_close(line.evaluate(2.0).unwrap(), 0.0, INTERP_EPS))
  assert_true(interp_close(line.evaluate(0.25).unwrap(), 1.25, INTERP_EPS))
}

///|
/// The natural spline through the golden table: SciPy's samples at the three
/// midpoints, and the moments it solves for.
test "a natural cubic spline" {
  let spline = interp_spline(interp_xs(), interp_ys(), SplineBoundary::Natural)
  assert_true(
    interp_close(spline.evaluate(0.5).unwrap(), 1.9500000000000002, INTERP_EPS),
  )
  assert_true(
    interp_close(spline.evaluate(1.5).unwrap(), 0.77499999999999991, INTERP_EPS),
  )
  assert_true(
    interp_close(spline.evaluate(2.5).unwrap(), 1.3250000000000002, INTERP_EPS),
  )
  let moments = spline.moments()
  assert_true(interp_close(moments[0], 0.0, INTERP_EPS))
  assert_true(interp_close(moments[1], -7.2, INTERP_EPS))
  assert_true(interp_close(moments[2], 10.8, INTERP_EPS))
  assert_true(interp_close(moments[3], 0.0, INTERP_EPS))
}

///|
/// The clamped spline through the same table with both end slopes zero: the
/// values SciPy returns for it, and the slopes the fit has to have at the two
/// ends.
test "a clamped cubic spline" {
  let spline = interp_spline(
    interp_xs(),
    interp_ys(),
    SplineBoundary::Clamped(0.0, 0.0),
  )
  assert_true(interp_close(spline.evaluate(0.5).unwrap(), 1.65, INTERP_EPS))
  assert_true(interp_close(spline.evaluate(1.5).unwrap(), 0.625, INTERP_EPS))
  assert_true(interp_close(spline.evaluate(2.5).unwrap(), 2.225, INTERP_EPS))
  let rows = spline.coefficients()
  assert_true(interp_close(interp_slope(rows[0], 0.0), 0.0, INTERP_EPS))
  assert_true(interp_close(interp_slope(rows[2], 3.0), 0.0, INTERP_EPS))
}

///|
/// The not-a-knot spline: with no condition at the ends at all, the two outer
/// cubics are one polynomial.
test "a not-a-knot cubic spline" {
  let spline = interp_spline(interp_xs(), interp_ys(), SplineBoundary::NotAKnot)
  assert_true(interp_close(spline.evaluate(0.5).unwrap(), 2.4375, INTERP_EPS))
  assert_true(
    interp_close(spline.evaluate(1.5).unwrap(), 0.81249999999999989, INTERP_EPS),
  )
  assert_true(
    interp_close(spline.evaluate(2.5).unwrap(), 0.68750000000000022, INTERP_EPS),
  )
  let rows = spline.coefficients()
  assert_true(
    interp_close(
      interp_curvature(rows[0], 1.0),
      interp_curvature(rows[1], 1.0),
      INTERP_EPS,
    ),
  )
}

///|
/// The same table with the middle knot moved: the bands of the tridiagonal
/// system are no longer equal, and the answer is SciPy's for that table.
test "a spline through unevenly spaced knots" {
  let natural = interp_spline(
    interp_uneven_xs(),
    interp_ys(),
    SplineBoundary::Natural,
  )
  assert_true(
    interp_close(natural.evaluate(0.5).unwrap(), 1.7578125, INTERP_EPS),
  )
  assert_true(interp_close(natural.evaluate(2.0).unwrap(), 0.4375, INTERP_EPS))
  assert_true(
    interp_close(natural.evaluate(3.5).unwrap(), 1.6015625, INTERP_EPS),
  )
  let free = interp_spline(
    interp_uneven_xs(),
    interp_ys(),
    SplineBoundary::NotAKnot,
  )
  assert_true(interp_close(free.evaluate(0.5).unwrap(), 2.03125, INTERP_EPS))
  assert_true(
    interp_close(free.evaluate(2.0).unwrap(), 0.50000000000000022, INTERP_EPS),
  )
  assert_true(
    interp_close(free.evaluate(3.5).unwrap(), 1.2187499999999998, INTERP_EPS),
  )
}

///|
/// What makes a spline a spline: at every interior knot the two cubics that
/// meet there agree in value, in slope and in second derivative.
test "a spline is twice differentiable across every knot" {
  let boundaries : Array[SplineBoundary] = [
    SplineBoundary::Natural,
    SplineBoundary::Clamped(0.0, 0.0),
    SplineBoundary::NotAKnot,
  ]
  let xs = interp_xs()
  for boundary in boundaries {
    let rows = interp_spline(xs, interp_ys(), boundary).coefficients()
    let mut i = 0
    while i < rows.length() - 1 {
      let knot = xs[i + 1]
      assert_true(
        interp_close(
          interp_polynomial(rows[i], knot),
          interp_polynomial(rows[i + 1], knot),
          INTERP_EPS,
        ),
      )
      assert_true(
        interp_close(
          interp_slope(rows[i], knot),
          interp_slope(rows[i + 1], knot),
          INTERP_EPS,
        ),
      )
      assert_true(
        interp_close(
          interp_curvature(rows[i], knot),
          interp_curvature(rows[i + 1], knot),
          INTERP_EPS,
        ),
      )
      i = i + 1
    }
  }
}

///|
/// The moments are the second derivatives the cubics are built from, at the
/// knots they are indexed by.
test "the moments are the second derivatives at the knots" {
  let xs = interp_xs()
  let spline = interp_spline(xs, interp_ys(), SplineBoundary::Natural)
  let rows = spline.coefficients()
  let moments = spline.moments()
  let mut i = 1
  while i < xs.length() - 1 {
    assert_true(
      interp_close(interp_curvature(rows[i - 1], xs[i]), moments[i], INTERP_EPS),
    )
    i = i + 1
  }
}

///|
/// A spline passes through the points it was built from, exactly, at every
/// knot — the second derivatives it solves for are chosen to make that so.
test "a spline passes through its knots" {
  let xs = interp_xs()
  let ys = interp_ys()
  let spline = interp_spline(xs, ys, SplineBoundary::Natural)
  let mut i = 0
  while i < xs.length() {
    assert_eq(spline.evaluate(xs[i]).unwrap(), ys[i])
    i = i + 1
  }
}

///|
/// A straight line is the one table every boundary condition agrees on: the
/// cubics collapse to the line, whatever the ends are told to do.
test "a table on a straight line is interpolated by that line" {
  let xs = interp_xs()
  let ys : Array[Double] = [1.0, 3.0, 5.0, 7.0]
  let boundaries : Array[SplineBoundary] = [
    SplineBoundary::Natural,
    SplineBoundary::Clamped(2.0, 2.0),
    SplineBoundary::NotAKnot,
  ]
  for boundary in boundaries {
    let spline = interp_spline(xs, ys, boundary)
    assert_true(interp_close(spline.evaluate(1.5).unwrap(), 4.0, INTERP_EPS))
    assert_true(interp_close(spline.evaluate(2.25).unwrap(), 5.5, INTERP_EPS))
    for row in spline.coefficients() {
      assert_true(interp_close(row[0], 0.0, INTERP_EPS))
      assert_true(interp_close(row[1], 0.0, INTERP_EPS))
    }
  }
}

///|
/// A table of two points is one interval: a line under either of the two
/// boundaries that need no third knot, and no not-a-knot spline at all.
test "a table of two points" {
  let xs : Array[Double] = [0.0, 1.0]
  let ys : Array[Double] = [1.0, 3.0]
  let natural = interp_spline(xs, ys, SplineBoundary::Natural)
  assert_true(interp_close(natural.evaluate(0.25).unwrap(), 1.5, INTERP_EPS))
  let clamped = interp_spline(xs, ys, SplineBoundary::Clamped(2.0, 2.0))
  assert_true(interp_close(clamped.evaluate(0.75).unwrap(), 2.5, INTERP_EPS))
  assert_true(
    interp_too_few(
      () => ignore(CubicSpline::new(xs, ys, SplineBoundary::NotAKnot)),
      2,
    ),
  )
}

///|
/// Outside the range of the table there is no curve to read: a line describes
/// the interval it was built from, not the line beyond it.
test "a point outside the table has no value" {
  let line = LinearInterp::new(interp_xs(), interp_ys()) catch {
    _ => abort("the golden table is a valid one")
  }
  let spline = interp_spline(interp_xs(), interp_ys(), SplineBoundary::Natural)
  assert_true(line.evaluate(-0.000001) is None)
  assert_true(line.evaluate(3.000001) is None)
  assert_true(spline.evaluate(-0.000001) is None)
  assert_true(spline.evaluate(3.000001) is None)
  assert_true(line.evaluate(0.0) is Some(_))
  assert_true(spline.evaluate(3.0) is Some(_))
}

///|
/// A table that is not a table: too short, two arrays of different lengths,
/// knots out of order, and a knot repeated.
test "a table that cannot be interpolated" {
  assert_true(interp_too_few(() => ignore(LinearInterp::new([0.0], [1.0])), 1))
  assert_true(
    interp_uneven(
      () => ignore(LinearInterp::new([0.0, 1.0, 2.0], [1.0, 2.0])),
      3,
      2,
    ),
  )
  assert_true(
    interp_unsorted(
      () => ignore(LinearInterp::new([0.0, 2.0, 1.0], [1.0, 2.0, 3.0])),
      1.0,
    ),
  )
  assert_true(
    interp_duplicate(
      () => ignore(LinearInterp::new([0.0, 1.0, 1.0], [1.0, 2.0, 3.0])),
      1.0,
    ),
  )
  assert_true(
    interp_too_few(
      () => ignore(CubicSpline::new([0.0], [1.0], SplineBoundary::Natural)),
      1,
    ),
  )
  assert_true(
    interp_uneven(
      () => {
        ignore(
          CubicSpline::new([0.0, 1.0, 2.0], [1.0, 2.0], SplineBoundary::Natural),
        )
      },
      3,
      2,
    ),
  )
  assert_true(
    interp_unsorted(
      () => {
        ignore(
          CubicSpline::new(
            [0.0, 2.0, 1.0],
            [1.0, 2.0, 3.0],
            SplineBoundary::Natural,
          ),
        )
      },
      1.0,
    ),
  )
  assert_true(
    interp_duplicate(
      () => {
        ignore(
          CubicSpline::new(
            [0.0, 1.0, 1.0],
            [1.0, 2.0, 3.0],
            SplineBoundary::Natural,
          ),
        )
      },
      1.0,
    ),
  )
}

///|
/// A table of the exponential function, with the two end slopes it really has:
/// the spline through four samples of it, and how far from the curve it is
/// between them — a quarter of a percent at worst on the widest interval.
test "a spline through samples of the exponential" {
  let xs : Array[Double] = [0.0, 0.25, 0.5, 1.0]
  let ys = interp_exp_ys()
  let spline = interp_spline(
    xs,
    ys,
    SplineBoundary::Clamped(1.0, @math.exp(1.0)),
  )
  assert_true(
    interp_close(spline.evaluate(0.125).unwrap(), 1.1331261323836186, 1.0e-10),
  )
  assert_true(
    interp_close(spline.evaluate(0.75).unwrap(), 2.1165572806439754, 1.0e-10),
  )
  assert_true(
    (spline.evaluate(0.125).unwrap() - @math.exp(0.125)).abs() < 1.0e-4,
  )
  assert_true((spline.evaluate(0.75).unwrap() - @math.exp(0.75)).abs() < 1.0e-3)
  assert_true(spline.evaluate(0.125).unwrap() > ys[0])
  assert_true(spline.evaluate(0.75).unwrap() < ys[3])
}
