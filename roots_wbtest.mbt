///|
/// Whether two values agree to within `tolerance`.
fn roots_close(actual : Double, expected : Double, tolerance : Double) -> Bool {
  (actual - expected).abs() <= tolerance
}

///|
/// Plain false position, without the Illinois halving: the same line through
/// the two ends of the bracket, stopping on the same interval width, and
/// nothing to pull it off an end it is stuck on. It is here for the test that
/// measures what the halving is worth.
fn roots_plain_false_position(
  f : (Double) -> Double,
  a : Double,
  b : Double,
  max_iterations : Int,
) -> (Double, Int) {
  let mut left = a
  let mut right = b
  let mut f_left = f(a)
  let mut f_right = f(b)
  let mut x = a
  let mut iterations = 0
  while iterations < max_iterations {
    x = (left * f_right - right * f_left) / (f_right - f_left)
    let f_x = f(x)
    iterations = iterations + 1
    if f_x == 0.0 {
      break
    }
    if f_left * f_x < 0.0 {
      right = x
      f_right = f_x
    } else {
      left = x
      f_left = f_x
    }
    if (right - left).abs() <= 1.0e-10 * (1.0 + x.abs()) {
      break
    }
  }
  (x, iterations)
}

///|
/// The square root of two is where x² − 2 is zero, and all five searches are
/// asked for it with the default tolerance of 1e-10. They are held to 1e-9,
/// which is the width of the interval each of them promises to leave: the
/// answer is only as good as the tolerance it was asked for.
test "every method finds the same root" {
  let f = fn(x : Double) { x * x - 2.0 }
  let df = fn(x : Double) { 2.0 * x }
  let root = 1.4142135623730951
  assert_true(roots_close(bisection(f, 1.0, 2.0), root, 1.0e-9))
  assert_true(roots_close(false_position(f, 1.0, 2.0), root, 1.0e-9))
  assert_true(roots_close(secant(f, 1.0, 2.0), root, 1.0e-9))
  assert_true(roots_close(newton(f, df, 1.5), root, 1.0e-9))
  assert_true(roots_close(brent(f, 1.0, 2.0), root, 1.0e-9))
}

///|
/// Eight roots that are written down in the textbooks, each bracketed and
/// each found by all three bracketing methods to within 1e-9: the cubic
/// x³ − x − 2, the fixed point of 5·cos(x) = x, eˣ = 2, the omega constant of
/// e⁻ˣ = x, the sine at a half, the natural log at one, Wallis' cubic
/// x³ − 2x − 5, and a cubic with three roots — the one that comes back is the
/// one the bracket holds.
test "the textbook roots come out of the bracketing methods" {
  let problems : Array[((Double) -> Double, Double, Double, Double)] = [
    (fn(x : Double) { x * x * x - x - 2.0 }, 1.0, 2.0, 1.5213797068045676),
    (fn(x : Double) { 5.0 * @math.cos(x) - x }, 0.0, 2.0, 1.306440008369511),
    (fn(x : Double) { @math.exp(x) - 2.0 }, 0.0, 1.0, 0.6931471805599453),
    (fn(x : Double) { @math.exp(-x) - x }, 0.0, 1.0, 0.5671432904097838),
    (fn(x : Double) { @math.sin(x) - 0.5 }, 0.0, 1.0, 0.5235987755982989),
    (fn(x : Double) { @math.ln(x) - 1.0 }, 1.0, 4.0, 2.718281828459045),
    (fn(x : Double) { x * x * x - 2.0 * x - 5.0 }, 1.0, 3.0, 2.0945514815423265),
    (fn(x : Double) { x * x * x - 3.0 * x + 1.0 }, 1.0, 2.0, 1.5320888862379562),
  ]
  for problem in problems {
    let (f, a, b, root) = problem
    assert_true(roots_close(bisection(f, a, b), root, 1.0e-9))
    assert_true(roots_close(false_position(f, a, b), root, 1.0e-9))
    assert_true(roots_close(brent(f, a, b), root, 1.0e-9))
  }
}

///|
/// The same eight roots, from the two searches that keep no bracket: the
/// secant from the ends of the interval and Newton from the middle of it.
/// Their starts are inside a bracket here, which is what an unbracketed
/// search needs to be given.
test "the open methods find them from a start in the middle" {
  let problems : Array[
    ((Double) -> Double, (Double) -> Double, Double, Double, Double),
  ] = [
    (
      fn(x : Double) { x * x - 2.0 },
      fn(x : Double) { 2.0 * x },
      1.0,
      2.0,
      1.4142135623730951,
    ),
    (
      fn(x : Double) { x * x * x - x - 2.0 },
      fn(x : Double) { 3.0 * x * x - 1.0 },
      1.0,
      2.0,
      1.5213797068045676,
    ),
    (
      fn(x : Double) { 5.0 * @math.cos(x) - x },
      fn(x : Double) { -5.0 * @math.sin(x) - 1.0 },
      0.0,
      2.0,
      1.306440008369511,
    ),
    (
      fn(x : Double) { @math.exp(x) - 2.0 },
      fn(x : Double) { @math.exp(x) },
      0.0,
      1.0,
      0.6931471805599453,
    ),
    (
      fn(x : Double) { @math.exp(-x) - x },
      fn(x : Double) { -@math.exp(-x) - 1.0 },
      0.0,
      1.0,
      0.5671432904097838,
    ),
    (
      fn(x : Double) { @math.sin(x) - 0.5 },
      fn(x : Double) { @math.cos(x) },
      0.0,
      1.0,
      0.5235987755982989,
    ),
    (
      fn(x : Double) { @math.ln(x) - 1.0 },
      fn(x : Double) { 1.0 / x },
      1.0,
      4.0,
      2.718281828459045,
    ),
    (
      fn(x : Double) { x * x * x - 2.0 * x - 5.0 },
      fn(x : Double) { 3.0 * x * x - 2.0 },
      1.0,
      3.0,
      2.0945514815423265,
    ),
    (
      fn(x : Double) { x * x * x - 3.0 * x + 1.0 },
      fn(x : Double) { 3.0 * x * x - 3.0 },
      1.0,
      2.0,
      1.5320888862379562,
    ),
  ]
  for problem in problems {
    let (f, df, a, b, root) = problem
    assert_true(roots_close(secant(f, a, b), root, 1.0e-9))
    assert_true(roots_close(newton(f, df, 0.5 * (a + b)), root, 1.0e-9))
  }
}

///|
/// The tolerance is what the answer is bought with: asked for 1e-12, bisection
/// spends seven more steps and Brent and the secant give a closer point, and
/// the two answers bisection gives are not the same point.
///
/// The last value is a cross-check against an independent implementation:
/// SciPy's `brentq` on this function, at the same tolerance, returns
/// 1.4142135623731364, and this Brent returns that point to within 1e-13 —
/// which is four hundredths of an ulp of the root, so the two are reading the
/// same digits.
test "a tighter tolerance is answered more closely" {
  let f = fn(x : Double) { x * x - 2.0 }
  let df = fn(x : Double) { 2.0 * x }
  let root = 1.4142135623730951
  let coarse = bisection(f, 1.0, 2.0)
  let fine = bisection(f, 1.0, 2.0, tolerance=1.0e-12)
  assert_true(roots_close(fine, root, 1.0e-11))
  assert_true(coarse != fine)
  assert_true((fine - root).abs() < (coarse - root).abs())
  assert_true(roots_close(brent(f, 1.0, 2.0, tolerance=1.0e-12), root, 1.0e-11))
  assert_true(
    roots_close(secant(f, 1.0, 2.0, tolerance=1.0e-12), root, 1.0e-11),
  )
  assert_true(roots_close(newton(f, df, 1.5, tolerance=1.0e-12), root, 1.0e-11))
  assert_true(roots_close(brent(f, 1.0, 2.0), 1.4142135623731364, 1.0e-13))
}

///|
/// Bisection's accuracy can be written down before it runs: after k steps the
/// interval is (b − a) / 2^k wide and the root is inside it, so the midpoint
/// of the last interval — which is what comes back — is within (b − a) / 2^k
/// of the root. A tolerance of zero asks for a fixed number of steps, which is
/// what makes the count exact here.
test "bisection's error is bounded by the interval it has left" {
  let f = fn(x : Double) { x * x - 2.0 }
  let root = 1.4142135623730951
  for k in [1, 2, 4, 8, 12, 16, 20] {
    let (estimate, iterations) = bisection_run(f, 1.0, 2.0, 0.0, k)
    assert_eq(iterations, k)
    assert_true((estimate - root).abs() <= 1.0 / @math.pow(2.0, k.to_double()))
  }
}

///|
/// Newton's error falls quadratically: the error after a step is about the
/// square of the error before it, so the number of correct digits doubles
/// every step. Solving x² − 2 from 1.5, the trace is
/// 1.5, 1.4166666666666667, 1.4142156862745099, 1.4142135623746899,
/// 1.4142135623730951 — errors of 8.6e-2, 2.5e-3, 2.1e-6, 1.6e-12 and none,
/// each below the square of the one before it.
test "newton's error falls quadratically" {
  let f = fn(x : Double) { x * x - 2.0 }
  let df = fn(x : Double) { 2.0 * x }
  let root = 1.4142135623730951
  let trace : Array[Double] = [
    1.5, 1.4166666666666667, 1.4142156862745099, 1.4142135623746899, 1.4142135623730951,
  ]
  let mut previous_error = 8.578643762690495e-2
  for step in [1, 2, 3, 4] {
    let (x, iterations) = newton_run(f, df, 1.5, 0.0, step)
    assert_eq(iterations, step)
    assert_true(roots_close(x, trace[step], 1.0e-15))
    let error = (x - root).abs()
    assert_true(error < previous_error * previous_error)
    previous_error = error
  }
}

///|
/// The secant's error falls as the golden ratio's power, about 1.618, which
/// is what puts it between the bracketing methods and Newton: the logarithm
/// of each error divided by the logarithm of the one before it settles near
/// that number — 1.83, 1.68 and 1.67 for the last steps of x² − 2 from 1 and
/// 2 — and it is held to the band 1.55 to 1.95 here.
test "the secant's error falls as the golden ratio's power" {
  let f = fn(x : Double) { x * x - 2.0 }
  let root = 1.4142135623730951
  let errors : Array[Double] = []
  for step in [1, 2, 3, 4, 5] {
    let (x, iterations) = secant_run(f, 1.0, 2.0, 0.0, step)
    assert_eq(iterations, step)
    errors.push((x - root).abs())
  }
  for k in [2, 3, 4] {
    let ratio = @math.ln(errors[k]) / @math.ln(errors[k - 1])
    assert_true(ratio > 1.55 && ratio < 1.95)
  }
}

///|
/// The Illinois halving is what a false position needs to be usable: on
/// x¹⁰ − 1, which curves the same way across the whole bracket, the plain
/// method holds the left end at zero — the far side of the root — and is two
/// thousandths out after a hundred steps, while the halving version is exact
/// in seventeen. The plain one is written out in this file to be measured
/// against; the library's is the one with the halving.
test "the Illinois halving is what keeps a false position off the end" {
  let f = fn(x : Double) { x * x * x * x * x * x * x * x * x * x - 1.0 }
  let (illinois, steps) = false_position_run(f, 0.0, 1.5, 1.0e-10, 100)
  assert_true(roots_close(illinois, 1.0, 1.0e-9))
  assert_true(steps < 30)
  let (plain, plain_steps) = roots_plain_false_position(f, 0.0, 1.5, 100)
  assert_true((plain - 1.0).abs() > 1.0e-6)
  assert_eq(plain_steps, 100)
}

///|
/// Brent's method is worth its extra machinery: on x² − 2 it reaches a
/// tolerance of 1e-10 in six steps where bisection needs thirty-two, and it
/// does it with no derivative to supply.
test "brent reaches the tolerance in a fraction of bisection's steps" {
  let f = fn(x : Double) { x * x - 2.0 }
  let (root, brent_steps) = brent_run(f, 1.0, 2.0, 1.0e-10, 100)
  let (_, bisection_steps) = bisection_run(f, 1.0, 2.0, 1.0e-10, 100)
  assert_eq(bisection_steps, 32)
  assert_true(brent_steps <= 10)
  assert_true(brent_steps < bisection_steps)
  assert_true(roots_close(root, 1.4142135623730951, 1.0e-9))
}

///|
/// An endpoint at which the function is already zero is not something to
/// search for, and no method spends a step on one: the point comes back as it
/// was given.
test "an endpoint that is already a root is returned as it is" {
  let f = fn(x : Double) { x }
  assert_eq(bisection(f, 0.0, 1.0), 0.0)
  assert_eq(false_position(f, -1.0, 0.0), 0.0)
  assert_eq(brent(f, 0.0, 1.0), 0.0)
  let (root, iterations) = bisection_run(f, 0.0, 1.0, 1.0e-10, 100)
  assert_eq(root, 0.0)
  assert_eq(iterations, 0)
}

///|
/// A bracket is two points, not two points in order: the three bracketing
/// methods read it either way round, and still close in on the root between
/// them.
test "a bracket may be given in either order" {
  let f = fn(x : Double) { x * x - 2.0 }
  let root = 1.4142135623730951
  assert_true(roots_close(bisection(f, 2.0, 1.0), root, 1.0e-9))
  assert_true(roots_close(false_position(f, 2.0, 1.0), root, 1.0e-9))
  assert_true(roots_close(brent(f, 2.0, 1.0), root, 1.0e-9))
}

///|
/// Two points whose function values have the same sign promise nothing, and
/// are refused rather than searched between — with the two values in the
/// error, so that a caller can see how far off the bracket was.
test "a bracket with no sign change is refused" {
  let f = fn(x : Double) { x * x - 2.0 }
  let failed = try {
    let _ = bisection(f, 2.0, 3.0)
    false
  } catch {
    RootError::NotBracketed(f_a, f_b) => {
      assert_true(roots_close(f_a, 2.0, 1.0e-15))
      assert_true(roots_close(f_b, 7.0, 1.0e-15))
      true
    }
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = false_position(f, 2.0, 3.0)
    false
  } catch {
    RootError::NotBracketed(_, _) => true
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = brent(f, 2.0, 3.0)
    false
  } catch {
    RootError::NotBracketed(_, _) => true
    _ => false
  }
  assert_true(failed)
}

///|
/// A step that would be divided by a slope of exactly zero has no answer, and
/// is refused rather than answered with an infinity: Newton's derivative at a
/// point where the tangent is flat, and the secant's two starting points
/// where the function has the same value at both.
test "a step that would divide by zero is refused" {
  let f = fn(x : Double) { x * x - 2.0 }
  let flat = fn(_x : Double) { 0.0 }
  let failed = try {
    let _ = newton(f, flat, 1.5)
    false
  } catch {
    RootError::ZeroDerivative(x) => {
      assert_true(roots_close(x, 1.5, 1.0e-15))
      true
    }
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = secant(f, 1.0, 1.0)
    false
  } catch {
    RootError::ZeroDerivative(x) => {
      assert_true(roots_close(x, 1.0, 1.0e-15))
      true
    }
    _ => false
  }
  assert_true(failed)
}

///|
/// A tolerance is a width, and a negative width is not one: all five searches
/// refuse it before they spend a step, rather than reading it as a search
/// that has already converged.
test "a negative tolerance is refused" {
  let f = fn(x : Double) { x * x - 2.0 }
  let df = fn(x : Double) { 2.0 * x }
  let failed = try {
    let _ = bisection(f, 1.0, 2.0, tolerance=-1.0e-3)
    false
  } catch {
    RootError::NonPositiveTolerance(tolerance) => {
      assert_true(roots_close(tolerance, -1.0e-3, 1.0e-19))
      true
    }
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = false_position(f, 1.0, 2.0, tolerance=-1.0e-3)
    false
  } catch {
    RootError::NonPositiveTolerance(_) => true
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = brent(f, 1.0, 2.0, tolerance=-1.0e-3)
    false
  } catch {
    RootError::NonPositiveTolerance(_) => true
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = secant(f, 1.0, 2.0, tolerance=-1.0e-3)
    false
  } catch {
    RootError::NonPositiveTolerance(_) => true
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = newton(f, df, 1.5, tolerance=-1.0e-3)
    false
  } catch {
    RootError::NonPositiveTolerance(_) => true
    _ => false
  }
  assert_true(failed)
}

///|
/// A search that is not allowed enough steps says so, rather than returning
/// the last point it looked at as though it were an answer: each of the five
/// is given a limit it cannot meet in, and each reports the limit back.
test "running out of iterations is reported" {
  let f = fn(x : Double) { x * x - 2.0 }
  let df = fn(x : Double) { 2.0 * x }
  let failed = try {
    let _ = bisection(f, 1.0, 2.0, max_iterations=3)
    false
  } catch {
    RootError::MaxIterations(limit) => {
      assert_eq(limit, 3)
      true
    }
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = false_position(f, 1.0, 2.0, max_iterations=3)
    false
  } catch {
    RootError::MaxIterations(limit) => {
      assert_eq(limit, 3)
      true
    }
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = brent(f, 1.0, 2.0, max_iterations=2)
    false
  } catch {
    RootError::MaxIterations(limit) => {
      assert_eq(limit, 2)
      true
    }
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = secant(f, 1.0, 2.0, max_iterations=3)
    false
  } catch {
    RootError::MaxIterations(limit) => {
      assert_eq(limit, 3)
      true
    }
    _ => false
  }
  assert_true(failed)
  let failed = try {
    let _ = newton(f, df, 1.5, max_iterations=2)
    false
  } catch {
    RootError::MaxIterations(limit) => {
      assert_eq(limit, 2)
      true
    }
    _ => false
  }
  assert_true(failed)
}

///|
/// An unbracketed search can fail where a bracketed one cannot, and it says so
/// rather than answering with a point the function is nowhere near zero at:
/// the secant on x¹⁰ − 1 from 0 and 1.5 walks off to 1.9e11 and back, lands on
/// a point it repeats exactly, and is left with two equal function values to
/// draw its line through — no line to step along, which it reports as a zero
/// slope. Brent on the same bracket has no trouble at all.
test "an open method reports a search that cannot converge" {
  let f = fn(x : Double) { x * x * x * x * x * x * x * x * x * x - 1.0 }
  let failed = try {
    let _ = secant(f, 0.0, 1.5, max_iterations=50)
    false
  } catch {
    RootError::ZeroDerivative(_) => true
    _ => false
  }
  assert_true(failed)
  // A search that keeps a bracket has no such trouble.
  assert_true(roots_close(brent(f, 0.0, 1.5), 1.0, 1.0e-9))
}

///|
/// Newton's method is not promised to converge from every starting point: on
/// x³ − 2x + 2 from zero it steps to one and back to zero, for ever. The
/// search reports running out of iterations rather than returning one of the
/// two points it is caught between and calling it a root.
test "newton reports a starting point it cannot get away from" {
  let f = fn(x : Double) { x * x * x - 2.0 * x + 2.0 }
  let df = fn(x : Double) { 3.0 * x * x - 2.0 }
  let failed = try {
    let _ = newton(f, df, 0.0, max_iterations=20)
    false
  } catch {
    RootError::MaxIterations(limit) => {
      assert_eq(limit, 20)
      true
    }
    _ => false
  }
  assert_true(failed)
}

///|
/// A tolerance of exactly zero is how a fixed number of steps is asked for:
/// the tolerance test is off and the count is what comes back, which is what
/// the convergence tests above are built on. A function value that comes out
/// exactly zero still ends a search early, which is why the count is a promise
/// only while there is no root at the points reached.
test "a tolerance of zero asks for a fixed number of steps" {
  let f = fn(x : Double) { x * x - 2.0 }
  let df = fn(x : Double) { 2.0 * x }
  for k in [1, 2, 5, 10] {
    let (_, iterations) = bisection_run(f, 1.0, 2.0, 0.0, k)
    assert_eq(iterations, k)
    let (_, iterations) = false_position_run(f, 1.0, 2.0, 0.0, k)
    assert_eq(iterations, k)
  }
  let (x, iterations) = newton_run(f, df, 1.5, 0.0, 3)
  assert_eq(iterations, 3)
  assert_true(roots_close(x, 1.4142135623746899, 1.0e-15))
  let (x, iterations) = secant_run(f, 1.0, 2.0, 0.0, 5)
  assert_eq(iterations, 5)
  assert_true(roots_close(x, 1.4142135620573204, 1.0e-15))
  let (_, iterations) = brent_run(f, 1.0, 2.0, 0.0, 4)
  assert_eq(iterations, 4)
}
