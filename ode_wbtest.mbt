///|
/// Whether two values agree to within `tolerance`.
fn ode_close(actual : Double, expected : Double, tolerance : Double) -> Bool {
  (actual - expected).abs() <= tolerance
}

///|
/// y' = y, whose solution is eᵗ: the equation every step size is measured on.
fn ode_growth(_t : Double, y : Double) -> Double {
  y
}

///|
/// y' = −y, whose solution is e⁻ᵗ: the same equation with the sign that makes
/// a stepper's error cancel rather than compound.
fn ode_decay(_t : Double, y : Double) -> Double {
  -y
}

///|
/// y' = t, whose solution is t²/2: the right-hand side does not depend on y,
/// so a Runge-Kutta step is a quadrature rule and `rk4` is Simpson's rule.
fn ode_ramp(t : Double, _y : Double) -> Double {
  t
}

///|
/// The three fixed-step steppers on y' = y over [0, 1] from 1, ten steps each,
/// against e = 2.718281828459045. Euler is the one whose error is visible in
/// the second decimal, the midpoint method in the third, and `rk4` in the
/// sixth.
test "euler, midpoint and rk4 on a doubling equation" {
  assert_true(
    ode_close(
      euler(ode_growth, 1.0, 0.0, 1.0, 10).end_value(),
      2.5937424601,
      1.0e-12,
    ),
  )
  assert_true(
    ode_close(
      midpoint(ode_growth, 1.0, 0.0, 1.0, 10).end_value(),
      2.714080846608224,
      1.0e-12,
    ),
  )
  assert_true(
    ode_close(
      rk4(ode_growth, 1.0, 0.0, 1.0, 10).end_value(),
      2.718279744135166,
      1.0e-12,
    ),
  )
  assert_true(
    ode_close(
      rk4(ode_growth, 1.0, 0.0, 1.0, 20).end_value(),
      2.718281692656335,
      1.0e-12,
    ),
  )
}

///|
/// The order of each stepper, measured rather than claimed: halving the step
/// should divide Euler's error by 2, the midpoint method's by 4 and `rk4`'s by
/// 16. The bands are wide because the measured ratios are 1.92, 3.85 and
/// 15.35 — the step is not yet small enough for the leading term alone to
/// describe the error — but they are narrow enough that a wrong order would
/// miss them.
test "halving the step divides the error by the order" {
  let exact = @math.exp(1.0)
  let ratios = [
    (
      (euler(ode_growth, 1.0, 0.0, 1.0, 10).end_value() - exact).abs(),
      (euler(ode_growth, 1.0, 0.0, 1.0, 20).end_value() - exact).abs(),
      1.8,
      2.2,
    ),
    (
      (midpoint(ode_growth, 1.0, 0.0, 1.0, 10).end_value() - exact).abs(),
      (midpoint(ode_growth, 1.0, 0.0, 1.0, 20).end_value() - exact).abs(),
      3.6,
      4.4,
    ),
    (
      (rk4(ode_growth, 1.0, 0.0, 1.0, 10).end_value() - exact).abs(),
      (rk4(ode_growth, 1.0, 0.0, 1.0, 20).end_value() - exact).abs(),
      13.0,
      19.0,
    ),
  ]
  for band in ratios {
    let (coarse, fine, low, high) = band
    let ratio = coarse / fine
    assert_true(ratio >= low && ratio <= high)
  }
}

///|
/// A right-hand side that does not depend on y is integrated by a Runge-Kutta
/// step as a quadrature rule, and Simpson's rule is exact on a cubic: y' = t
/// over [0, 1] from 0 comes out as exactly a half from `rk4`, to the last
/// rounding. Euler's rule is the rectangle rule and is not exact — 0.45
/// instead, which is what an order-one method spends to get there.
test "rk4 integrates a polynomial exactly and euler does not" {
  assert_true(
    ode_close(rk4(ode_ramp, 0.0, 0.0, 1.0, 10).end_value(), 0.5, 1.0e-15),
  )
  assert_true(
    ode_close(euler(ode_ramp, 0.0, 0.0, 1.0, 10).end_value(), 0.45, 1.0e-15),
  )
}

///|
/// The same method on an equation whose solution decays rather than grows,
/// written down so that a change to the weights shows up as this number moving
/// rather than as a difference too small to notice.
test "rk4 on a decaying equation" {
  assert_true(
    ode_close(
      rk4(ode_decay, 1.0, 0.0, 1.0, 10).end_value(),
      0.3678797744124984,
      1.0e-12,
    ),
  )
}

///|
/// A trajectory is the samples the method took: one more of them than there
/// were steps, starting at `t0` and ending at `t1`, spaced by the step size,
/// with `end_value` the last value and the two arrays the same length.
test "a trajectory is the samples the run took" {
  let trajectory = rk4(ode_growth, 1.0, 0.0, 2.0, 8)
  assert_eq(trajectory.length(), 9)
  assert_eq(trajectory.times().length(), 9)
  assert_eq(trajectory.values().length(), 9)
  assert_true(ode_close(trajectory.times()[0], 0.0, 0.0))
  assert_true(ode_close(trajectory.times()[8], 2.0, 0.0))
  assert_true(ode_close(trajectory.times()[1], 0.25, 1.0e-15))
  assert_true(ode_close(trajectory.values()[0], 1.0, 0.0))
  assert_true(ode_close(trajectory.end_value(), trajectory.values()[8], 0.0))
  let mut index = 1
  while index < 9 {
    let gap = trajectory.times()[index] - trajectory.times()[index - 1]
    assert_true(ode_close(gap, 0.25, 1.0e-15))
    index = index + 1
  }
}

///|
/// The same walk in the other direction: y' = y from (1, e) backwards to 0
/// lands on 1 again, because the equation is the same one either way and only
/// the sign of the step changes. The trajectory is then sampled from 1 down.
test "a trajectory can be walked backwards" {
  let exact = @math.exp(1.0)
  let backwards = rk4(ode_growth, exact, 1.0, 0.0, 10)
  assert_true(ode_close(backwards.end_value(), 1.0, 1.0e-5))
  assert_true(ode_close(backwards.times()[10], 0.0, 0.0))
  assert_true(ode_close(backwards.end_value(), 1.0, 1.0e-5))
  let adaptive = rk45(ode_growth, exact, 1.0, 0.0)
  assert_true(ode_close(adaptive.end_value(), 1.0, 1.0e-7))
}

///|
/// The adaptive stepper on y' = y: a tolerance of 1e-8 is met with an error of
/// 7.8e-9, in eleven steps and 77 evaluations. The step count and the
/// evaluation count are the ones the step-size control arrives at from its
/// opening step of a hundredth of the interval, so a change to the controller
/// moves them and this test says so.
test "the adaptive stepper meets the tolerance it was given" {
  let (trajectory, evaluations) = rk45_run(
    ode_growth, 1.0, 0.0, 1.0, 1.0e-8, 100000,
  )
  let exact = @math.exp(1.0)
  assert_true(ode_close(trajectory.end_value(), 2.7182818362088534, 1.0e-12))
  assert_true((trajectory.end_value() - exact).abs() <= 1.0e-8)
  assert_eq(trajectory.length(), 12)
  assert_eq(evaluations, 77)
}

///|
/// The step size the controller chose: the opening step is a hundredth of the
/// interval — small enough for the first error estimate to mean something —
/// and it is grown from there while the error stays under the tolerance. The
/// grid is not even, which is the point of an adaptive method.
test "the adaptive stepper sizes its own steps" {
  let trajectory = rk45(ode_growth, 1.0, 0.0, 1.0)
  assert_true(
    ode_close(trajectory.times()[1] - trajectory.times()[0], 0.01, 1.0e-15),
  )
  assert_true(
    ode_close(trajectory.times()[2] - trajectory.times()[1], 0.05, 1.0e-15),
  )
  assert_true(trajectory.times()[6] - trajectory.times()[5] > 0.1)
}

///|
/// A tighter tolerance is answered more closely, and the error stays under the
/// tolerance it was asked for: 1.9e-5 for 1e-4, 5.6e-7 for 1e-6, 7.8e-9 for
/// 1e-8 and 8.7e-11 for 1e-10. Meeting the tolerance is not what a step is
/// allowed to do — the step size aims at nine tenths of it — but the run as a
/// whole lands under it.
test "a tighter tolerance is answered more closely" {
  let exact = @math.exp(1.0)
  let tolerances : Array[Double] = [1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10]
  let mut previous = 1.0
  for tolerance in tolerances {
    let error = (rk45(ode_growth, 1.0, 0.0, 1.0, tolerance~).end_value() - exact).abs()
    assert_true(error < previous)
    assert_true(error <= tolerance)
    previous = error
  }
}

///|
/// What the adaptivity buys: at an error of 8e-9 the adaptive stepper spends 77
/// evaluations, where `rk4` needs forty steps — 160 evaluations — to reach the
/// same error. Both numbers are measured here, not quoted.
test "the adaptive stepper is cheaper than the fixed one" {
  let exact = @math.exp(1.0)
  let (adaptive, evaluations) = rk45_run(
    ode_growth, 1.0, 0.0, 1.0, 1.0e-8, 100000,
  )
  let adaptive_error = (adaptive.end_value() - exact).abs()
  let fixed_error = (rk4(ode_growth, 1.0, 0.0, 1.0, 40).end_value() - exact).abs()
  assert_true(fixed_error <= adaptive_error * 1.5)
  assert_true(evaluations < 160)
  assert_eq(evaluations, 77)
}

///|
/// A polynomial right-hand side is differentiated exactly, so the error
/// estimate comes out as zero and the controller grows the step as far as it
/// is allowed to: y' = t over [0, 1] is done in four steps and 28 evaluations,
/// and lands on a half to the last rounding.
test "the adaptive stepper is exact on a polynomial" {
  let (trajectory, evaluations) = rk45_run(
    ode_ramp, 0.0, 0.0, 1.0, 1.0e-8, 100000,
  )
  assert_true(ode_close(trajectory.end_value(), 0.5, 1.0e-15))
  assert_eq(trajectory.length(), 5)
  assert_eq(evaluations, 28)
}

///|
/// An equation that is hard everywhere: y' = −1000y collapses to nothing
/// within the first few thousandths of the interval, and a stepper without
/// step-size control would have to take the whole interval at the step the
/// collapse needs. The adaptive one spends its steps where the solution is
/// changing and coasts where it is not, and reaches the end of the interval.
test "a stiff equation is still walked to the end" {
  let stiff = fn(_t : Double, y : Double) { -1000.0 * y }
  let trajectory = rk45(stiff, 1.0, 0.0, 1.0)
  assert_eq(trajectory.times()[trajectory.length() - 1], 1.0)
  assert_true(trajectory.end_value().abs() <= 1.0e-6)
}

///|
/// A run of no steps has no trajectory to return, and each fixed-step stepper
/// says so with the number it was given rather than returning the starting
/// point as if it had done something.
test "a run of no steps is refused" {
  let counts : Array[Int] = [0, -3]
  for steps in counts {
    let failed = try {
      let _ = euler(ode_growth, 1.0, 0.0, 1.0, steps)
      false
    } catch {
      OdeError::NonPositiveSteps(count) => count == steps
      _ => false
    }
    assert_true(failed)
    let failed = try {
      let _ = midpoint(ode_growth, 1.0, 0.0, 1.0, steps)
      false
    } catch {
      OdeError::NonPositiveSteps(count) => count == steps
      _ => false
    }
    assert_true(failed)
    let failed = try {
      let _ = rk4(ode_growth, 1.0, 0.0, 1.0, steps)
      false
    } catch {
      OdeError::NonPositiveSteps(count) => count == steps
      _ => false
    }
    assert_true(failed)
  }
}

///|
/// The adaptive stepper has no fixed number of steps to be given, and a
/// tolerance it cannot meet — zero, or a negative one — is refused rather than
/// run until it gives up. This is the one place the tolerance of this module
/// differs from the tolerance of `roots.mbt`, where zero means "take the fixed
/// number of steps".
test "an adaptive run needs a positive tolerance" {
  let tolerances : Array[Double] = [0.0, -1.0e-3]
  for tolerance in tolerances {
    let failed = try {
      let _ = rk45(ode_growth, 1.0, 0.0, 1.0, tolerance~)
      false
    } catch {
      OdeError::NonPositiveTolerance(refused) => refused == tolerance
      _ => false
    }
    assert_true(failed)
  }
}

///|
/// The budget is on evaluations of the right-hand side, and each step costs
/// seven of them: a budget of ten cannot finish even two steps, and the run
/// stops with the budget in the error rather than returning a trajectory that
/// stops short of `t1`.
test "a budget that runs out is refused" {
  let failed = try {
    let _ = rk45_run(ode_growth, 1.0, 0.0, 1.0, 1.0e-8, 10)
    false
  } catch {
    OdeError::MaxEvaluations(limit) => limit == 10
    _ => false
  }
  assert_true(failed)
}

///|
/// The public `rk45` passes its optional arguments through: the default
/// tolerance is the one the documentation states, and a labelled call that
/// names a smaller budget runs out of it.
test "the adaptive stepper takes its arguments by label" {
  assert_true(
    ode_close(
      rk45(ode_growth, 1.0, 0.0, 1.0).end_value(),
      2.7182818362088534,
      1.0e-12,
    ),
  )
  let failed = try {
    let _ = rk45(ode_growth, 1.0, 0.0, 1.0, max_evaluations=10)
    false
  } catch {
    OdeError::MaxEvaluations(limit) => limit == 10
    _ => false
  }
  assert_true(failed)
}
