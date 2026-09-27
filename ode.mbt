///|
/// Ordinary differential equations: one equation, `y'(t) = f(t, y)`, walked
/// forward from a starting value.
///
/// Four steppers, each handed the right-hand side, the starting value `y0`,
/// the interval `t0` to `t1`, and either a number of steps to take or a
/// tolerance to meet. All four return a `Trajectory`: the times they stopped
/// at and the values they found there, the first of them the starting point
/// itself.
///
/// The three fixed-step steppers differ only in how much of the Taylor series
/// they keep, and that is where their accuracy comes from — halving the step
/// divides the error of `euler` by two, of `midpoint` by four, and of `rk4` by
/// sixteen:
///
/// | Stepper | Order | Evaluations per step | Error when the step halves |
/// | --- | --- | --- | --- |
/// | `euler` | 1 | 1 | ÷ 2 |
/// | `midpoint` | 2 | 2 | ÷ 4 |
/// | `rk4` | 4 | 4 | ÷ 16 |
///
/// `euler` is here as the thing the others are measured against rather than as
/// a method to use: its error is large enough that `rk4` beats it at equal
/// cost by a wide margin.
///
/// `rk45` is the odd one out. It is a fifth-order stepper paired with a
/// fourth-order one, and the difference between the two answers it gets out of
/// the same seven evaluations estimates the error of the step. A step whose
/// estimate is within the tolerance is taken, the next step is sized from that
/// estimate, and the run ends with as many steps as the accuracy asked for
/// rather than as many as the caller asked for.
///
/// The tolerance of `rk45` is not the tolerance of `roots.mbt`: it is the
/// relative error a step is allowed, measured as `|y5 - y4| / (1 + |y5|)`,
/// where `y5` and `y4` are the two orders. It has to be positive — a tolerance
/// of zero can never be met, and is refused rather than run forever.
///
/// Integration backwards is allowed: `t1` may be below `t0`, and the steps
/// then have the sign of `t1 - t0`.
///
/// The steppers are scalar. A system of equations would be a `Double` holding
/// a vector, and there is no such type here.

///|
/// A trajectory that could not be run, or that ran out of its budget.
pub suberror OdeError {
  /// The number of steps is not positive, so there is no trajectory to walk.
  /// The payload is that number. The fixed-step steppers take it as a required
  /// argument rather than a defaulted one: a protocol printed at a hundred
  /// points is a hundred, and there is nothing sensible to guess.
  NonPositiveSteps(Int)
  /// The tolerance is not positive. The payload is the tolerance. The
  /// fixed-step steppers take no tolerance and never raise this; the adaptive
  /// one refuses zero as well as a negative value, because a step is only
  /// taken when its error estimate is at or below the tolerance, and an
  /// estimate of exactly zero is not something to count on.
  NonPositiveTolerance(Double)
  /// The adaptive stepper spent the evaluations it was allowed without
  /// reaching `t1`. The payload is that budget. A step costs one evaluation
  /// per stage and there are seven stages, so a budget is spent in sevens.
  MaxEvaluations(Int)
}

///|
/// Where a stepper went: the times it stopped at and the value it found at
/// each, with `times[0] == t0` and the two arrays the same length.
///
/// The times are the points the method actually used. The fixed-step steppers
/// space them evenly; `rk45` does not, and the values of its `times` are the
/// ones its step-size control chose — that is the part of its trajectory worth
/// reading.
pub struct Trajectory {
  priv times : Array[Double]
  priv values : Array[Double]
}

///|
/// The times the trajectory was sampled at, oldest first.
pub fn Trajectory::times(self : Trajectory) -> Array[Double] {
  self.times
}

///|
/// The value at each of `times`, in the same order.
pub fn Trajectory::values(self : Trajectory) -> Array[Double] {
  self.values
}

///|
/// The last value, which is the answer the run was for. A trajectory always
/// has at least one sample — the starting point — so there is one.
pub fn Trajectory::end_value(self : Trajectory) -> Double {
  self.values[self.values.length() - 1]
}

///|
/// How many samples the trajectory holds: one more than the number of steps
/// taken.
pub fn Trajectory::length(self : Trajectory) -> Int {
  self.times.length()
}

///|
/// A trajectory of one sample: the starting point, before anything is walked.
fn Trajectory::at(t : Double, y : Double) -> Trajectory {
  { times: [t], values: [y], }
}

///|
/// Walk `steps` steps of `f` from `(t0, y0)` to `t1`, a step being whatever
/// `step` computes.
///
/// The step size is `(t1 - t0) / steps`, and each step is taken from the time
/// the last one ended at rather than from `t0 + i · h`: that is the arithmetic
/// the methods are written in, and it is what the published error constants are
/// measured against. The time recorded for the last sample is `t1` itself, so
/// that a trajectory ends where it was asked to and not a rounding away from
/// it.
fn integrate(
  f : (Double, Double) -> Double,
  y0 : Double,
  t0 : Double,
  t1 : Double,
  steps : Int,
  step : ((Double, Double) -> Double, Double, Double, Double) -> Double,
) -> Trajectory raise OdeError {
  if steps <= 0 {
    raise OdeError::NonPositiveSteps(steps)
  }
  let h = (t1 - t0) / steps.to_double()
  let trajectory = Trajectory::at(t0, y0)
  let times = trajectory.times
  let values = trajectory.values
  let mut t = t0
  let mut y = y0
  let mut taken = 0
  while taken < steps {
    y = step(f, t, y, h)
    taken = taken + 1
    t = t + h
    values.push(y)
    times.push(if taken == steps { t1 } else { t })
  }
  trajectory
}

///|
/// Euler's method: follow the slope at the point you are standing on for the
/// whole step.
///
/// One evaluation per step, first-order accurate, and the method the other
/// three are improvements on.
pub fn euler(
  f : (Double, Double) -> Double,
  y0 : Double,
  t0 : Double,
  t1 : Double,
  steps : Int,
) -> Trajectory raise OdeError {
  integrate(f, y0, t0, t1, steps, fn(f, t, y, h) { y + h * f(t, y) })
}

///|
/// The midpoint method: take a Euler step, read the slope there, and use that
/// slope for the whole step instead.
///
/// Two evaluations per step for second-order accuracy: the slope is sampled in
/// the middle of the step rather than at its start, which cancels the term
/// Euler's error comes from.
pub fn midpoint(
  f : (Double, Double) -> Double,
  y0 : Double,
  t0 : Double,
  t1 : Double,
  steps : Int,
) -> Trajectory raise OdeError {
  integrate(f, y0, t0, t1, steps, fn(f, t, y, h) {
    let k1 = f(t, y)
    let k2 = f(t + 0.5 * h, y + 0.5 * h * k1)
    y + h * k2
  })
}

///|
/// The classic fourth-order Runge-Kutta: four slopes per step — at the start,
/// twice in the middle, at the end — averaged with weights `1/6`, `2/6`,
/// `2/6`, `1/6`.
///
/// The weights are Simpson's rule, which is why a right-hand side that does
/// not depend on `y` is integrated exactly: a polynomial of degree three or
/// less is what this method reproduces without error.
pub fn rk4(
  f : (Double, Double) -> Double,
  y0 : Double,
  t0 : Double,
  t1 : Double,
  steps : Int,
) -> Trajectory raise OdeError {
  integrate(f, y0, t0, t1, steps, fn(f, t, y, h) {
    let k1 = f(t, y)
    let k2 = f(t + 0.5 * h, y + 0.5 * h * k1)
    let k3 = f(t + 0.5 * h, y + 0.5 * h * k2)
    let k4 = f(t + h, y + h * k3)
    y + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
  })
}

///|
/// The evaluations of `f` a run is allowed, and how many it has spent.
///
/// The budget is checked before each evaluation rather than after a step, so
/// that a run which cannot finish stops where the budget ran out and not a
/// step later.
priv struct EvaluationBudget {
  mut spent : Int
  limit : Int
}

///|
/// One evaluation of `f`, counted, or a raise if the budget is spent.
fn EvaluationBudget::evaluate(
  self : EvaluationBudget,
  f : (Double, Double) -> Double,
  t : Double,
  y : Double,
) -> Double raise OdeError {
  if self.spent >= self.limit {
    raise OdeError::MaxEvaluations(self.limit)
  }
  self.spent = self.spent + 1
  f(t, y)
}

///|
/// The fraction of the interval the adaptive stepper opens with: small enough
/// that the first error estimate means something, and grown within a couple of
/// steps if the equation turns out to be easy.
let initial_step_fraction : Double = 100.0

///|
/// The most a step is grown by in one go, so that a misleading error estimate
/// cannot throw the run across the interval in a single step.
let maximum_growth : Double = 5.0

///|
/// The least a step is shrunk by in one go. Below it the estimate is taken to
/// be unreliable rather than the step to be too long, and the step is cut by
/// this much and tried again.
let minimum_shrink : Double = 0.2

///|
/// The margin between the error a step is allowed and the error the next step
/// is sized for: the next step aims at nine tenths of the tolerance, because a
/// step that only just met it would be rejected on the next try.
let safety_factor : Double = 0.9

///|
/// One Dormand-Prince 5(4) step: seven evaluations, two answers, and the
/// difference between them.
///
/// Returns the fifth-order answer and the error estimate `|y5 - y4|`, the
/// latter already divided by `1 + |y5|` so that a value a long way from zero
/// and a value near it are held to the same relative standard.
fn dormand_prince_step(
  f : (Double, Double) -> Double,
  t : Double,
  y : Double,
  h : Double,
  budget : EvaluationBudget,
) -> (Double, Double) raise OdeError {
  let k1 = budget.evaluate(f, t, y)
  let k2 = budget.evaluate(f, t + 1.0 / 5.0 * h, y + h * (1.0 / 5.0 * k1))
  let k3 = budget.evaluate(
    f,
    t + 3.0 / 10.0 * h,
    y + h * (3.0 / 40.0 * k1 + 9.0 / 40.0 * k2),
  )
  let k4 = budget.evaluate(
    f,
    t + 4.0 / 5.0 * h,
    y + h * (44.0 / 45.0 * k1 + -56.0 / 15.0 * k2 + 32.0 / 9.0 * k3),
  )
  let k5 = budget.evaluate(
    f,
    t + 8.0 / 9.0 * h,
    y +
    h *
    (
      19372.0 / 6561.0 * k1 +
      -25360.0 / 2187.0 * k2 +
      64448.0 / 6561.0 * k3 +
      -212.0 / 729.0 * k4
    ),
  )
  let k6 = budget.evaluate(
    f,
    t + h,
    y +
    h *
    (
      9017.0 / 3168.0 * k1 +
      -355.0 / 33.0 * k2 +
      46732.0 / 5247.0 * k3 +
      49.0 / 176.0 * k4 +
      -5103.0 / 18656.0 * k5
    ),
  )
  let k7 = budget.evaluate(
    f,
    t + h,
    y +
    h *
    (
      35.0 / 384.0 * k1 +
      500.0 / 1113.0 * k3 +
      125.0 / 192.0 * k4 +
      -2187.0 / 6784.0 * k5 +
      11.0 / 84.0 * k6
    ),
  )
  let mut fifth = y
  fifth = fifth + h * (35.0 / 384.0) * k1
  fifth = fifth + h * (500.0 / 1113.0) * k3
  fifth = fifth + h * (125.0 / 192.0) * k4
  fifth = fifth - h * (2187.0 / 6784.0) * k5
  fifth = fifth + h * (11.0 / 84.0) * k6
  let mut fourth = y
  fourth = fourth + h * (5179.0 / 57600.0) * k1
  fourth = fourth + h * (7571.0 / 16695.0) * k3
  fourth = fourth + h * (393.0 / 640.0) * k4
  fourth = fourth - h * (92097.0 / 339200.0) * k5
  fourth = fourth + h * (187.0 / 2100.0) * k6
  fourth = fourth + h * (1.0 / 40.0) * k7
  (fifth, (fifth - fourth).abs() / (1.0 + fifth.abs()))
}

///|
/// The adaptive stepper, and how many evaluations of `f` it took.
fn rk45_run(
  f : (Double, Double) -> Double,
  y0 : Double,
  t0 : Double,
  t1 : Double,
  tolerance : Double,
  max_evaluations : Int,
) -> (Trajectory, Int) raise OdeError {
  if tolerance <= 0.0 {
    raise OdeError::NonPositiveTolerance(tolerance)
  }
  let span = t1 - t0
  let budget : EvaluationBudget = { spent: 0, limit: max_evaluations, }
  let trajectory = Trajectory::at(t0, y0)
  let times = trajectory.times
  let values = trajectory.values
  let mut t = t0
  let mut y = y0
  let mut h = span / initial_step_fraction
  while t != t1 {
    // A step that would land past `t1` is cut back to it, which is also what
    // makes the last time exact. A step too small to move `t` at all is
    // replaced by the whole remaining interval rather than looping forever.
    if (t + h - t1) * span > 0.0 || t + h == t {
      h = t1 - t
    }
    let (fifth, error) = dormand_prince_step(f, t, y, h, budget)
    if error <= tolerance {
      y = fifth
      t = t + h
      values.push(y)
      times.push(t)
    }
    // The next step is `h · (tolerance / error)^(1/5)` — the power that turns a
    // fifth-order error into a length — with the safety factor and the two
    // clamps on top of it. A step that came out exact cannot be sized that way
    // and is grown as far as the clamp allows.
    let growth = if error == 0.0 {
      maximum_growth
    } else {
      let factor = safety_factor * @math.pow(tolerance / error, 1.0 / 5.0)
      if factor > maximum_growth {
        maximum_growth
      } else if factor < minimum_shrink {
        minimum_shrink
      } else {
        factor
      }
    }
    h = h * growth
  }
  (trajectory, budget.spent)
}

///|
/// The adaptive stepper: Dormand-Prince 5(4)7M, with the step size chosen by
/// the tolerance rather than given.
///
/// Seven evaluations of `f` per step produce both a fifth-order answer and a
/// fourth-order one; their difference is the error estimate, and the step is
/// taken only if that estimate, divided by `1 + |y|`, is at or below
/// `tolerance`. An accepted step is followed by one sized from the same
/// estimate, so a right-hand side that is easy in one stretch and hard in
/// another is spent steps where it is hard.
///
/// The trajectory it returns is sampled where it chose to stop, which is not
/// evenly spaced and is further apart than a fixed-step method's: for `y' = y`
/// on `[0, 1]` a tolerance of `1e-8` reaches an error of 8e-9 in eleven steps
/// and 77 evaluations, where `rk4` needs forty steps and 160 evaluations for
/// an error of the same size.
///
/// `max_evaluations` is a budget on calls to `f` rather than on steps, and each
/// step costs seven of them. Running out of it raises rather than returning a
/// trajectory that stops short of `t1`.
pub fn rk45(
  f : (Double, Double) -> Double,
  y0 : Double,
  t0 : Double,
  t1 : Double,
  tolerance? : Double = 1.0e-8,
  max_evaluations? : Int = 100000,
) -> Trajectory raise OdeError {
  rk45_run(f, y0, t0, t1, tolerance, max_evaluations).0
}
