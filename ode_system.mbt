///|
/// Ordinary differential equations with a vector state:
/// `y'(t) = f(t, y)`, where `y` is an `Array[Double]`.
///
/// The scalar steppers in `ode.mbt` remain available for one-equation
/// problems. This module adds the system form expected by general-purpose ODE
/// packages, plus a single `solve_ivp` entry point that selects the method.
/// The algorithms are standard explicit Runge-Kutta methods; this package is
/// deliberately a small pure-MoonBit implementation rather than a replacement
/// for mature ecosystems such as SciPy's `solve_ivp`.
///
/// Errors are measured with both an absolute and a relative tolerance:
/// `atol + rtol * max(|y_old|, |y_new|)` per component. The `rk45_system`
/// controller uses the root-mean-square of those scaled component errors.

///|
/// A vector ODE trajectory: the times the solver stopped at and the state at
/// each one, in matching order. The first sample is always the initial state.
pub struct VectorTrajectory {
  priv times : Array[Double]
  priv states : Array[Array[Double]]
  priv derivatives : Array[Array[Double]]
}

///|
/// The solver times, oldest first.
pub fn VectorTrajectory::times(self : VectorTrajectory) -> Array[Double] {
  self.times
}

///|
/// The state vectors at `times`, in the same order.
pub fn VectorTrajectory::states(
  self : VectorTrajectory,
) -> Array[Array[Double]] {
  self.states
}

///|
/// The derivatives at the trajectory samples. These are stored so the
/// trajectory can interpolate between accepted steps without asking the
/// caller for the right-hand side again.
pub fn VectorTrajectory::derivatives(
  self : VectorTrajectory,
) -> Array[Array[Double]] {
  self.derivatives
}

///|
/// The final state vector, which is the result of the run.
pub fn VectorTrajectory::end_state(self : VectorTrajectory) -> Array[Double] {
  self.states[self.states.length() - 1]
}

///|
/// The number of samples in the trajectory.
pub fn VectorTrajectory::length(self : VectorTrajectory) -> Int32 {
  self.times.length()
}

///|
/// The number of components in each state vector.
pub fn VectorTrajectory::dimension(self : VectorTrajectory) -> Int32 {
  self.states[0].length()
}

///|
/// A one-sample trajectory containing the initial state and its derivative.
fn VectorTrajectory::at(
  t : Double,
  y : Array[Double],
  derivative : Array[Double],
) -> VectorTrajectory {
  { times: [t], states: [y.copy()], derivatives: [derivative.copy()], }
}

///|
/// Locate the interval containing `t` in a monotonically ordered trajectory.
fn trajectory_interval(times : Array[Double], t : Double) -> Int32 {
  let last = times.length() - 1
  let ascending = times[last] >= times[0]
  let mut low = 0
  let mut high = last
  while low + 1 < high {
    let middle = (low + high) / 2
    if (ascending && times[middle] <= t) || (!ascending && times[middle] >= t) {
      low = middle
    } else {
      high = middle
    }
  }
  low
}

///|
/// Interpolate the trajectory at one time with a cubic Hermite segment.
///
/// The value and derivative at both ends of the segment are known, so the
/// cubic matches `y` and `y'` at the knots. Requesting a time outside the
/// trajectory raises `OutOfRange`; requesting a knot returns its stored state
/// without interpolation.
pub fn VectorTrajectory::interpolate_at(
  self : VectorTrajectory,
  t : Double,
) -> Array[Double] raise OdeError {
  let last = self.times.length() - 1
  let first_time = self.times[0]
  let last_time = self.times[last]
  let lower = if first_time < last_time { first_time } else { last_time }
  let upper = if first_time < last_time { last_time } else { first_time }
  if t < lower || t > upper {
    raise OdeError::OutOfRange(lower, upper, t)
  }
  if t == first_time {
    return self.states[0].copy()
  }
  if t == last_time {
    return self.states[last].copy()
  }
  let index = trajectory_interval(self.times, t)
  let left_time = self.times[index]
  let right_time = self.times[index + 1]
  let h = right_time - left_time
  if h == 0.0 {
    return self.states[index].copy()
  }
  let s = (t - left_time) / h
  let s2 = s * s
  let s3 = s2 * s
  let h00 = 2.0 * s3 - 3.0 * s2 + 1.0
  let h10 = s3 - 2.0 * s2 + s
  let h01 = -2.0 * s3 + 3.0 * s2
  let h11 = s3 - s2
  let left = self.states[index]
  let right = self.states[index + 1]
  let left_derivative = self.derivatives[index]
  let right_derivative = self.derivatives[index + 1]
  let result : Array[Double] = Array(capacity=left.length())
  let mut i = 0
  while i < left.length() {
    result.push(
      h00 * left[i] +
      h10 * h * left_derivative[i] +
      h01 * right[i] +
      h11 * h * right_derivative[i],
    )
    i = i + 1
  }
  result
}

///|
/// The derivative of the cubic Hermite segment at one time.
fn VectorTrajectory::interpolate_derivative_at(
  self : VectorTrajectory,
  t : Double,
) -> Array[Double] raise OdeError {
  let last = self.times.length() - 1
  let first_time = self.times[0]
  let last_time = self.times[last]
  let lower = if first_time < last_time { first_time } else { last_time }
  let upper = if first_time < last_time { last_time } else { first_time }
  if t < lower || t > upper {
    raise OdeError::OutOfRange(lower, upper, t)
  }
  if t == first_time {
    return self.derivatives[0].copy()
  }
  if t == last_time {
    return self.derivatives[last].copy()
  }
  let index = trajectory_interval(self.times, t)
  let left_time = self.times[index]
  let right_time = self.times[index + 1]
  let h = right_time - left_time
  if h == 0.0 {
    return self.derivatives[index].copy()
  }
  let s = (t - left_time) / h
  let s2 = s * s
  let dh00 = 6.0 * s2 - 6.0 * s
  let dh10 = 3.0 * s2 - 4.0 * s + 1.0
  let dh01 = -6.0 * s2 + 6.0 * s
  let dh11 = 3.0 * s2 - 2.0 * s
  let left = self.states[index]
  let right = self.states[index + 1]
  let left_derivative = self.derivatives[index]
  let right_derivative = self.derivatives[index + 1]
  let result : Array[Double] = Array(capacity=left.length())
  let mut i = 0
  while i < left.length() {
    result.push(
      (
        dh00 * left[i] +
        dh10 * h * left_derivative[i] +
        dh01 * right[i] +
        dh11 * h * right_derivative[i]
      ) /
      h,
    )
    i = i + 1
  }
  result
}

///|
/// Copy the first `count` samples of a trajectory.
fn VectorTrajectory::copy_prefix(
  self : VectorTrajectory,
  count : Int32,
) -> VectorTrajectory {
  let times : Array[Double] = Array(capacity=count)
  let states : Array[Array[Double]] = Array(capacity=count)
  let derivatives : Array[Array[Double]] = Array(capacity=count)
  let mut i = 0
  while i < count {
    times.push(self.times[i])
    states.push(self.states[i].copy())
    derivatives.push(self.derivatives[i].copy())
    i = i + 1
  }
  { times, states, derivatives, }
}

///|
/// Truncate a trajectory at a time inside its range.
///
/// The samples after `t` are removed and a final interpolated sample is added
/// at `t`, so the result can represent a terminal event returned by
/// `first_event`. The endpoint keeps the derivative of the Hermite segment,
/// which means it can still be sampled or truncated again.
pub fn VectorTrajectory::truncate_at(
  self : VectorTrajectory,
  t : Double,
) -> VectorTrajectory raise OdeError {
  let last = self.times.length() - 1
  let first_time = self.times[0]
  let last_time = self.times[last]
  let lower = if first_time < last_time { first_time } else { last_time }
  let upper = if first_time < last_time { last_time } else { first_time }
  if t < lower || t > upper {
    raise OdeError::OutOfRange(lower, upper, t)
  }
  if t == first_time {
    return self.copy_prefix(1)
  }
  if t == last_time {
    return self
  }
  let index = trajectory_interval(self.times, t)
  if t == self.times[index] {
    return self.copy_prefix(index + 1)
  }
  if t == self.times[index + 1] {
    return self.copy_prefix(index + 2)
  }
  let result = self.copy_prefix(index + 1)
  result.times.push(t)
  result.states.push(self.interpolate_at(t))
  result.derivatives.push(self.interpolate_derivative_at(t))
  result
}

///|
/// Sample the trajectory at each requested time, in the order requested.
/// The times need not be sorted; each one is located independently.
pub fn VectorTrajectory::sample_at(
  self : VectorTrajectory,
  sample_times : Array[Double],
) -> Array[Array[Double]] raise OdeError {
  let samples : Array[Array[Double]] = Array(capacity=sample_times.length())
  for t in sample_times {
    samples.push(self.interpolate_at(t))
  }
  samples
}

///|
/// The method used by `solve_ivp`.
pub(all) enum OdeMethod {
  Euler
  Midpoint
  RK4
  RK45
}

///|
/// Copy a state vector.
fn vector_copy(y : Array[Double]) -> Array[Double] {
  let result : Array[Double] = Array(capacity=y.length())
  let mut i = 0
  while i < y.length() {
    result.push(y[i])
    i = i + 1
  }
  result
}

///|
/// Add `scale * direction` to `target` in place.
fn vector_add_scaled(
  target : Array[Double],
  scale : Double,
  direction : Array[Double],
) -> Unit {
  let mut i = 0
  while i < target.length() {
    target[i] = target[i] + scale * direction[i]
    i = i + 1
  }
}

///|
/// `base + scale * direction`, without modifying either input.
fn vector_scale_add(
  base : Array[Double],
  scale : Double,
  direction : Array[Double],
) -> Array[Double] {
  let result = vector_copy(base)
  vector_add_scaled(result, scale, direction)
  result
}

///|
/// Evaluate the right-hand side and enforce the state dimension.
fn evaluate_system(
  f : (Double, Array[Double]) -> Array[Double],
  t : Double,
  y : Array[Double],
) -> Array[Double] raise OdeError {
  let derivative = f(t, y)
  if derivative.length() != y.length() {
    raise OdeError::DimensionMismatch(y.length(), derivative.length())
  }
  derivative
}

///|
/// Walk a fixed number of equal steps for a vector right-hand side.
fn integrate_system(
  f : (Double, Array[Double]) -> Array[Double],
  y0 : Array[Double],
  t0 : Double,
  t1 : Double,
  steps : Int32,
  step : (
    (Double, Array[Double]) -> Array[Double],
    Double,
    Array[Double],
    Double,
  ) -> Array[Double] raise OdeError,
) -> VectorTrajectory raise OdeError {
  if y0.length() == 0 {
    raise OdeError::EmptyState
  }
  if steps <= 0 {
    raise OdeError::NonPositiveSteps(steps)
  }
  let dimension = y0.length()
  let h = (t1 - t0) / steps.to_double()
  let initial_derivative = evaluate_system(f, t0, y0)
  let trajectory = VectorTrajectory::at(t0, y0, initial_derivative)
  let times = trajectory.times
  let states = trajectory.states
  let derivatives = trajectory.derivatives
  let mut t = t0
  let mut y = y0.copy()
  let mut taken = 0
  while taken < steps {
    let next = step(f, t, y, h)
    if next.length() != dimension {
      raise OdeError::DimensionMismatch(dimension, next.length())
    }
    y = next
    taken = taken + 1
    t = t + h
    states.push(y.copy())
    derivatives.push(evaluate_system(f, t, y))
    times.push(if taken == steps { t1 } else { t })
  }
  trajectory
}

///|
/// Euler's method for a vector state: one derivative evaluation per step.
pub fn euler_system(
  f : (Double, Array[Double]) -> Array[Double],
  y0 : Array[Double],
  t0 : Double,
  t1 : Double,
  steps : Int32,
) -> VectorTrajectory raise OdeError {
  integrate_system(f, y0, t0, t1, steps, fn(f, t, y, h) raise OdeError {
    vector_scale_add(y, h, evaluate_system(f, t, y))
  })
}

///|
/// The explicit midpoint method for a vector state: two derivative evaluations
/// per step and second-order accuracy.
pub fn midpoint_system(
  f : (Double, Array[Double]) -> Array[Double],
  y0 : Array[Double],
  t0 : Double,
  t1 : Double,
  steps : Int32,
) -> VectorTrajectory raise OdeError {
  integrate_system(f, y0, t0, t1, steps, fn(f, t, y, h) raise OdeError {
    let k1 = evaluate_system(f, t, y)
    let k2 = evaluate_system(f, t + 0.5 * h, vector_scale_add(y, 0.5 * h, k1))
    vector_scale_add(y, h, k2)
  })
}

///|
/// The classic fourth-order Runge-Kutta method for a vector state.
pub fn rk4_system(
  f : (Double, Array[Double]) -> Array[Double],
  y0 : Array[Double],
  t0 : Double,
  t1 : Double,
  steps : Int32,
) -> VectorTrajectory raise OdeError {
  integrate_system(f, y0, t0, t1, steps, fn(f, t, y, h) raise OdeError {
    let k1 = evaluate_system(f, t, y)
    let k2 = evaluate_system(f, t + 0.5 * h, vector_scale_add(y, 0.5 * h, k1))
    let k3 = evaluate_system(f, t + 0.5 * h, vector_scale_add(y, 0.5 * h, k2))
    let k4 = evaluate_system(f, t + h, vector_scale_add(y, h, k3))
    let result = vector_copy(y)
    vector_add_scaled(result, h / 6.0, k1)
    vector_add_scaled(result, h / 3.0, k2)
    vector_add_scaled(result, h / 3.0, k3)
    vector_add_scaled(result, h / 6.0, k4)
    result
  })
}

///|
/// The evaluations of a vector right-hand side a run is allowed.
priv struct VectorEvaluationBudget {
  mut spent : Int32
  limit : Int32
  dimension : Int32
}

///|
/// One counted vector evaluation, with the dimension checked before it is
/// accepted.
fn VectorEvaluationBudget::evaluate(
  self : VectorEvaluationBudget,
  f : (Double, Array[Double]) -> Array[Double],
  t : Double,
  y : Array[Double],
) -> Array[Double] raise OdeError {
  if self.spent >= self.limit {
    raise OdeError::MaxEvaluations(self.limit)
  }
  let derivative = f(t, y)
  if derivative.length() != self.dimension {
    raise OdeError::DimensionMismatch(self.dimension, derivative.length())
  }
  self.spent = self.spent + 1
  derivative
}

///|
/// The root-mean-square scaled error used by the vector adaptive controller.
fn vector_error_norm(
  high : Array[Double],
  low : Array[Double],
  rtol : Double,
  atol : Double,
) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < high.length() {
    let high_magnitude = high[i].abs()
    let low_magnitude = low[i].abs()
    let magnitude = if high_magnitude > low_magnitude {
      high_magnitude
    } else {
      low_magnitude
    }
    let scale = atol + rtol * magnitude
    let ratio = (high[i] - low[i]) / scale
    sum = sum + ratio * ratio
    i = i + 1
  }
  (sum / high.length().to_double()).sqrt()
}

///|
/// One Dormand-Prince 5(4) step for a vector state.
fn dormand_prince_system_step(
  f : (Double, Array[Double]) -> Array[Double],
  t : Double,
  y : Array[Double],
  h : Double,
  budget : VectorEvaluationBudget,
  rtol : Double,
  atol : Double,
) -> (Array[Double], Double, Array[Double]) raise OdeError {
  let k1 = budget.evaluate(f, t, y)
  let k2 = budget.evaluate(
    f,
    t + 1.0 / 5.0 * h,
    vector_scale_add(y, h * (1.0 / 5.0), k1),
  )
  let k3 = budget.evaluate(
    f,
    t + 3.0 / 10.0 * h,
    vector_scale_add(
      vector_scale_add(y, h * (3.0 / 40.0), k1),
      h * (9.0 / 40.0),
      k2,
    ),
  )
  let k4 = budget.evaluate(
    f,
    t + 4.0 / 5.0 * h,
    vector_scale_add(
      vector_scale_add(
        vector_scale_add(y, h * (44.0 / 45.0), k1),
        h * (-56.0 / 15.0),
        k2,
      ),
      h * (32.0 / 9.0),
      k3,
    ),
  )
  let k5 = budget.evaluate(
    f,
    t + 8.0 / 9.0 * h,
    vector_scale_add(
      vector_scale_add(
        vector_scale_add(
          vector_scale_add(y, h * (19372.0 / 6561.0), k1),
          h * (-25360.0 / 2187.0),
          k2,
        ),
        h * (64448.0 / 6561.0),
        k3,
      ),
      h * (-212.0 / 729.0),
      k4,
    ),
  )
  let k6 = budget.evaluate(
    f,
    t + h,
    vector_scale_add(
      vector_scale_add(
        vector_scale_add(
          vector_scale_add(
            vector_scale_add(y, h * (9017.0 / 3168.0), k1),
            h * (-355.0 / 33.0),
            k2,
          ),
          h * (46732.0 / 5247.0),
          k3,
        ),
        h * (49.0 / 176.0),
        k4,
      ),
      h * (-5103.0 / 18656.0),
      k5,
    ),
  )
  let k7 = budget.evaluate(
    f,
    t + h,
    vector_scale_add(
      vector_scale_add(
        vector_scale_add(
          vector_scale_add(y, h * (35.0 / 384.0), k1),
          h * (500.0 / 1113.0),
          k3,
        ),
        h * (125.0 / 192.0),
        k4,
      ),
      h * (-2187.0 / 6784.0),
      k5,
    ),
  )
  let fifth = vector_copy(y)
  vector_add_scaled(fifth, h * (35.0 / 384.0), k1)
  vector_add_scaled(fifth, h * (500.0 / 1113.0), k3)
  vector_add_scaled(fifth, h * (125.0 / 192.0), k4)
  vector_add_scaled(fifth, h * (-2187.0 / 6784.0), k5)
  vector_add_scaled(fifth, h * (11.0 / 84.0), k6)
  let fourth = vector_copy(y)
  vector_add_scaled(fourth, h * (5179.0 / 57600.0), k1)
  vector_add_scaled(fourth, h * (7571.0 / 16695.0), k3)
  vector_add_scaled(fourth, h * (393.0 / 640.0), k4)
  vector_add_scaled(fourth, h * (-92097.0 / 339200.0), k5)
  vector_add_scaled(fourth, h * (187.0 / 2100.0), k6)
  vector_add_scaled(fourth, h * (1.0 / 40.0), k7)
  (fifth, vector_error_norm(fifth, fourth, rtol, atol), k7)
}

///|
/// The adaptive vector stepper and its evaluation count.
fn rk45_system_run(
  f : (Double, Array[Double]) -> Array[Double],
  y0 : Array[Double],
  t0 : Double,
  t1 : Double,
  rtol : Double,
  atol : Double,
  max_evaluations : Int32,
) -> (VectorTrajectory, Int32) raise OdeError {
  if y0.length() == 0 {
    raise OdeError::EmptyState
  }
  if rtol <= 0.0 {
    raise OdeError::NonPositiveTolerance(rtol)
  }
  if atol <= 0.0 {
    raise OdeError::NonPositiveTolerance(atol)
  }
  let span = t1 - t0
  let budget : VectorEvaluationBudget = {
    spent: 0,
    limit: max_evaluations,
    dimension: y0.length(),
  }
  let initial_derivative = budget.evaluate(f, t0, y0)
  let trajectory = VectorTrajectory::at(t0, y0, initial_derivative)
  let times = trajectory.times
  let states = trajectory.states
  let derivatives = trajectory.derivatives
  let mut t = t0
  let mut y = y0.copy()
  let mut h = span / 100.0
  while t != t1 {
    if (t + h - t1) * span > 0.0 || t + h == t {
      h = t1 - t
    }
    let (fifth, error, derivative) = dormand_prince_system_step(
      f, t, y, h, budget, rtol, atol,
    )
    if error <= 1.0 {
      y = fifth
      t = t + h
      states.push(y.copy())
      derivatives.push(derivative)
      times.push(t)
    }
    let growth = if error == 0.0 {
      5.0
    } else {
      let factor = 0.9 * @math.pow(1.0 / error, 1.0 / 5.0)
      if factor > 5.0 {
        5.0
      } else if factor < 0.2 {
        0.2
      } else {
        factor
      }
    }
    h = h * growth
  }
  (trajectory, budget.spent)
}

///|
/// Dormand-Prince 5(4) for a vector state, with adaptive step-size control.
///
/// `rtol` and `atol` are both positive. Each component is scaled by
/// `atol + rtol * max(|y_old|, |y_new|)`, and the accepted step must have a
/// root-mean-square scaled error no greater than one.
pub fn rk45_system(
  f : (Double, Array[Double]) -> Array[Double],
  y0 : Array[Double],
  t0 : Double,
  t1 : Double,
  rtol~ : Double = 1.0e-8,
  atol~ : Double = 1.0e-10,
  max_evaluations~ : Int32 = 100000,
) -> VectorTrajectory raise OdeError {
  rk45_system_run(f, y0, t0, t1, rtol, atol, max_evaluations).0
}

///|
/// Solve a vector initial-value problem with a common, `solve_ivp`-style
/// entry point.
///
/// `method` selects a fixed-step stepper (`Euler`, `Midpoint`, `RK4`) or the
/// adaptive `RK45` method. Fixed-step methods read `steps`; `RK45` reads
/// `rtol`, `atol`, and `max_evaluations`. Integration backwards is allowed.
pub fn solve_ivp(
  f : (Double, Array[Double]) -> Array[Double],
  y0 : Array[Double],
  t0 : Double,
  t1 : Double,
  solver~ : OdeMethod = OdeMethod::RK45,
  steps~ : Int32 = 100,
  rtol~ : Double = 1.0e-8,
  atol~ : Double = 1.0e-10,
  max_evaluations~ : Int32 = 100000,
) -> VectorTrajectory raise OdeError {
  match solver {
    OdeMethod::Euler => euler_system(f, y0, t0, t1, steps)
    OdeMethod::Midpoint => midpoint_system(f, y0, t0, t1, steps)
    OdeMethod::RK4 => rk4_system(f, y0, t0, t1, steps)
    OdeMethod::RK45 =>
      rk45_system(f, y0, t0, t1, rtol~, atol~, max_evaluations~)
  }
}
