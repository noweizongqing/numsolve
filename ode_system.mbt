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
/// The final state vector, which is the result of the run.
pub fn VectorTrajectory::end_state(self : VectorTrajectory) -> Array[Double] {
  self.states[self.states.length() - 1]
}

///|
/// The number of samples in the trajectory.
pub fn VectorTrajectory::length(self : VectorTrajectory) -> Int {
  self.times.length()
}

///|
/// The number of components in each state vector.
pub fn VectorTrajectory::dimension(self : VectorTrajectory) -> Int {
  self.states[0].length()
}

///|
/// A one-sample trajectory containing the initial state.
fn VectorTrajectory::at(t : Double, y : Array[Double]) -> VectorTrajectory {
  { times: [t], states: [y.copy()], }
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
  steps : Int,
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
  let trajectory = VectorTrajectory::at(t0, y0)
  let times = trajectory.times
  let states = trajectory.states
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
  steps : Int,
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
  steps : Int,
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
  steps : Int,
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
  mut spent : Int
  limit : Int
  dimension : Int
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
) -> (Array[Double], Double) raise OdeError {
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
  (fifth, vector_error_norm(fifth, fourth, rtol, atol))
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
  max_evaluations : Int,
) -> (VectorTrajectory, Int) raise OdeError {
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
  let trajectory = VectorTrajectory::at(t0, y0)
  let times = trajectory.times
  let states = trajectory.states
  let mut t = t0
  let mut y = y0.copy()
  let mut h = span / 100.0
  while t != t1 {
    if (t + h - t1) * span > 0.0 || t + h == t {
      h = t1 - t
    }
    let (fifth, error) = dormand_prince_system_step(
      f, t, y, h, budget, rtol, atol,
    )
    if error <= 1.0 {
      y = fifth
      t = t + h
      states.push(y.copy())
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
  rtol? : Double = 1.0e-8,
  atol? : Double = 1.0e-10,
  max_evaluations? : Int = 100000,
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
  solver? : OdeMethod = OdeMethod::RK45,
  steps? : Int = 100,
  rtol? : Double = 1.0e-8,
  atol? : Double = 1.0e-10,
  max_evaluations? : Int = 100000,
) -> VectorTrajectory raise OdeError {
  match solver {
    OdeMethod::Euler => euler_system(f, y0, t0, t1, steps)
    OdeMethod::Midpoint => midpoint_system(f, y0, t0, t1, steps)
    OdeMethod::RK4 => rk4_system(f, y0, t0, t1, steps)
    OdeMethod::RK45 =>
      rk45_system(f, y0, t0, t1, rtol~, atol~, max_evaluations~)
  }
}
