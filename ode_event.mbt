///|
/// Event detection over a vector ODE trajectory.
///
/// The event function is zero at the event. The detector scans the accepted
/// trajectory intervals for a sign change, filters by direction, and then
/// locates the crossing with Brent's method on the dense-output interpolation.
/// It reports the first event; it does not silently continue past it.

///|
/// Which sign changes count as an event.
pub(all) enum EventDirection {
  Any
  Increasing
  Decreasing
}

///|
/// The first event found in a trajectory.
pub struct Event {
  priv time : Double
  priv state : Array[Double]
  priv index : Int32
}

///|
/// The time at which the event occurred.
pub fn Event::time(self : Event) -> Double {
  self.time
}

///|
/// The interpolated state at the event time.
pub fn Event::state(self : Event) -> Array[Double] {
  self.state
}

///|
/// The index of the trajectory interval in which the event was found.
pub fn Event::index(self : Event) -> Int32 {
  self.index
}

///|
/// Construct an event result.
fn Event::at(time : Double, state : Array[Double], index : Int32) -> Event {
  { time, state, index, }
}

///|
/// Whether a sign change in the requested direction counts as an event.
fn event_crossed(
  left : Double,
  right : Double,
  direction : EventDirection,
) -> Bool {
  match direction {
    EventDirection::Any =>
      (left < 0.0 && right > 0.0) || (left > 0.0 && right < 0.0)
    EventDirection::Increasing => left < 0.0 && right > 0.0
    EventDirection::Decreasing => left > 0.0 && right < 0.0
  }
}

///|
/// Find the first event in a trajectory, or `None` if there is no crossing.
///
/// `event(t, y)` is evaluated at the accepted samples. A zero at a sample is
/// returned directly. A sign change is refined with Brent's method on the
/// cubic Hermite interpolation, so the event time lies inside the accepted
/// step and no extra integration step is taken.
pub fn first_event(
  trajectory : VectorTrajectory,
  event : (Double, Array[Double]) -> Double,
  direction~ : EventDirection = EventDirection::Any,
  tolerance~ : Double = 1.0e-10,
) -> Event? raise OdeError {
  if tolerance <= 0.0 {
    raise OdeError::NonPositiveTolerance(tolerance)
  }
  let times = trajectory.times()
  let states = trajectory.states()
  if times.length() < 2 {
    return None
  }
  let mut left_value = event(times[0], states[0])
  if left_value == 0.0 {
    return Some(Event::at(times[0], states[0].copy(), 0))
  }
  let mut i = 0
  while i + 1 < times.length() {
    let right_value = event(times[i + 1], states[i + 1])
    if right_value == 0.0 {
      return Some(Event::at(times[i + 1], states[i + 1].copy(), i + 1))
    }
    if event_crossed(left_value, right_value, direction) {
      let left_time = times[i]
      let right_time = times[i + 1]
      let lower = if left_time < right_time { left_time } else { right_time }
      let upper = if left_time < right_time { right_time } else { left_time }
      let event_at = fn(t : Double) -> Double {
        let state = trajectory.interpolate_at(t) catch { _ => states[i].copy() }
        event(t, state)
      }
      let root = brent(event_at, lower, upper, tolerance~, max_iterations=100) catch {
        _ => return None
      }
      let state = trajectory.interpolate_at(root) catch {
        _ => states[i].copy()
      }
      return Some(Event::at(root, state, i))
    }
    left_value = right_value
    i = i + 1
  }
  None
}
