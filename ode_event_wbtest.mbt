///|
test "first_event locates a scalar event inside an accepted step" {
  let trajectory = rk4_system(
    fn(_t, _y : Array[Double]) { [1.0] },
    [0.0],
    0.0,
    2.0,
    100,
  ) catch {
    _ => abort("the test equation is valid")
  }
  let event = first_event(trajectory, fn(t, _y) { t - 1.0 }) catch {
    _ => abort("the event search is valid")
  }
  match event {
    Some(found) => {
      assert_true((found.time() - 1.0).abs() <= 1.0e-9)
      assert_true((found.state()[0] - 1.0).abs() <= 1.0e-8)
    }
    None => abort("the event at t = 1 must be found")
  }
}

///|
test "event direction filters a crossing" {
  let trajectory = rk4_system(
    fn(_t, _y : Array[Double]) { [1.0] },
    [0.0],
    0.0,
    2.0,
    100,
  ) catch {
    _ => abort("the test equation is valid")
  }
  let increasing = first_event(
    trajectory,
    fn(t, _y) { t - 1.0 },
    direction=EventDirection::Increasing,
  ) catch {
    _ => abort("the event search is valid")
  }
  let decreasing = first_event(
    trajectory,
    fn(t, _y) { t - 1.0 },
    direction=EventDirection::Decreasing,
  ) catch {
    _ => abort("the event search is valid")
  }
  assert_true(increasing is Some(_))
  assert_true(decreasing is None)
}

///|
test "first_event works with an adaptive oscillator trajectory" {
  let trajectory = rk45_system(
    harmonic_oscillator,
    [1.0, 0.0],
    0.0,
    3.141592653589793,
    rtol=1.0e-8,
    atol=1.0e-10,
    max_evaluations=1000000,
  ) catch {
    _ => abort("the oscillator is valid")
  }
  let event = first_event(
    trajectory,
    fn(_t, y) { y[0] },
    direction=EventDirection::Decreasing,
  ) catch {
    _ => abort("the event search is valid")
  }
  match event {
    Some(found) => {
      assert_true((found.time() - 1.5707963267948966).abs() <= 1.0e-6)
      assert_true(found.state()[0].abs() <= 1.0e-6)
      assert_true((found.state()[1] + 1.0).abs() <= 1.0e-6)
    }
    None => abort("the oscillator must cross zero at pi/2")
  }
}

///|
test "first_event returns none when no crossing exists" {
  let trajectory = rk4_system(
    fn(_t, _y : Array[Double]) { [1.0] },
    [0.0],
    0.0,
    1.0,
    20,
  ) catch {
    _ => abort("the test equation is valid")
  }
  let event = first_event(trajectory, fn(t, _y) { t - 3.0 }) catch {
    _ => abort("the event search is valid")
  }
  assert_true(event is None)
}

///|
test "first_event validates its tolerance" {
  let trajectory = rk4_system(
    fn(_t, _y : Array[Double]) { [1.0] },
    [0.0],
    0.0,
    1.0,
    20,
  ) catch {
    _ => abort("the test equation is valid")
  }
  let failed = try {
    let _ = first_event(trajectory, fn(t, _y) { t - 1.0 }, tolerance=0.0)
    false
  } catch {
    OdeError::NonPositiveTolerance(value) => value == 0.0
    _ => false
  }
  assert_true(failed)
}
