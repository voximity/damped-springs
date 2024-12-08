# damped-springs

Implementation of damped springs for smooth and springy motion.

Adapted from
[Ryan Juckett's Damped Springs](https://www.ryanjuckett.com/damped-springs/).
License included in source.

## Usage

See the [examples](./examples/).

First, start by configuring your spring:

```rs
// all spring types are generic over `f32` and `f64`!
let config = SpringConfig::new(5.0, 0.5);
```

Then, pre-compute some spring math:

```rs
let params = SpringParams::from(config);
```

Create a spring or two:

```rs
let spring_x = Spring::from_equilibrium(1.0);
let spring_y = Spring::from_equilibrium(2.0);
```

Then, update your springs either on a fixed time step or on a varying interval
(if you're using a game engine):

```rs
// fixed time step of 0.1s
let time_step = SpringTimeStep::new(params, 0.1);
loop {
    spring_x.update(time_step);
    spring_y.update(time_step);
    println!("boing! ({}, {})", spring_x.position, spring_y.position);
}

// varying time step
loop {
    let time_step = SpringTimeStep::new(params, delta_time);
    spring_x.update(time_step);
    spring_y.update(time_step);
    println!("boing! ({}, {})", spring_x.position, spring_y.position);
}
```

Or, if you're only updating one spring at a varying time interval...

```rs
loop {
    spring.update_single(params, delta_time);
    // this creates the `SpringTimeStep` for you!
}
```

### Why so many types?

- A `Spring` is simply the current state of a spring that is being simulated.
- A `SpringConfig` describes the user-level parameters of (a) spring(s), like
  its angular frequency and damping ratio.
- `SpringParams` includes pre-computed coefficients to perform efficient spring
  updating later. It is pre-computed up to the point of delta time.
- Finally, `SpringTimeStep` also has pre-computed coefficients, but is now
  specific to a particular time step interval.
