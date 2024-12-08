use damped_springs::prelude::*;

const ITERATIONS: i32 = 32;

fn main() {
    let config = SpringConfig::new(5.0, 0.5);
    let params = SpringParams::from(config);

    let mut x = Spring::from_equilibrium(1.0);
    let mut y = Spring::from_equilibrium(2.0);

    // you may need to create new `SpringTimeStep` values
    // if your delta time is not constant
    let time_step = SpringTimeStep::new(params, 0.1);
    for _ in 0..ITERATIONS {
        x.update(time_step);
        y.update(time_step);
        println!("boing! ({}, {})", x.position, y.position);
    }
}
