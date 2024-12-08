use damped_springs::prelude::*;

const ITERATIONS: i32 = 32;

fn main() {
    let config = SpringConfig::new(5.0, 0.5);
    let params = SpringParams::from(config);

    // manually-created time step values
    // do this when:
    // - your delta time is constant, or
    // - you are updating multiple spring simultaneously
    println!("manual time step");
    let time_step = SpringTimeStep::new(params, 0.1);
    let mut spring = Spring::from_equilibrium(1.0);
    for _ in 0..ITERATIONS {
        spring.update(time_step);
        println!("boing! {}", spring.position);
    }

    // automatically-created time step values
    // - do this when your delta time is variable, and
    // - you are only updating a single spring
    println!("auto time step");
    let mut spring = Spring::from_equilibrium(1.0);
    for _ in 0..ITERATIONS {
        spring.update_single(params, 0.1);
        println!("boing! {}", spring.position);
    }
}
