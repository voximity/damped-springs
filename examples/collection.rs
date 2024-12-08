use damped_springs::prelude::*;

const ITERATIONS: i32 = 32;

fn main() {
    // A `SpringCollection` is a convenient way to manage a number of
    // identically-configured springs at once.

    let config = SpringConfig::new(5.0, 0.5);
    let mut collection = SpringCollection::from_equilibriums(config, [1.0, 2.0]);

    for _ in 0..ITERATIONS {
        collection.update(0.1);
        let [x, y] = collection.positions();
        println!("boing! ({x}, {y})");
    }
}
