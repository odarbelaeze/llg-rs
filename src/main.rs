extern crate llg;
extern crate nalgebra;

use llg::{HeunIntegrator};
use nalgebra::{Vector3};

fn main() {
    let integrator = HeunIntegrator::new(0.1);
    let mut spin = Vector3::new(1.0, 0.0, 0.0);
    let field = Vector3::new(0.0, 0.0, 0.05);
    for _ in 0..1_000_000 {
        println!("{} {} {}", spin.x, spin.y, spin.z);
        spin = integrator.step(&spin, &field, 1.0e-2);
    }
}