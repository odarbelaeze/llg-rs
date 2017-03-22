extern crate nalgebra;

use nalgebra::{Vector3, cross, norm};


#[derive(Debug)]
pub struct Integrator {
    // Gilbert damping parameter, should be 1 for critical
    damping: f64,
    // Absolute value of the gyromagnetic ratio
    giromagnetic_ratio: f64,
    // This prefactor is always used so might as well cahe it
    prefactor: f64,
}


impl Integrator {
     pub fn new(damping: f64) -> Self {
         let giromagnetic_ratio = 1.0; // 1.760859644e-11; // rad s^-1 T^-1
         Self {
             damping: damping,
             giromagnetic_ratio: giromagnetic_ratio,
             prefactor: giromagnetic_ratio / (1.0f64 + damping * damping),
         }
     }

     pub fn critical() -> Integrator {
         Self::new(1.0f64)
     }

     /// Little wrapper for that ugly delta function
     fn delta(&self, spin: &Vector3<f64>, field: &Vector3<f64>) -> Vector3<f64> {
        - self.prefactor * (
            cross(spin, field) +
            self.damping * cross(spin, &cross(spin, field))
        )
     }

     /// Performs a heun step over a spin, assuming the field
     /// independent of the spin change.
     pub fn step(&self, spin: &Vector3<f64>, field: &Vector3<f64>, time: f64) -> Vector3<f64> {
        // This is the predictor step
        let delta_spin = self.delta(spin, field);
        let spin_prima = spin + time * delta_spin;
        let spin_prima = spin_prima / norm(&spin_prima);
        let delta_spin_prima = self.delta(&spin_prima, field);
        let final_spin = spin + 0.5 * time * (delta_spin + delta_spin_prima);
        final_spin / norm(&final_spin)
     }
}


#[cfg(test)]
mod tests {
    use nalgebra::{Vector3, cross, norm, approx_eq};
    use super::{Integrator};

    #[test]
    fn cross_product_works() {
        let i = Vector3::new(1f64, 0f64, 0f64);
        let j = Vector3::new(0f64, 1f64, 0f64);
        let k = Vector3::new(0f64, 0f64, 1f64);
        assert_eq!(cross(&i, &j), k);
        assert_eq!(cross(&j, &k), i);
        assert_eq!(cross(&k, &i), j);
    }

    #[test]
    fn norm_works() {
        let i = Vector3::new(1f64, 0f64, 0f64);
        let j = Vector3::new(0f64, 1f64, 0f64);
        let k = Vector3::new(0f64, 0f64, 1f64);
        assert_eq!(norm(&i), 1.0);
        assert_eq!(norm(&j), 1.0);
        assert_eq!(norm(&k), 1.0);
        assert_eq!(norm(&(i + j)), (2.0f64).sqrt());
        assert_eq!(norm(&(i + j + k)), (3.0f64).sqrt());
    }

    #[test]
    fn normalization_works() {
        let i = Vector3::new(1f64, 0f64, 0f64);
        let not_normalized = Vector3::new(100f64, 0f64, 0f64);
        assert_eq!(not_normalized / norm(&not_normalized), i)
    }

    #[test]
    fn can_create_integrators() {
        let _integrator = Integrator::new(1.0);
        let _integrator = Integrator::critical();
    }

    #[test]
    fn can_perform_spin_step() {
        let integrator = Integrator::critical();
        let old_spin = Vector3::new(1f64, 0f64, 0f64);
        let field = Vector3::new(0f64, 0f64, 1f64);
        let delta_time = 1.0e-16;
        let _new_spin = integrator.step(&old_spin, &field, delta_time);
    }

    #[test]
    fn step_moves_spin_towards_field() {
        let integrator = Integrator::critical();
        let old_spin = Vector3::new(1f64, 0f64, 0f64);
        let field = Vector3::new(0f64, 0f64, 1f64);
        let delta_time = 1.0e-16;
        let new_spin = integrator.step(&old_spin, &field, delta_time);
        assert!(new_spin.y > 0.0);
        assert!(new_spin.z > 0.0);
    }

    #[test]
    fn many_steps_make_the_spin_go_towards_field() {
        let integrator = Integrator::critical();
        let field = Vector3::new(0f64, 0f64, 1f64);
        let delta_time = 1.0e-4;
        let mut spin = Vector3::new(1f64, 0f64, 0f64);
        for _ in 0..1000_000 {
            spin = integrator.step(&spin, &field, delta_time);
        }
        println!("{:?}", spin);
        assert!(approx_eq(&spin, &field));
    }
}
