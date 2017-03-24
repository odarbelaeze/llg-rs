extern crate nalgebra;

use nalgebra::{Vector3, cross, norm};

#[derive(Debug)]
pub struct HeunIntegrator {
    gamma: f64,
    lambda: f64,
    prefactor: f64,
}


impl HeunIntegrator {
    pub fn new(lambda: f64) -> HeunIntegrator {
        HeunIntegrator {
            gamma: 1.0,
            lambda: lambda,
            prefactor: 1.0 / (1.0 + lambda * lambda)
        }
    }

    fn delta(&self, spin: &Vector3<f64>, field: &Vector3<f64>) -> Vector3<f64> {
        - self.prefactor * (
            cross(spin, field) +
            self.lambda * &cross(spin, &cross(spin, field))
        )
    }

    pub fn step(&self, spin: &Vector3<f64>, field: &Vector3<f64>, time: f64) -> Vector3<f64> {
        let delta = self.delta(spin, field);
        let spin_prima = spin + delta * time;
        let spin_prima = spin_prima / norm(&spin_prima);
        // I should have computed field_prima
        let delta_prima = self.delta(&spin_prima, field);
        let new_spin = spin + 0.5 * (delta + delta_prima) * time;
        new_spin / norm(&new_spin)
    }
}


#[cfg(test)]
mod test {
    use nalgebra::{Vector3, cross, norm};
    use super::{HeunIntegrator};

    #[test]
    fn cross_product_works() {
        let i = Vector3::new(1.0, 0.0, 0.0);
        let j = Vector3::new(0.0, 1.0, 0.0);
        let k = Vector3::new(0.0, 0.0, 1.0);
        assert_eq!(cross(&i, &j), k);
        assert_eq!(cross(&j, &k), i);
    }

    #[test]
    fn norm_works_as_expected() {
        assert_eq!(norm(&Vector3::new(0.0, 0.0, 1.0)), 1.0);
        assert_eq!(norm(&Vector3::new(0.0, 1.0, 1.0)), (2.0f64).sqrt());
    }

    #[test]
    fn can_create_heun_integrators() {
        let _integrator = HeunIntegrator::new(1.0);
    }

    #[test]
    fn heun_integrator_is_norm_conserving() {
        let integrator = HeunIntegrator::new(1.0);
        let initial_spin = Vector3::new(1.0, 0.0, 0.0);
        let field = Vector3::new(0.0, 0.0, 1.0);
        let new_spin = integrator.step(&initial_spin, &field, 1.0e-5);
        assert_eq!(norm(&new_spin), 1.0);
    }

    #[test]
    fn heun_integrator_advances_spin() {
        let integrator = HeunIntegrator::new(1.0);
        let initial_spin = Vector3::new(1.0, 0.0, 0.0);
        let field = Vector3::new(0.0, 0.0, 1.0);
        let new_spin = integrator.step(&initial_spin, &field, 1.0e-5);
        assert!(new_spin.x > 0.0);
        assert!(new_spin.y > 0.0);
    }
}