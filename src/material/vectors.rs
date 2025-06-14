use nalgebra::Vector3;

use rand::Rng;

fn random_in_unit_sphere() -> Vector3<f64> {
    let mut rng = rand::rng();
    loop {
        let vec = Vector3::new(
            rng.random_range(-1.0..1.0),
            rng.random_range(-1.0..1.0),
            rng.random_range(-1.0..1.0),
        );

        if vec.magnitude_squared() < 1.0 {
            break vec;
        }
    }
}

pub fn random_unit_vector() -> Vector3<f64> {
    random_in_unit_sphere().normalize()
}

fn random_on_hemisphere(normal: &Vector3<f64>) -> Vector3<f64> {
    let on_unit_sphere = random_unit_vector();
    if on_unit_sphere.dot(normal) > 0.0 {
        // In the same hemisphere as the normal
        on_unit_sphere
    } else {
        -on_unit_sphere
    }
}
