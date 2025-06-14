use nalgebra::Vector3;

pub fn reflect(v: &Vector3<f64>, n: &Vector3<f64>) -> Vector3<f64> {
    v - 2.0 * v.dot(n) * n
}
