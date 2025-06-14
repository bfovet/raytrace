use crate::hittable::Hittable;
use crate::material::Scattered;
use lerp::Lerp;
use nalgebra::Vector3;

pub struct Ray {
    pub origin: Vector3<f64>,
    pub direction: Vector3<f64>,
}

impl Ray {
    pub(crate) fn at(&self, t: f64) -> Vector3<f64> {
        self.origin + self.direction * t
    }

    pub(crate) fn color<T>(&self, depth: i32, world: &T) -> Vector3<f64>
    where
        T: Hittable + std::marker::Sync,
    {
        if depth <= 0 {
            return Vector3::zeros();
        }
        if let Some(rec) = world.hit(self, (0.001)..f64::INFINITY) {
            if let Some(Scattered {
                attenuation,
                scattered,
            }) = rec.material.scatter(self, rec.clone())
            {
                return attenuation.component_mul(&scattered.color(depth - 1, world));
            }
            return Vector3::zeros();
        }

        let unit_direction = self.direction.normalize();
        let a = 0.5 * (unit_direction.y + 1.0);
        Vector3::new(1.0, 1.0, 1.0).lerp(Vector3::new(0.5, 0.7, 1.0), a)
    }
}
