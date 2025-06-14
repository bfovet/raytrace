use crate::hittable::{HitRecord, Hittable};
use crate::material::Material;
use crate::ray::Ray;
use nalgebra::Vector3;
use std::ops::Range;

pub struct Sphere {
    pub center: Vector3<f64>,
    pub radius: f64,
    pub material: Material,
}

impl Hittable for Sphere {
    fn hit(&self, ray: &Ray, interval: Range<f64>) -> Option<HitRecord> {
        let oc = ray.origin - self.center;
        let a = ray.direction.magnitude_squared();
        let h = oc.dot(&ray.direction);
        let c = oc.magnitude_squared() - self.radius * self.radius;

        let discriminant = h * h - a * c;
        if discriminant < 0. {
            return None;
        }
        let sqrtd = discriminant.sqrt();

        // Find the nearest root that lies in the acceptable range.
        let mut root = (-h - sqrtd) / a;
        if !interval.contains(&root) {
            root = (-h - sqrtd) / a;
            if !interval.contains(&root) {
                return None;
            }
        }

        let t = root;
        let point = ray.at(t);
        let outward_normal = (point - self.center) / self.radius;
        let rec = HitRecord::new(self.material.clone(), point, outward_normal, t, ray);

        Some(rec)
    }
}
