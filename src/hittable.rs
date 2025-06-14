use crate::material::Material;
use crate::ray::Ray;
use nalgebra::Vector3;
use std::ops::Range;

pub trait Hittable {
    fn hit(&self, ray: &Ray, interval: Range<f64>) -> Option<HitRecord>;
}

#[derive(Clone)]
pub struct HitRecord {
    pub point: Vector3<f64>,
    pub normal: Vector3<f64>,
    pub t: f64,
    pub front_face: bool,
    pub material: Material,
}

impl HitRecord {
    pub(crate) fn new(
        material: Material,
        point: Vector3<f64>,
        outward_normal: Vector3<f64>,
        t: f64,
        ray: &Ray,
    ) -> Self {
        let front_face = ray.direction.dot(&outward_normal) < 0.0;
        let normal = if front_face {
            outward_normal
        } else {
            -outward_normal
        };
        HitRecord {
            material,
            point,
            normal,
            t,
            front_face,
        }
    }

    fn set_face_normal(&mut self, ray: &Ray, outward_normal: &Vector3<f64>) {
        self.front_face = ray.direction.dot(&self.normal) < 0.0;
        self.normal = if self.front_face {
            *outward_normal
        } else {
            -*outward_normal
        }
    }
}

pub struct HittableList {
    pub hittables: Vec<Box<dyn Hittable + Sync>>,
}

impl HittableList {
    fn clear(&mut self) {
        self.hittables.clear();
    }

    pub fn add<T>(&mut self, hittable: T)
    where
        T: Hittable + 'static + Sync,
    {
        self.hittables.push(Box::new(hittable));
    }
}

impl Hittable for HittableList {
    fn hit(&self, ray: &Ray, interval: Range<f64>) -> Option<HitRecord> {
        let (_closest, hit_record) =
            self.hittables
                .iter()
                .fold((interval.end, None), |acc, item| {
                    if let Some(temp_rec) = item.hit(ray, interval.start..acc.0) {
                        (temp_rec.t, Some(temp_rec))
                    } else {
                        acc
                    }
                });
        hit_record
    }
}
