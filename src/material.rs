use crate::hittable::HitRecord;
use crate::ray::Ray;
use approx::AbsDiffEq;
use nalgebra::Vector3;
mod reflections;
use reflections::*;
mod vectors;
use vectors::*;

#[non_exhaustive]
#[derive(Clone)]
pub enum Material {
    Lambertian { albedo: Vector3<f64> },
    Metal { albedo: Vector3<f64>, fuzz: f64 },
}

pub struct Scattered {
    pub attenuation: Vector3<f64>,
    pub scattered: Ray,
}

impl Material {
    pub(crate) fn scatter(&self, r_in: &Ray, hit_record: HitRecord) -> Option<Scattered> {
        match self {
            Material::Lambertian { albedo } => {
                let mut scatter_direction = hit_record.normal + random_unit_vector();

                // Catch degenerate scatter direction
                if scatter_direction.abs_diff_eq(&Vector3::new(0., 0., 0.), 1e-8) {
                    scatter_direction = hit_record.normal;
                }

                let scattered = Ray {
                    origin: hit_record.point,
                    direction: scatter_direction,
                };

                Some(Scattered {
                    attenuation: *albedo,
                    scattered,
                })
            }
            Material::Metal { albedo, fuzz } => {
                let reflected: Vector3<f64> =
                    reflect(&r_in.direction.normalize(), &hit_record.normal);
                let scattered = Ray {
                    origin: hit_record.point,
                    direction: reflected + *fuzz * random_unit_vector(),
                };
                if scattered.direction.dot(&hit_record.normal) > 0.0 {
                    Some(Scattered {
                        attenuation: *albedo,
                        scattered,
                    })
                } else {
                    None
                }
            }
            _ => None,
        }
    }
}
