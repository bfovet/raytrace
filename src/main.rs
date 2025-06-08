use indicatif::ProgressIterator;
use itertools::Itertools;
use lerp::Lerp;
use nalgebra::Vector3;
use std::{fs, io};

const IMAGE_WIDTH: u32 = 400;
const MAX_VALUE: u8 = 255;
const ASPECT_RATIO: f64 = 16.0 / 9.0;
const IMAGE_HEIGHT: u32 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u32;
const VIEWPORT_HEIGHT: f64 = 2.0;
const VIEWPORT_WIDTH: f64 = VIEWPORT_HEIGHT * (IMAGE_WIDTH as f64 / IMAGE_HEIGHT as f64);
const FOCAL_LENGTH: f64 = 1.0;
const CAMERA_CENTER: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);

// Calculate the vectors across the horizontal and down the vertical viewport edges
const VIEWPORT_U: Vector3<f64> = Vector3::new(VIEWPORT_WIDTH, 0.0, 0.0);
const VIEWPORT_V: Vector3<f64> = Vector3::new(0.0, -VIEWPORT_HEIGHT, 0.0);

fn generate_gradient(width: u32, height: u32) -> String {
    // Calculate the horizontal and vertical delta vectors from pixel to pixel
    let pixel_delta_u = VIEWPORT_U / width as f64;
    let pixel_delta_v = VIEWPORT_V / height as f64;

    // Calculate the location of the upper left pixel
    let viewport_upper_left =
        CAMERA_CENTER - Vector3::new(0.0, 0.0, FOCAL_LENGTH) - VIEWPORT_U / 2.0 - VIEWPORT_V / 2.0;
    let pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

    (0..height)
        .cartesian_product(0..width)
        .progress_count(height as u64 * width as u64)
        .map(|(y, x)| {
            let pixel_center =
                pixel00_loc + (x as f64 * pixel_delta_u) + (y as f64 * pixel_delta_v);
            let ray_direction = pixel_center - CAMERA_CENTER;
            let ray = Ray {
                origin: CAMERA_CENTER,
                direction: ray_direction,
            };

            let pixel_color = ray.color() * MAX_VALUE as f64;

            format!(
                "{} {} {}",
                pixel_color.x as u8, pixel_color.y as u8, pixel_color.z as u8
            )
        })
        .join("\n")
}

/// Write a ppm file
fn write_ppm_file(filename: &str) -> io::Result<()> {
    let pixels = generate_gradient(IMAGE_WIDTH, IMAGE_HEIGHT);

    fs::write(
        filename,
        format!(
            "P3
{IMAGE_WIDTH} {IMAGE_HEIGHT}
{MAX_VALUE}
{pixels}"
        ),
    )?;

    Ok(())
}

fn main() {
    let filename = "image.ppm";

    if let Err(e) = write_ppm_file(filename) {
        eprintln!("Error writing PPM file: {e}");
    } else {
        println!("PPM file '{filename}' created successfully.");
    }
}

struct Ray {
    origin: Vector3<f64>,
    direction: Vector3<f64>,
}

impl Ray {
    fn at(&self, t: f64) -> Vector3<f64> {
        self.origin + self.direction * t
    }

    fn color(&self) -> Vector3<f64> {
        let t = hit_sphere(&Vector3::new(0.0, 0.0, -1.0), 0.5, self);
        if t > 0.0 {
            let n = (self.at(t) - Vector3::new(0.0, 0.0, -1.0)).normalize();
            return 0.5 * n.add_scalar(1.0);
        }

        let unit_direction = self.direction.normalize();
        let a = 0.5 * (unit_direction.y + 1.0);
        Vector3::new(1.0, 1.0, 1.0).lerp(Vector3::new(0.5, 0.7, 1.0), a)
    }
}

fn hit_sphere(center: &Vector3<f64>, radius: f64, ray: &Ray) -> f64 {
    let oc = center - ray.origin;
    let a = ray.direction.magnitude_squared();
    let h = ray.direction.dot(&oc);
    let c = oc.magnitude_squared() - radius * radius;
    let discriminant = h * h - a * c;

    if discriminant < 0.0 {
        -1.0
    } else {
        (h - discriminant.sqrt()) / a
    }
}

trait Hittable {
    fn hit(&self, ray: &Ray, ray_tmin: f64, ray_tmax: f64) -> Option<HitRecord>;
}

struct HitRecord {
    point: Vector3<f64>,
    normal: Vector3<f64>,
    t: f64,
    front_face: bool,
}

impl HitRecord {
    fn new(point: Vector3<f64>, outward_normal: Vector3<f64>, t: f64, ray: &Ray) -> Self {
        let front_face = ray.direction.dot(&outward_normal) < 0.0;
        let normal = if front_face {
            outward_normal
        } else {
            -outward_normal
        };
        HitRecord {
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

struct Sphere {
    center: Vector3<f64>,
    radius: f64,
}

impl Hittable for Sphere {
    fn hit(&self, ray: &Ray, ray_tmin: f64, ray_tmax: f64) -> Option<HitRecord> {
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
        let mut root = (h - sqrtd) / a;
        if root <= ray_tmin || ray_tmax <= root {
            root = (h - sqrtd) / a;
            if root <= ray_tmin || ray_tmax <= root {
                return None;
            }
        }

        let t = root;
        let point = ray.at(t);
        let outward_normal = (point - self.center) / self.radius;
        let rec = HitRecord::new(
            point,
            outward_normal,
            t,
            ray,
        );

        Some(rec)
    }
}

struct HittableList {
    hittables: Vec<Box<dyn Hittable>>,
}

impl HittableList {
    fn clear(&mut self) {
        self.hittables.clear();
    }
    
    fn add<T>(&mut self, hittable: T) where T: Hittable + 'static {
        self.hittables.push(Box::new(hittable));
    }
}

impl Hittable for HittableList {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let (closest, hit_record) = self.hittables.iter().fold((t_max, None), |acc, item| {
            if let Some(temp_rec) = item.hit(ray, t_min, acc.0) {
                (temp_rec.t, Some(temp_rec))
            } else {
                acc
            }
        });
        hit_record
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_gradient() {
        let pixels = generate_gradient(4, 2);
        assert_eq!(
            pixels,
            "172 205 255\n164 200 255\n164 200 255\n172 205 255\n209 227 255\n217 232 255\n217 232 255\n209 227 255"
        );
    }
}
