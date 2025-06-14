use indicatif::ProgressIterator;
use itertools::Itertools;
use lerp::Lerp;
use nalgebra::Vector3;
use rand::Rng;
use std::ops::Range;
use std::{fs, io};

fn main() -> io::Result<()> {
    let mut world = HittableList { hittables: vec![] };

    world.add(Sphere {
        center: Vector3::new(0.0, 0.0, -1.0),
        radius: 0.5,
    });
    world.add(Sphere {
        center: Vector3::new(0.0, -100.5, -1.0),
        radius: 100.0,
    });

    let camera = Camera::new(400, 16.0 / 9.0);
    camera.render_to_disk(world)?;

    Ok(())
}

struct Camera {
    image_width: u32,
    image_height: u32,
    max_value: u8,
    aspect_ratio: f64,
    center: Vector3<f64>,
    pixel_delta_u: Vector3<f64>,
    pixel_delta_v: Vector3<f64>,
    pixel00_loc: Vector3<f64>,
    samples_per_pixel: u32,
    max_depth: u32,
}

impl Camera {
    fn new(image_width: u32, aspect_ratio: f64) -> Self {
        let max_value: u8 = 255;
        let image_height: u32 = (image_width as f64 / aspect_ratio) as u32;
        let viewport_height: f64 = 2.0;
        let viewport_width: f64 = viewport_height * (image_width as f64 / image_height as f64);
        let focal_length: f64 = 1.0;
        let center: Vector3<f64> = Vector3::zeros();

        // Calculate the vectors across the horizontal and down the vertical viewport edges.
        let viewport_u: Vector3<f64> = Vector3::new(viewport_width, 0., 0.);
        let viewport_v: Vector3<f64> = Vector3::new(0., -viewport_height, 0.);

        // Calculate the horizontal and vertical delta vectors from pixel to pixel.
        let pixel_delta_u: Vector3<f64> = viewport_u / image_width as f64;
        let pixel_delta_v: Vector3<f64> = viewport_v / image_height as f64;

        // Calculate the location of the upper left pixel.
        let viewport_upper_left: Vector3<f64> =
            center - Vector3::new(0., 0., focal_length) - viewport_u / 2. - viewport_v / 2.;
        let pixel00_loc: Vector3<f64> = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

        Self {
            image_width,
            image_height,
            max_value,
            aspect_ratio,
            center,
            pixel_delta_u,
            pixel_delta_v,
            pixel00_loc,
            samples_per_pixel: 100,
            max_depth: 50,
        }
    }

    fn get_ray(&self, i: i32, j: i32) -> Ray {
        // Get a randomly sampled camera ray for the pixel at location i,j.

        let pixel_center =
            self.pixel00_loc + (i as f64 * self.pixel_delta_u) + (j as f64 * self.pixel_delta_v);
        let pixel_sample = pixel_center + self.pixel_sample_square();

        let ray_origin = self.center;
        let ray_direction = pixel_sample - ray_origin;

        Ray {
            origin: self.center,
            direction: ray_direction,
        }
    }

    fn pixel_sample_square(&self) -> Vector3<f64> {
        let mut rng = rand::rng();
        // Returns a random point in the square surrounding a pixel at the origin.
        let px = -0.5 + rng.random::<f64>();
        let py = -0.5 + rng.random::<f64>();
        (px * self.pixel_delta_u) + (py * self.pixel_delta_v)
    }

    fn render_to_disk<T>(&self, world: T) -> io::Result<()>
    where
        T: Hittable,
    {
        let pixels = (0..self.image_height)
            .cartesian_product(0..self.image_width)
            .progress_count(self.image_height as u64 * self.image_width as u64)
            .map(|(y, x)| {
                let scale_factor = (self.samples_per_pixel as f64).recip();

                let multisampled_pixel_color = (0..self.samples_per_pixel)
                    .into_iter()
                    .map(|_| self.get_ray(x as i32, y as i32).color(self.max_depth as i32, &world) * 255.0 * scale_factor)
                    .sum::<Vector3<f64>>();

                format!(
                    "{} {} {}",
                    multisampled_pixel_color.x as u8,
                    multisampled_pixel_color.y as u8,
                    multisampled_pixel_color.z as u8
                )
            })
            .join("\n");

        fs::write(
            "output.ppm",
            format!(
                "P3
{} {}
{}
{pixels}",
                self.image_width, self.image_height, self.max_value
            ),
        )
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

    fn color<T>(&self, depth: i32, world: &T) -> Vector3<f64>
    where
        T: Hittable,
    {
        if depth <= 0 {
            return Vector3::zeros();
        }
        if let Some(rec) = world.hit(self, (0.001)..f64::INFINITY) {
            let direction = rec.normal + random_unit_vector();
            let ray = Ray {
                origin: rec.point,
                direction,
            };
            return 0.5 * ray.color(depth - 1, world);
        }

        let unit_direction = self.direction.normalize();
        let a = 0.5 * (unit_direction.y + 1.0);
        Vector3::new(1.0, 1.0, 1.0).lerp(Vector3::new(0.5, 0.7, 1.0), a)
    }
}

trait Hittable {
    fn hit(&self, ray: &Ray, interval: Range<f64>) -> Option<HitRecord>;
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
        let rec = HitRecord::new(point, outward_normal, t, ray);

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

    fn add<T>(&mut self, hittable: T)
    where
        T: Hittable + 'static,
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

fn random_in_unit_sphere() -> Vector3<f64> {
    let mut rng = rand::rng();
    loop {
        let vec = Vector3::new(rng.random_range(-1.0..1.0), rng.random_range(-1.0..1.0), rng.random_range(-1.0..1.0));

        if vec.magnitude_squared() < 1.0 {
            break vec;
        }
    }
}

fn random_unit_vector() -> Vector3<f64> {
    return random_in_unit_sphere().normalize();
}

fn random_on_hemisphere(normal: &Vector3<f64>) -> Vector3<f64> {
    let on_unit_sphere = random_unit_vector();
    if on_unit_sphere.dot(normal) > 0.0 { // In the same hemisphere as the normal
        on_unit_sphere
    } else {
        -on_unit_sphere
    }
}