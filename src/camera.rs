use crate::hittable::Hittable;
use crate::ray::Ray;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use nalgebra::Vector3;
use rand::Rng;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use std::{fs, io};

pub struct Camera {
    pub image_width: u32,
    pub image_height: u32,
    pub max_value: u8,
    pub aspect_ratio: f64,
    pub center: Vector3<f64>,
    pub pixel_delta_u: Vector3<f64>,
    pub pixel_delta_v: Vector3<f64>,
    pub pixel00_loc: Vector3<f64>,
    pub samples_per_pixel: u32,
    pub max_depth: u32,
}

impl Camera {
    pub fn new(image_width: u32, aspect_ratio: f64) -> Self {
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

    pub fn render_to_disk<T>(&self, world: T) -> io::Result<()>
    where
        T: Hittable + std::marker::Sync,
    {
        let pixels = (0..self.image_height)
            .cartesian_product(0..self.image_width)
            .collect::<Vec<(u32, u32)>>()
            .into_par_iter()
            .progress_count(self.image_height as u64 * self.image_width as u64)
            .map(|(y, x)| {
                let scale_factor = (self.samples_per_pixel as f64).recip();

                let multisampled_pixel_color = (0..self.samples_per_pixel)
                    .into_iter()
                    .map(|_| {
                        let color = self
                            .get_ray(x as i32, y as i32)
                            .color(self.max_depth as i32, &world)
                            * 255.0
                            * scale_factor;
                        Vector3::new(
                            linear_to_gamma(color.x),
                            linear_to_gamma(color.y),
                            linear_to_gamma(color.z),
                        )
                    })
                    .sum::<Vector3<f64>>();

                format!(
                    "{} {} {}",
                    multisampled_pixel_color.x as u8,
                    multisampled_pixel_color.y as u8,
                    multisampled_pixel_color.z as u8
                )
            })
            .collect::<Vec<String>>()
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

fn linear_to_gamma(scalar: f64) -> f64 {
    scalar.sqrt()
}
