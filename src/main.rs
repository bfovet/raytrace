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

            format!("{} {} {}", pixel_color.x as u8, pixel_color.y as u8, pixel_color.z as u8)
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
        let unit_direction = self.direction.normalize();
        let a = 0.5 * (unit_direction.y + 1.0);
        Vector3::new(1.0, 1.0, 1.0).lerp(Vector3::new(0.5, 0.7, 1.0), a)
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
