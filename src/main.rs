use itertools::Itertools;
use indicatif::ProgressIterator;
use std::{fs, io};

const MAX_VALUE: u8 = 255;

fn generate_gradient(width: u32, height: u32) -> String {
    (0..height)
        .cartesian_product(0..width)
        .progress_count(height as u64 * width as u64)
        .map(|(y, x)| {
            let r = x as f64 / (width - 1) as f64;
            let g = y as f64 / (height - 1) as f64;
            let b = 0.0;

            format!("{} {} {}", r * 255.0, g * 255.0, b * 255.0)
        })
        .join("\n")
}

/// Write a ppm file
fn write_ppm_file(filename: &str, width: u32, height: u32) -> io::Result<()> {
    let pixels = generate_gradient(width, height);

    fs::write(filename, format!(
        "P3
{width} {height}
{MAX_VALUE}
{pixels}"))?;

    Ok(())
}

fn main() {
    let filename = "image.ppm";
    let width = 256;
    let height = 256;

    if let Err(e) = write_ppm_file(filename, width, height) {
        eprintln!("Error writing PPM file: {e}");
    } else {
        println!("PPM file '{filename}' created successfully.");
    }
}
