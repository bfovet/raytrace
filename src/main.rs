use itertools::Itertools;
use std::{fs, io};

const MAX_VALUE: u8 = 255;

/// Write a ppm file
fn write_ppm_file(filename: &str, width: u32, height: u32) -> io::Result<()> {
    let pixels = (0..height)
        .cartesian_product(0..width)
        .map(|(y, x)| {
            let r = y as f64 / (width - 1) as f64;
            let g = x as f64 / (height - 1) as f64;
            let b = 0.0;

            format!("{} {} {}", r * 255.0, g * 255.0, b * 255.0)
        })
        .join("\n");

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
