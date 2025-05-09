use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};

/// Write a ppm file
fn write_ppm_file(filename: &str, width: u32, height: u32) -> Result<(), Box<dyn Error>> {
    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "P3")?;
    writeln!(writer, "{width} {height}")?;
    writeln!(writer, "255")?;

    for j in 0..height {
        for i in 0..width {
            let r = i as f64 / (width - 1) as f64;
            let g = j as f64 / (height - 1) as f64;
            let b = 0.0;

            let ir = (255.999 * r) as i32;
            let ig = (255.999 * g) as i32;
            let ib = (255.999 * b) as i32;

            writeln!(writer, "{ir} {ig} {ib}")?;
        }
    }

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
