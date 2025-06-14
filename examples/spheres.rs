use nalgebra::Vector3;
use raytrace::camera::Camera;
use raytrace::hittable::HittableList;
use raytrace::material::Material;
use raytrace::sphere::Sphere;
use std::io;

fn main() -> io::Result<()> {
    let mut world = HittableList { hittables: vec![] };

    let material_ground = Material::Lambertian {
        albedo: Vector3::new(0.8, 0.8, 0.0),
    };
    let material_center = Material::Lambertian {
        albedo: Vector3::new(0.7, 0.3, 0.3),
    };
    let material_left = Material::Metal {
        albedo: Vector3::new(0.8, 0.8, 0.8),
        fuzz: 0.3,
    };
    let material_right = Material::Metal {
        albedo: Vector3::new(0.8, 0.6, 0.2),
        fuzz: 1.0,
    };

    world.add(Sphere {
        center: Vector3::new(0.0, 0.0, -1.0),
        radius: 0.5,
        material: material_ground,
    });
    world.add(Sphere {
        center: Vector3::new(0.0, -100.5, -1.0),
        radius: 100.0,
        material: material_center,
    });
    world.add(Sphere {
        center: Vector3::new(-1.0, 0.0, -1.0),
        radius: 0.5,
        material: material_left,
    });
    world.add(Sphere {
        center: Vector3::new(1.0, 0.0, -1.0),
        radius: 0.5,
        material: material_right,
    });

    let camera = Camera::new(400, 16.0 / 9.0);
    camera.render_to_disk(world)?;

    Ok(())
}
