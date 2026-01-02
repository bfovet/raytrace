#pragma once

#include "hittable.hpp"
#include "material.hpp"

class camera
{
public:
  double aspect_ratio = 1.0;  // Ratio of image width over height
  int image_width = 100;      // Rendered image width in pixel count
  int samples_per_pixel = 10; // Count of random samples for each pixel
  int max_depth = 10;         // Maximum number of ray bounces into scene

  double vfov = 90;                  // Vertical view angle (field of view)
  point3 lookfrom = point3(0, 0, 0); // Point camera is looking from
  point3 lookat = point3(0, 0, -1);  // Point camera is looking at
  vec3 vup = vec3(0, 1, 0);          // Camera-relative "up" direction

  void render(const hittable& world)
  {
    initialize();

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = 0; j < image_height; j++) {
      std::clog << "\rScanlines remaining: " << (image_height - j) << ' '
                << std::flush;
      for (int i = 0; i < image_width; i++) {
        color pixel_color(0, 0, 0);
        for (int sample = 0; sample < samples_per_pixel; sample++) {
          ray r = get_ray(i, j);
          pixel_color += ray_color(r, max_depth, world);
        }
        write_color(std::cout, pixel_samples_scale * pixel_color);
      }
    }

    std::clog << "\rDone.                 \n";
  }

private:
  int image_height = 0; // Rendered image height
  double pixel_samples_scale =
    0;                // Color scale factor for a sum of pixel samples
  point3 center;      // Camera center
  point3 pixel00_loc; // Location of pixel 0, 0
  vec3 pixel_delta_u; // Offset to pixel to the right
  vec3 pixel_delta_v; // Offset to pixel below
  vec3 u, v, w;       // Camera frame basis vectors

  void initialize()
  {
    image_height = static_cast<int>(image_width / aspect_ratio);
    image_height = (image_height < 1) ? 1 : image_height;

    pixel_samples_scale = 1.0 / samples_per_pixel;

    center = lookfrom;

    // Determine viewport dimensions.
    const auto focal_length = (lookfrom - lookat).length();
    const auto theta = degrees_to_radians(vfov);
    const auto h = std::tan(theta / 2);
    const auto viewport_height = 2 * h * focal_length;
    const auto viewport_width =
      viewport_height * (static_cast<double>(image_width) / image_height);

    // Calculate the u,v,w unit basis vectors for the camera coordinate frame.
    w = unit_vector(lookfrom - lookat);
    u = unit_vector(cross(vup, w));
    v = cross(w, u);

    // Calculate the vectors across the horizontal and down the vertical
    // viewport edges.
    const vec3 viewport_u =
      viewport_width * u; // Vector across viewport horizontal edge
    const vec3 viewport_v =
      viewport_height * -v; // Vector down viewport vertical edge

    // Calculate the horizontal and vertical delta vectors from pixel to pixel.
    pixel_delta_u = viewport_u / image_width;
    pixel_delta_v = viewport_v / image_height;

    // Calculate the location of the upper left pixel.
    const auto viewport_upper_left =
      center - (focal_length * w) - viewport_u / 2 - viewport_v / 2;
    pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);
  }

  [[nodiscard]] ray get_ray(int i, int j) const
  {
    // Construct a camera ray originating from the origin and directed at
    // randomly sampled point around the pixel location i, j.

    const auto offset = sample_square();
    const auto pixel_sample = pixel00_loc + ((i + offset.x()) * pixel_delta_u) +
      ((j + offset.y()) * pixel_delta_v);

    const auto ray_origin = center;
    const auto ray_direction = pixel_sample - ray_origin;

    return {ray_origin, ray_direction};
  }

  [[nodiscard]] vec3 sample_square() const
  {
    // Returns the vector to a random point in the [-.5,-.5]-[+.5,+.5] unit
    // square.
    return {random_double() - 0.5, random_double() - 0.5, 0};
  }

  [[nodiscard]] color ray_color(const ray& r,
    const int depth,
    const hittable& world) const
  {
    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0) {
      return {0, 0, 0};
    }

    if (hit_record rec; world.hit(r, interval(0.001, infinity), rec)) {
      ray scattered;
      if (color attenuation; rec.mat->scatter(r, rec, attenuation, scattered))
        return attenuation * ray_color(scattered, depth - 1, world);
      return {0, 0, 0};
    }

    const vec3 unit_direction = unit_vector(r.direction());
    const auto a = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - a) * color(1.0, 1.0, 1.0) + a * color(0.5, 0.7, 1.0);
  }
};