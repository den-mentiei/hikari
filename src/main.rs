#![allow(dead_code)]

use std::error::Error;
use std::fmt::Write;
use std::fs;
use std::time::Instant;

use ultraviolet::{Vec3, Lerp};

const ASPECT_RATIO: f32 = 16.0 / 9.0;
const WIDTH: usize      = 600;
const HEIGHT: usize     = ((WIDTH as f32) / ASPECT_RATIO) as usize;

fn main() -> Result<(), Box<dyn Error>> {
	println!("Hello, sailor!");

	// Camera

	let viewport_height = 2.0f32;
	let viewport_width  = ASPECT_RATIO * viewport_height;
	let focal_length    = 1.0f32;

	let origin     = Point3::zero();
	let horizontal = Vec3::new(viewport_width, 0.0, 0.0);
	let vertical   = Vec3::new(0.0, viewport_height, 0.0);

	let lower_left_corner = origin - horizontal / 2.0 - vertical / 2.0 - Vec3::new(0.0, 0.0, focal_length);

	// Render

	let mut data = String::with_capacity(WIDTH * HEIGHT * 3);

	// P3 means colours are in ascii, then columns and rows,
	// then 255 for max color, then rgb triplets.
	_ = writeln!(&mut data, "P3");
	_ = writeln!(&mut data, "{WIDTH} {HEIGHT}\n255");

	let t = Timer::new("ppm in-memory formatting");
	for y in (0..HEIGHT).rev() {
		for x in 0..WIDTH {
			let u   = (x as f32) / ((WIDTH - 1) as f32);
			let v   = (y as f32) / ((HEIGHT - 1) as f32);
			let dir = lower_left_corner + u * horizontal + v * vertical - origin;
			let ray = Ray::new(origin, dir);
			let c   = ray_color(ray);
			write_color(&mut data, c);
		}
	}
	drop(t);

	let t = Timer::new("output writing");
	fs::write("image.ppm", &data)?;
	drop(t);

	println!("ppm buffer capacity={} used={}", data.capacity(), data.len());

	Ok(())
}

fn ray_color(ray: Ray) -> Color3 {
	let t = hit_sphere(Point3::new(0.0, 0.0, -1.0), 0.5, ray);
	if t > 0.0 {
		// hit_point - C
		let n = (ray.at(t) - Vec3::new(0.0, 0.0, -1.0)).normalized();
		return 0.5 * (n + Vec3::broadcast(1.0));
	}
	let dir = ray.dir.normalized();
	let t   = 0.5 * (dir.y + 1.0);
	Color3::broadcast(1.0).lerp(Color3::new(0.5, 0.7, 1.0), t)
}

// (P - C)*(P - C) = r^2
// where P(t) is the ray and (C, r) is a sphere
fn hit_sphere(center: Point3, radius: f32, ray: Ray) -> f32 {
	let oc     = ray.origin - center;
	let a      = ray.dir.mag_sq();
	let half_b = oc.dot(ray.dir);
	let c      = oc.mag_sq() - radius * radius;
	let d      = half_b * half_b - a * c;
	if d < 0.0 {
		-1.0
	} else {
		(-half_b - d.sqrt()) / a
	}
}

type Point3 = Vec3;
type Color3 = Vec3;

fn write_color<W: Write>(w: &mut W, mut c: Color3) {
	c *= 255.999;
	_ = writeln!(w, "{} {} {}", c.x as u8, c.y as u8, c.z as u8);
}

struct Timer {
	label: &'static str,
	start: Instant,
}

impl Timer {
	fn new(label: &'static str) -> Self {
		Self {
			label,
			start: Instant::now(),
		}
	}
}

impl Drop for Timer {
	fn drop(&mut self) {
		let elapsed = Instant::now().duration_since(self.start);
		println!("{} took {:?}", self.label, elapsed);
	}
}

#[derive(Debug, Copy, Clone)]
struct Ray {
	origin: Point3,
	dir:    Vec3,
}

impl Ray {
	fn new(origin: Point3, dir: Vec3) -> Self {
		Self { origin, dir }
	}

	#[inline]
	fn at(&self, t: f32) -> Point3 {
		self.origin + self.dir * t
	}
}
