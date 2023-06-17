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

	let mut world = World::new();
	world.objects.push(Box::new(Sphere {
		center: Point3::new(0.0, 0.0, -1.0),
		radius: 0.5,
	}));
	world.objects.push(Box::new(Sphere {
		center: Point3::new(0.0, -100.5, -1.0),
		radius: 100.0,
	}));

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
			let c   = ray_color(&world, ray);
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

fn ray_color(world: &World, ray: Ray) -> Color3 {
	if let Some(hit) = world.hit(ray, 0.0, f32::INFINITY) {
		return 0.5 * (hit.n + Vec3::broadcast(1.0));
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

#[derive(Debug, Copy, Clone)]
struct Hit {
	p: Point3,
	n: Vec3,
	t: f32,
	front_face: bool,
}

trait Hittable {
	fn hit(&self, ray: Ray, t_min: f32, t_max: f32) -> Option<Hit>;
}

struct Sphere {
	center: Point3,
	radius: f32,
}

impl Hittable for Sphere {
    fn hit(&self, ray: Ray, t_min: f32, t_max: f32) -> Option<Hit> {
		let oc     = ray.origin - self.center;
		let a      = ray.dir.mag_sq();
		let half_b = oc.dot(ray.dir);
		let c      = oc.mag_sq() - self.radius * self.radius;
		let d      = half_b * half_b - a * c;
		if d < 0.0 {
			return None;
		}

		let sqrtd = d.sqrt();

		// Nearest root in [t_min, t_max]:
		let root = (-half_b - sqrtd) / a;
		if root < t_min || t_max < root {
			let root = (-half_b + sqrtd) / a;
			if root < t_min || t_max < root {
				return None;
			}
		}

		let t = root;
		let p = ray.at(t);

		let outward_n  = (p - self.center) / self.radius;
		let front_face = ray.dir.dot(outward_n) < 0.0;

		let n = if front_face { outward_n } else { -outward_n };

		Some(Hit { p, n, t, front_face })
    }
}

struct World {
	// @Speed Replace this with a bunch of homogenous lists.
	// We actually need a nice acceleration structure here.
	objects: Vec<Box<dyn Hittable>>,
}

impl World {
	fn new() -> Self {
		Self {
			objects: Vec::new(),
		}
	}
}

impl Hittable for World {
	fn hit(&self, ray: Ray, t_min: f32, t_max: f32) -> Option<Hit> {
		let mut hit     = None;
		let mut closest = t_max;

		for o in self.objects.iter() {
			if let Some(object_hit) = o.hit(ray, t_min, closest) {
				closest = object_hit.t;
				hit     = Some(object_hit);
			}
		}

		return hit;
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
