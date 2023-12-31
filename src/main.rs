#![allow(dead_code)]

use std::error::Error;
use std::fmt::Write;
use std::fs;
use std::time::Instant;

use rand::prelude::*;
use ultraviolet::{Vec3, Lerp};

const ASPECT_RATIO: f32 = 3.0 / 2.0;
const WIDTH: usize      = 800;
const HEIGHT: usize     = ((WIDTH as f32) / ASPECT_RATIO) as usize;

const SAMPLES_PER_PIXEL: u32 = 500;
const MAX_DEPTH: i32         = 50;

fn main() -> Result<(), Box<dyn Error>> {
	println!("Hello, sailor!");

	let look_from      = Point3::new(13.0, 2.0, 3.0);
	let look_at        = Point3::zero();
	let up             = Vec3::new(0.0, 1.0, 0.0);
	let focus_distance = 10.0;
	let aperture       = 0.1;

	let camera = DummyCamera::new(
		look_from,
		look_at,
		up,
		20.0,
		ASPECT_RATIO,
		aperture,
		focus_distance,
	);

	let mut world = World::new();

	let material = MaterialIndex(world.materials.len());
	world.materials.push(Box::new(Lambertian::new(Color3::new(0.5, 0.5, 0.5))));
	world.objects.push(Box::new(Sphere {
		center: Point3::new(0.0, -1000.0, 0.0),
		radius: 1000.0,
		material,
	}));

	let mut rng = rand::thread_rng();
	for a in -11..11 {
		for b in -11..11 {
			let a = a as f32;
			let b = b as f32;

			let choose_material: f32 = rng.gen();
			let center = Point3::new(a + 0.9 * rng.gen::<f32>(), 0.2, b + 0.9 * rng.gen::<f32>());

			if (center - Point3::new(4.0, 0.2, 0.0)).mag_sq() > 0.9 {
				let material = MaterialIndex(world.materials.len());
				if choose_material < 0.8 {
					// Diffuse.
					let albedo   = random_color() * random_color();
					world.materials.push(Box::new(Lambertian::new(albedo)));
				} else if choose_material < 0.95 {
					// Metal.
					let albedo     = random_color();
					let rough: f32 = rng.gen_range(0.0..=0.5);
					world.materials.push(Box::new(Metal::new(albedo, rough)));
				} else {
					// Glass.
					world.materials.push(Box::new(Dielectric::new(1.5)));
				}
				world.objects.push(Box::new(Sphere {
					center,
					radius: 0.2,
					material,
				}));
			}
		}
	}

	let material = MaterialIndex(world.materials.len());
	world.materials.push(Box::new(Dielectric::new(1.5)));
	world.objects.push(Box::new(Sphere {
		center: Point3::new(0.0, 1.0, 0.0),
		radius: 1.0,
		material,
	}));
	let material = MaterialIndex(world.materials.len());
	world.materials.push(Box::new(Lambertian::new(Color3::new(0.4, 0.2, 0.1))));
	world.objects.push(Box::new(Sphere {
		center: Point3::new(-4.0, 1.0, 0.0),
		radius: 1.0,
		material,
	}));
	let material = MaterialIndex(world.materials.len());
	world.materials.push(Box::new(Metal::new(Color3::new(0.7, 0.6, 0.5), 0.0)));
	world.objects.push(Box::new(Sphere {
		center: Point3::new(4.0, 1.0, 0.0),
		radius: 1.0,
		material,
	}));

	let mut data = String::with_capacity(WIDTH * HEIGHT * 3);

	// P3 means colours are in ascii, then columns and rows,
	// then 255 for max color, then rgb triplets.
	_ = writeln!(&mut data, "P3");
	_ = writeln!(&mut data, "{WIDTH} {HEIGHT}\n255");

	let mut rng = rand::thread_rng();

	let t = Timer::new("ppm in-memory formatting");
	for y in (0..HEIGHT).rev() {
		for x in 0..WIDTH {
			let mut c = Color3::zero();
			for _ in 0..SAMPLES_PER_PIXEL {
				let dx: f32 = rng.gen();
				let dy: f32 = rng.gen();

				let u   = (x as f32 + dx) / ((WIDTH - 1) as f32);
				let v   = (y as f32 + dy) / ((HEIGHT - 1) as f32);
				let ray = camera.ray(u, v);
				c += ray_color(&world, ray, MAX_DEPTH);
			}
			write_color(&mut data, c, SAMPLES_PER_PIXEL);
		}
		// println!("rendering row #{}/{}...", y + 1, HEIGHT);
	}
	drop(t);

	let t = Timer::new("output writing");
	fs::write("image.ppm", &data)?;
	drop(t);

	println!("ppm buffer capacity={} used={}", data.capacity(), data.len());

	Ok(())
}

// @Speed Remove recursion.
fn ray_color(world: &World, ray: Ray, depth: i32) -> Color3 {
	if depth <= 0 {
		return Color3::zero();
	}

	// t_min = 0.001 to fix shadow acne.
	if let Some(hit) = world.hit(ray, 0.001, f32::INFINITY) {
		// let target = hit.p + hit.n + random_unit_vector();
		// return 0.5 * ray_color(world, Ray::new(hit.p, target - hit.p), depth - 1);
		let material = &world.materials[hit.material.0];
		if let Some((scattered, attenuation)) = material.scatter(&ray, &hit) {
			return attenuation * ray_color(world, scattered, depth - 1);
		} else {
			return Color3::zero();
		}
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

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
#[repr(transparent)]
struct MaterialIndex(usize);

#[derive(Debug, Copy, Clone)]
struct Hit {
	p: Point3,
	n: Vec3,
	t: f32,
	front_face: bool,
	material: MaterialIndex,
}

trait Hittable {
	fn hit(&self, ray: Ray, t_min: f32, t_max: f32) -> Option<Hit>;
}

trait Material {
	fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(Ray, Color3)>;
}

struct Lambertian {
	color: Color3,
}

impl Lambertian {
	fn new(color: Color3) -> Self {
		Self { color }
	}
}

impl Material for Lambertian {
	fn scatter(&self, _ray: &Ray, hit: &Hit) -> Option<(Ray, Color3)> {
		let mut scatter_dir = hit.n + random_unit_vector();
		// Degenerate scatter direction (aka opposite to the hit
		// normal).
		if is_near_zero(scatter_dir) {
			scatter_dir = hit.n;
		}
		let scattered   = Ray::new(hit.p, scatter_dir);
		Some((scattered, self.color))
	}
}

struct Metal {
	color: Color3,
	rough: f32,
}

impl Metal {
	fn new(color: Color3, rough: f32) -> Self {
		Self { color, rough }
	}
}

impl Material for Metal {
	fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(Ray, Color3)> {
		let reflected   = ray.dir.normalized().reflected(hit.n);
		let scattered   = Ray::new(hit.p, reflected + self.rough * random_in_unit_sphere());
		let attenuation = self.color;
		if scattered.dir.dot(hit.n) > 0.0 {
			Some((scattered, attenuation))
		} else {
			None
		}
	}
}

struct Dielectric {
	index_of_refraction: f32,
}

impl Dielectric {
	fn new(index_of_refraction: f32) -> Self {
		Self { index_of_refraction }
	}
}

impl Material for Dielectric {
	fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(Ray, Color3)> {
		let attenuation      = Color3::broadcast(1.0);
		let refraction_ratio = if hit.front_face {
			1.0 / self.index_of_refraction
		} else {
			self.index_of_refraction
		};


		let unit_dir     = ray.dir.normalized();
		let cos_theta    = (-unit_dir).dot(hit.n).min(1.0);
		let sin_theta    = (1.0 - cos_theta * cos_theta).sqrt();
		let cant_refract = refraction_ratio * sin_theta > 1.0;

		let mut rng   = rand::thread_rng();
		let direction = if cant_refract || reflectance(cos_theta, refraction_ratio) > rng.gen() {
			unit_dir.reflected(hit.n)
		} else {
			unit_dir.refracted(hit.n, refraction_ratio)
		};

		let scattered = Ray::new(hit.p, direction);

		Some((scattered, attenuation))
	}
}

// Schlick's approximation for reflectance.
fn reflectance(cosine: f32, ref_idx: f32) -> f32 {
	let r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
	let r0 = r0 * r0;
	r0 + (1.0 - r0) * (1.0 - cosine).powi(5)
}

struct Sphere {
	center: Point3,
	radius: f32,
	material: MaterialIndex,
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
		let mut root = (-half_b - sqrtd) / a;
		if root < t_min || t_max < root {
			root = (-half_b + sqrtd) / a;
			if root < t_min || t_max < root {
				return None;
			}
		}

		let t = root;
		let p = ray.at(t);

		let outward_n  = (p - self.center) / self.radius;
		let front_face = ray.dir.dot(outward_n) < 0.0;

		let n = if front_face { outward_n } else { -outward_n };

		let material = self.material;

		Some(Hit { p, n, t, front_face, material })
    }
}

struct World {
	// @Speed Replace this with a bunch of homogenous lists.
	// We actually need a nice acceleration structure here.
	objects:   Vec<Box<dyn Hittable>>,
	materials: Vec<Box<dyn Material>>,
}

impl World {
	fn new() -> Self {
		Self {
			objects:   Vec::new(),
			materials: Vec::new(),
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

trait Camera {
	fn ray(&self, u: f32, v: f32) -> Ray;
}

struct DummyCamera {
	origin: Point3,
	lower_left_corner: Point3,
	horizontal: Vec3,
	vertical: Vec3,
	u: Vec3,
	v: Vec3,
	w: Vec3,
	lens_radius: f32,
}

impl DummyCamera {
	fn new(
		look_from: Point3,
		look_at: Point3,
		up: Vec3,
		vertical_fov_in_degrees: f32,
		aspect_ratio: f32,
		aperture: f32,
		focus_distance: f32,
	) -> Self {
		let theta = vertical_fov_in_degrees.to_radians();
		let h     = (theta / 2.0).tan();

		let viewport_height = 2.0 * h;
		let viewport_width  = aspect_ratio * viewport_height;

		let w = (look_from - look_at).normalized();
		let u = up.cross(w).normalized();
		let v = w.cross(u);

		let origin     = look_from;
		let horizontal = focus_distance * viewport_width * u;
		let vertical   = focus_distance * viewport_height * v;

		let lower_left_corner = origin - horizontal / 2.0 - vertical / 2.0 - focus_distance * w;
		let lens_radius       = aperture / 2.0;

		Self {
			origin,
			lower_left_corner,
			horizontal,
			vertical,
			u,
			v,
			w,
			lens_radius,
		}
	}
}

impl Camera for DummyCamera {
	fn ray(&self, s: f32, t: f32) -> Ray {
		let rd     = self.lens_radius * random_in_unit_disk();
		let offset = self.u * rd.x + self.v * rd.y;
		let dir    = self.lower_left_corner + s * self.horizontal + t * self.vertical - self.origin - offset;
		Ray::new(self.origin + offset, dir)
	}
}

fn random_color() -> Color3 {
	random_color_range(0.0, 1.0)
}

fn random_color_range(lo: f32, hi: f32) -> Color3 {
	let mut rng = rand::thread_rng();
	let r = rng.gen_range(lo..=hi);
	let g = rng.gen_range(lo..=hi);
	let b = rng.gen_range(lo..=hi);
	Color3::new(r, g, b)
}

fn random_in_unit_disk() -> Vec3 {
	let mut rng = rand::thread_rng();
	loop {
		let x = rng.gen_range(-1.0..1.0);
		let y = rng.gen_range(-1.0..1.0);
		let p = Vec3::new(x, y, 0.0);
		if p.mag_sq() < 1.0 {
			return p;
		}
	}
}

type Point3 = Vec3;
type Color3 = Vec3;

fn write_color<W: Write>(w: &mut W, mut c: Color3, samples: u32) {
	let scale = 1.0 / samples as f32;
	c *= scale;

	// Gamma-correct for gamma of 2.0?
	c.x = c.x.sqrt();
	c.y = c.y.sqrt();
	c.z = c.z.sqrt();

	let min_color = Color3::broadcast(0.0);
	let max_color = Color3::broadcast(0.999);
	c.clamp(min_color, max_color);
	c *= 256.0;

	_ = writeln!(w, "{} {} {}", c.x as u8, c.y as u8, c.z as u8);
}

fn random_unit_vector() -> Vec3 {
	random_in_unit_sphere().normalized()
}

fn random_in_hemisphere(n: Vec3) -> Vec3 {
	let in_unit_sphere = random_in_unit_sphere();
	// If same hemisphere as normal.
	if in_unit_sphere.dot(n) > 0.0 {
		in_unit_sphere
	} else {
		-in_unit_sphere
	}
}

fn random_in_unit_sphere() -> Vec3 {
	let mut rng = rand::thread_rng();
	loop {
		let x = rng.gen_range(-1.0f32..1.0f32);
		let y = rng.gen_range(-1.0f32..1.0f32);
		let z = rng.gen_range(-1.0f32..1.0f32);
		let p = Vec3::new(x, y, z);
		if p.mag_sq() < 1.0 {
			return p;
		}
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

fn is_near_zero(v: Vec3) -> bool {
	const EPS: f32 = 1e-8;
	v.x.abs() < EPS && v.y.abs() < EPS && v.z.abs() < EPS
}
