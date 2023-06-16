use std::error::Error;
use std::fmt::Write;
use std::fs;

const WIDTH: usize  = 256;
const HEIGHT: usize = 256;

fn main() -> Result<(), Box<dyn Error>> {
	println!("Hello, sailor!");

	let mut data = String::with_capacity(WIDTH * HEIGHT * 3);

	// P3 means colours are in ascii, then columns and rows,
	// then 255 for max color, then rgb triplets.
	_ = writeln!(&mut data, "P3");
	_ = writeln!(&mut data, "{WIDTH} {HEIGHT}\n255");

	for y in 0..HEIGHT {
		for x in 0..WIDTH {
			let r = (x as f32) / ((WIDTH - 1) as f32);
			let g = (y as f32) / ((HEIGHT - 1) as f32);
			let b = 0.25;

			let ir = (255.999 * r) as u8;
			let ig = (255.999 * g) as u8;
			let ib = (255.999 * b) as u8;

			_ = writeln!(&mut data, "{ir} {ig} {ib}");
		}
	}

	println!("ppm buffer capacity={} used={}", data.capacity(), data.len());

	fs::write("image.ppm", data)?;

	Ok(())
}
