mod spa;

extern crate pyo3;
extern crate dict_derive;

use pyo3::prelude::{pymodule, PyModule, PyResult, Python};


#[pymodule]
fn spa_rust(_py: Python, m: &PyModule) -> PyResult<()> {

  fn spa_scalar(
    unixtime: u64,
    year: i16,
    month: u8,
    lat: f32,
    lon: f32,
    elev: f32,
    pressure: f32,
    temp: f32, 
    atmos_refract: f32
  ) -> PyResult<spa::Result> {
    let delta_t = spa::calculate_deltat(&year, &month);
    let res = spa::solar_position(&unixtime, &lat, &lon, &elev, &pressure, &temp,
                                  &delta_t, &atmos_refract);
    Ok(res)
  }

  #[pyfn(m, "spa")]
  fn spa_scalar_py(
    _py: Python,
    unixtime: u64,
    year: i16,
    month: u8,
    lat: f32,
    lon: f32,
    elev: f32,
    pressure: f32,
    temp: f32, 
    atmos_refract: f32
  ) -> PyResult<spa::Result> {
    spa_scalar(
      unixtime,
      year,
      month,
      lat,
      lon,
      elev,
      pressure,
      temp, 
      atmos_refract
    )
  }

  Ok(())
}

