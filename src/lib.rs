pub mod spa;
mod utils;

#[macro_use]
extern crate derive_builder;

use pyo3::prelude::*;
use pyo3::types::PyDict;

#[pymodule]
fn solpos(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(solar_position, m)?)?;
    Ok(())
}

#[pyfunction]
fn solar_position<'a>(
    _py: Python<'a>,
    unixtimes: Vec<f64>,
    request: &'a PyDict,
) -> PyResult<&'a PyDict> {
    let mut builder = spa::SolPosInputsBuilder::default();

    builder.unixtimes(unixtimes);

    if let Some(lat) = request.get_item("lat") {
        builder.latitude(lat.extract()?);
    }
    if let Some(lon) = request.get_item("lon") {
        builder.longitude(lon.extract()?);
    }
    if let Some(elev) = request.get_item("elev") {
        builder.elevation(elev.extract()?);
    } else {
        builder.elevation(0.);
    }
    if let Some(pressure) = request.get_item("pressure") {
        builder.pressure(pressure.extract()?);
    } else {
        builder.pressure(101325.);
    }
    if let Some(temp) = request.get_item("temp") {
        builder.temperature(temp.extract()?);
    } else {
        builder.temperature(12.);
    }
    if let Some(atmos_refract) = request.get_item("atmos_refract") {
        builder.atmos_refract(atmos_refract.extract()?);
    } else {
        builder.atmos_refract(0.5667);
    }
    if let Some(parallel_calcs) = request.get_item("parallel_calcs") {
        builder.parallel_calcs(parallel_calcs.extract()?);
    } else {
        builder.parallel_calcs(false);
    }
    let inputs = builder.build().unwrap(); // FIXME
    let results = _py.allow_threads(move || spa::calculate_solar_position(inputs));
    let theta_vec: Vec<f64> = results.iter().map(|r| r.theta).collect();
    let theta0_vec: Vec<f64> = results.iter().map(|r| r.theta0).collect();
    let e_vec: Vec<f64> = results.iter().map(|r| r.e).collect();
    let e0_vec: Vec<f64> = results.iter().map(|r| r.e0).collect();
    let phi_vec: Vec<f64> = results.iter().map(|r| r.phi).collect();
    let out = PyDict::new(_py);
    out.set_item("theta", theta_vec)?;
    out.set_item("theta0", theta0_vec)?;
    out.set_item("e", e_vec)?;
    out.set_item("e0", e0_vec)?;
    out.set_item("phi", phi_vec)?;
    Ok(out)
}
