mod spa;

extern crate assert_approx_eq;
extern crate dict_derive;
extern crate ndarray;
extern crate numpy;
extern crate pyo3;
extern crate rayon;

use ndarray::{Array2, ArrayView2};
use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::{pymodule, Py, PyModule, PyResult, Python};
use rayon::prelude::*;

struct Input {
    unixtime: u64,
    year: i16,
    month: u8,
    lat: f32,
    lon: f32,
    elev: f32,
    pressure: f32,
    temp: f32,
    atmos_refract: f32,
}

fn spa_vector(inputs: ArrayView2<f64>) -> Array2<f64> {
    // transform into vectors of inputs with correct types
    let mut v: Vec<Input> = Vec::new();
    for row in inputs.genrows() {
        v.push(Input {
            unixtime: row[0] as u64,
            year: row[1] as i16,
            month: row[2] as u8,
            lat: row[3] as f32,
            lon: row[4] as f32,
            elev: row[5] as f32,
            pressure: row[6] as f32,
            temp: row[7] as f32,
            atmos_refract: row[8] as f32,
        });
    }
    // parallel iterator
    let res_iter = v.par_iter().map(|x| {
        let delta_t = spa::calculate_deltat(&x.year, &x.month);
        let res = spa::solar_position(
            &x.unixtime,
            &x.lat,
            &x.lon,
            &x.elev,
            &x.pressure,
            &x.temp,
            &delta_t,
            &x.atmos_refract,
        );
        res
    });
    // collect results and prepare output
    let res_vec: Vec<spa::Result> = res_iter.collect();
    let mut res_array = Array2::<f64>::zeros((inputs.dim().0, 5));
    for (i, res) in res_vec.into_iter().enumerate() {
        res_array[[i, 0]] = res.theta;
        res_array[[i, 1]] = res.theta0;
        res_array[[i, 2]] = res.e;
        res_array[[i, 3]] = res.e0;
        res_array[[i, 4]] = res.phi;
    }
    res_array
}

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
        atmos_refract: f32,
    ) -> PyResult<spa::Result> {
        let delta_t = spa::calculate_deltat(&year, &month);
        let res = spa::solar_position(
            &unixtime,
            &lat,
            &lon,
            &elev,
            &pressure,
            &temp,
            &delta_t,
            &atmos_refract,
        );
        Ok(res)
    }

    #[pyfn(m, "spa_scalar")]
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
        atmos_refract: f32,
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
            atmos_refract,
        )
    }

    #[pyfn(m, "spa_vector")]
    fn spa_vector_py(py: Python, x: &PyArray2<f64>) -> Py<PyArray2<f64>> {
        let x = x.as_array();
        spa_vector(x).into_pyarray(py).to_owned()
    }

    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn test_spa_vector() {
        let inputs = vec![
            1538939951.,
            2018.,
            10.,
            30.29,
            -97.74,
            213.,
            1008.6,
            23.,
            0.5667,
        ];
        let mut two_inputs: Vec<f64> = Vec::new();
        two_inputs.extend(inputs.iter().cloned());
        two_inputs.extend(inputs.iter().cloned());
        let expected: Vec<f64> = vec![
            38.788075669522726,
            38.80104099275841,
            51.211924330477274,
            51.19895900724159,
            204.46709671428397,
        ];
        let a = Array2::from_shape_vec((2, 9), two_inputs).unwrap();
        let res = spa_vector(a.view());
        for row in res.genrows() {
            for (a, b) in row.to_vec().iter().zip(expected.iter()) {
                assert_approx_eq!(a, b);
            }
        }
    }
}
