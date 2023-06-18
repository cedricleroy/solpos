mod tables;

use crate::utils::get_year_and_month_from_unixtime;
use rayon::prelude::*;

const SECONDS_PER_DAY: f64 = 86400.0;
const DAYS_PER_CENTURY: f64 = 36525.0;
const JD_CONSTANT: f64 = 2440587.5;
const JC_CONSTANT: f64 = 2451545.0;
const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.;
const RAD_TO_DEG: f64 = 180. / std::f64::consts::PI;

pub struct SolPosResult {
    pub theta: f64,
    pub theta0: f64,
    pub e: f64,
    pub e0: f64,
    pub phi: f64,
}

#[derive(Default, Builder, Clone)]
pub struct SolPosInputs {
    pub unixtimes: Vec<f64>,
    pub latitude: f64,
    pub longitude: f64,
    pub elevation: f64,
    pub pressure: f64,
    pub temperature: f64,
    pub atmos_refract: f64,
    pub parallel_calcs: bool,
}

pub fn calculate_solar_position(inputs: SolPosInputs) -> Vec<SolPosResult> {
    let calculate = |u: &f64| {
        let year_month = get_year_and_month_from_unixtime(*u as i64);
        let delta_t = calculate_deltat(year_month.0, year_month.1);
        solar_position(
            *u,
            inputs.latitude,
            inputs.longitude,
            inputs.elevation,
            inputs.pressure,
            inputs.temperature,
            delta_t,
            inputs.atmos_refract,
        )
    };
    let results = match inputs.parallel_calcs {
        true => inputs.unixtimes.par_iter().map(calculate).collect(),
        false => inputs.unixtimes.iter().map(calculate).collect(),
    };
    results
}

pub fn solar_position(
    unixtime: f64,
    latitude: f64,
    longitude: f64,
    elevation: f64,
    pressure: f64,
    temperature: f64,
    delta_t: f64,
    atmos_refract: f64,
) -> SolPosResult {
    let jd = julian_day(unixtime);
    let jde = julian_ephemeris_day(jd, delta_t);
    let jc = julian_century(jd);
    let jce = julian_ephemeris_century(jde);
    let jme = julian_ephemeris_millennium(jce);
    let l = heliocentric_longitude(jme);
    let b = heliocentric_latitude(jme);
    let r = heliocentric_radius_vector(jme);
    let theta = geocentric_longitude(l);
    let beta = geocentric_latitude(b);
    let x0 = mean_elongation(jce);
    let x1 = mean_anomaly_sun(jce);
    let x2 = mean_anomaly_moon(jce);
    let x3 = moon_argument_latitude(jce);
    let x4 = moon_ascending_longitude(jce);
    let (delta_psi, delta_epsilon) = longitude_obliquity_nutation(jce, x0, x1, x2, x3, x4);
    let epsilon0 = mean_ecliptic_obliquity(jme);
    let epsilon = true_ecliptic_obliquity(epsilon0, delta_epsilon);
    let delta_tau = aberration_correction(r);
    let lambda = apparent_sun_longitude(theta, delta_epsilon, delta_tau);
    let v0 = mean_sidereal_time(jd, jc);
    let v = apparent_sidereal_time(v0, delta_psi, epsilon);
    let alpha = geocentric_sun_right_ascension(lambda, epsilon, beta);
    let delta = geocentric_sun_declination(lambda, epsilon, beta);
    let h = local_hour_angle(v, longitude, alpha);
    let xi = equatorial_horizontal_parallax(r);
    let u = u_term(latitude);
    let x = x_term(u, elevation, latitude);
    let y = y_term(u, elevation, latitude);
    let delta_alpha = parallax_sun_right_ascension(x, xi, h, delta);
    let delta_prime = topocentric_sun_declination(delta, x, y, xi, delta_alpha, h);
    let h_prime = topocentric_local_hour_angle(h, delta_alpha);
    let e0 = topocentric_elevation_angle_without_atmosphere(latitude, delta_prime, h_prime);
    let delta_e = atmospheric_refraction_correction(pressure, temperature, e0, atmos_refract);
    let e = topocentric_elevation_angle(e0, delta_e);
    let theta = topocentric_zenith_angle(e);
    let theta0 = topocentric_zenith_angle(e0);
    let gamma = topocentric_astronomers_azimuth(h_prime, delta_prime, latitude);
    let phi = topocentric_azimuth_angle(gamma);
    SolPosResult {
        theta,
        theta0,
        e,
        e0,
        phi,
    }
}

pub fn julian_day(unixtime: f64) -> f64 {
    // SPA paper is using year, month, day, etc. but using the unixtime directly is faster
    unixtime * 1.0 / SECONDS_PER_DAY + JD_CONSTANT
}

pub fn julian_ephemeris_day(jd: f64, delta_t: f64) -> f64 {
    jd + delta_t * 1. / SECONDS_PER_DAY
}

pub fn julian_century(jd: f64) -> f64 {
    (jd - JC_CONSTANT) / DAYS_PER_CENTURY
}

pub fn julian_ephemeris_century(jde: f64) -> f64 {
    (jde - JC_CONSTANT) / DAYS_PER_CENTURY
}

pub fn julian_ephemeris_millennium(jce: f64) -> f64 {
    jce / 10.
}

macro_rules! sum_table_rows {
    ($table:expr, $jme:expr) => {
        $table
            .iter()
            .map(|&row| row[0] * f64::cos(row[1] + row[2] * $jme))
            .sum::<f64>()
    };
}

pub fn heliocentric_longitude(jme: f64) -> f64 {
    let l0 = sum_table_rows!(tables::L0, jme);
    let l1 = sum_table_rows!(tables::L1, jme);
    let l2 = sum_table_rows!(tables::L2, jme);
    let l3 = sum_table_rows!(tables::L3, jme);
    let l4 = sum_table_rows!(tables::L4, jme);
    let l5: f64 = tables::L5[0] * f64::cos(tables::L5[1] + tables::L5[2] * jme);
    let l_rad = (l0
        + l1 * jme
        + l2 * f64::powf(jme, 2.)
        + l3 * f64::powf(jme, 3.)
        + l4 * f64::powf(jme, 4.)
        + l5 * f64::powf(jme, 5.))
        / f64::powf(10., 8.);

    let l_deg = l_rad * 180. / std::f64::consts::PI;
    l_deg.rem_euclid(360.0)
}

pub fn heliocentric_latitude(jme: f64) -> f64 {
    let b0 = sum_table_rows!(tables::B0, jme);
    let b1 = sum_table_rows!(tables::B1, jme);
    let b_rad = (b0 + b1 * jme) / f64::powf(10., 8.);
    let b_deg = b_rad * 180. / std::f64::consts::PI;
    b_deg
}

pub fn heliocentric_radius_vector(jme: f64) -> f64 {
    let r0 = sum_table_rows!(tables::R0, jme);
    let r1 = sum_table_rows!(tables::R1, jme);
    let r2 = sum_table_rows!(tables::R2, jme);
    let r3 = sum_table_rows!(tables::R3, jme);
    let r4: f64 = tables::R4[0] * f64::cos(tables::R4[1] + tables::R4[2] * jme);
    let r = (r0
        + r1 * jme
        + r2 * f64::powf(jme, 2.)
        + r3 * f64::powf(jme, 3.)
        + r4 * f64::powf(jme, 4.))
        / f64::powf(10., 8.);
    r
}

pub fn geocentric_longitude(l_deg: f64) -> f64 {
    let theta = l_deg + 180.;
    theta.rem_euclid(360.0)
}

pub fn geocentric_latitude(b_deg: f64) -> f64 {
    return -b_deg;
}

pub fn mean_elongation(jce: f64) -> f64 {
    297.85036 + 445267.11148 * jce - 0.0019142 * jce.powi(2) + jce.powi(3) / 189474.
}

pub fn mean_anomaly_sun(jce: f64) -> f64 {
    357.52772 + 35999.05034 * jce - 0.0001603 * jce.powi(2) - jce.powi(3) / 300000.
}

pub fn mean_anomaly_moon(jce: f64) -> f64 {
    134.96298 + 477198.867398 * jce + 0.0086972 * jce.powi(2) + jce.powi(3) / 56250.
}

pub fn moon_argument_latitude(jce: f64) -> f64 {
    93.27191 + 483202.017538 * jce - 0.0036825 * jce.powi(2) + jce.powi(3) / 327270.
}

pub fn moon_ascending_longitude(jce: f64) -> f64 {
    125.04452 - 1934.136261 * jce + 0.0020708 * jce.powi(2) + jce.powi(3) / 450000.
}

pub fn longitude_obliquity_nutation(
    jce: f64,
    x0: f64,
    x1: f64,
    x2: f64,
    x3: f64,
    x4: f64,
) -> (f64, f64) {
    let x = [x0, x1, x2, x3, x4];
    let (mut delta_psi_sum, mut delta_eps_sum) = (0., 0.);

    for (y, abcd) in tables::NUT.iter().zip(tables::ABCD.iter()) {
        let sum_xy = x
            .iter()
            .zip(y.iter())
            .map(|(&xi, &yi)| xi * yi as f64)
            .sum::<f64>()
            * std::f64::consts::PI
            / 180.;

        delta_psi_sum += (abcd[0] + abcd[1] * jce) * f64::sin(sum_xy);
        delta_eps_sum += (abcd[2] + abcd[3] * jce) * f64::cos(sum_xy);
    }

    (delta_psi_sum / 36000000., delta_eps_sum / 36000000.)
}

fn mean_ecliptic_obliquity(jme: f64) -> f64 {
    let u = jme / 10.;
    let e0 = 84381.448 - 4680.93 * u - 1.55 * u.powi(2) + 1999.25 * u.powi(3)
        - 51.38 * u.powi(4)
        - 249.67 * u.powi(5)
        - 39.05 * u.powi(6)
        + 7.12 * u.powi(7)
        + 27.87 * u.powi(8)
        + 5.79 * u.powi(9)
        + 2.45 * u.powi(10);
    e0
}

fn true_ecliptic_obliquity(epsilon0: f64, delta_epsilon: f64) -> f64 {
    epsilon0 / 3600. + delta_epsilon
}

fn aberration_correction(r: f64) -> f64 {
    -20.4898 / (3600. * r)
}

fn apparent_sun_longitude(theta: f64, delta_eps: f64, delta_tau: f64) -> f64 {
    theta + delta_eps + delta_tau
}

fn mean_sidereal_time(jd: f64, jc: f64) -> f64 {
    let v0 = 280.46061837 + 360.98564736629 * (jd - 2451545.) + 0.000387933 * jc.powi(2)
        - jc.powi(3) / 38710000.;
    v0.rem_euclid(360.0)
}

fn apparent_sidereal_time(v0: f64, delta_psi: f64, epsilon: f64) -> f64 {
    v0 + delta_psi * f64::cos(epsilon * DEG_TO_RAD)
}

fn geocentric_sun_right_ascension(lambda: f64, epsilon: f64, beta: f64) -> f64 {
    let lambda_rad = lambda * DEG_TO_RAD;
    let epsilon_rad = epsilon * DEG_TO_RAD;
    let beta_rad = beta * DEG_TO_RAD;
    let num =
        f64::sin(lambda_rad) * f64::cos(epsilon_rad) - f64::tan(beta_rad) * f64::sin(epsilon_rad);
    let alpha = f64::atan2(num, f64::cos(lambda_rad)) * RAD_TO_DEG;
    alpha.rem_euclid(360.0)
}

fn geocentric_sun_declination(lambda: f64, epsilon: f64, beta: f64) -> f64 {
    let lambda_rad = lambda * std::f64::consts::PI / 180.;
    let epsilon_rad = epsilon * std::f64::consts::PI / 180.;
    let beta_rad = beta * std::f64::consts::PI / 180.;
    let num = f64::sin(beta_rad) * f64::cos(epsilon_rad)
        + f64::cos(beta_rad) * f64::sin(epsilon_rad) * f64::sin(lambda_rad);
    let mut delta = f64::asin(num);
    delta = delta * 180. / std::f64::consts::PI;
    delta
}

fn local_hour_angle(v: f64, lon: f64, alpha: f64) -> f64 {
    let h = v + lon - alpha;
    h.rem_euclid(360.0)
}

fn equatorial_horizontal_parallax(r: f64) -> f64 {
    8.794 / (3600. * r)
}

fn u_term(lat: f64) -> f64 {
    let lat_rad = lat * DEG_TO_RAD;
    let u = f64::atan(0.99664719 * f64::tan(lat_rad));
    u
}

fn x_term(u: f64, elev: f64, lat: f64) -> f64 {
    let lat_rad = lat * DEG_TO_RAD;
    let x = f64::cos(u) + (elev / 6378140.) * f64::cos(lat_rad);
    x
}

fn y_term(u: f64, elev: f64, lat: f64) -> f64 {
    let lat_rad = lat * DEG_TO_RAD;
    let y = 0.99664719 * f64::sin(u) + (elev / 6378140.) * f64::sin(lat_rad);
    y
}

fn parallax_sun_right_ascension(x: f64, xi: f64, h: f64, delta: f64) -> f64 {
    let xi_rad = xi * DEG_TO_RAD;
    let h_rad = h * DEG_TO_RAD;
    let delta_rad = delta * std::f64::consts::PI / 180.;
    let num = -x * f64::sin(xi_rad) * f64::sin(h_rad);
    let denom = f64::cos(delta_rad) - x * f64::sin(xi_rad) * f64::cos(h_rad);
    f64::atan2(num, denom) * RAD_TO_DEG
}

fn topocentric_sun_declination(
    delta: f64,
    x: f64,
    y: f64,
    xi: f64,
    delta_alpha: f64,
    h: f64,
) -> f64 {
    let xi_rad = xi * DEG_TO_RAD;
    let h_rad = h * DEG_TO_RAD;
    let delta_rad = delta * std::f64::consts::PI / 180.;
    let delta_alpha_rad = delta_alpha * std::f64::consts::PI / 180.;
    let num = (f64::sin(delta_rad) - y * f64::sin(xi_rad)) * f64::cos(delta_alpha_rad);
    let denom = f64::cos(delta_rad) - x * f64::sin(xi_rad) * f64::cos(h_rad);
    let delta_prime = f64::atan2(num, denom);
    delta_prime * RAD_TO_DEG
}

fn topocentric_local_hour_angle(h: f64, delta_alpha: f64) -> f64 {
    h - delta_alpha
}

fn topocentric_elevation_angle_without_atmosphere(lat: f64, delta_prime: f64, h_prime: f64) -> f64 {
    let lat_rad = lat * DEG_TO_RAD;
    let delta_prime_rad = delta_prime * DEG_TO_RAD;
    let h_prime_rad = h_prime * DEG_TO_RAD;
    let e0 = f64::asin(
        f64::sin(lat_rad) * f64::sin(delta_prime_rad)
            + f64::cos(lat_rad) * f64::cos(delta_prime_rad) * f64::cos(h_prime_rad),
    );
    e0 * RAD_TO_DEG
}

fn atmospheric_refraction_correction(pressure: f64, temp: f64, e0: f64, atmos_refract: f64) -> f64 {
    if e0 >= (-0.26667 - atmos_refract) {
        let term = (e0 + 10.3 / (e0 + 5.11)) * DEG_TO_RAD;
        let delta_e = pressure / 1010. * 283. / (273. + temp) * 1.02 / (60. * f64::tan(term));
        return delta_e;
    }
    0.
}

fn topocentric_elevation_angle(e0: f64, delta_e: f64) -> f64 {
    e0 + delta_e
}

fn topocentric_zenith_angle(e: f64) -> f64 {
    90. - e
}

fn topocentric_astronomers_azimuth(h_prime: f64, delta_prime: f64, lat: f64) -> f64 {
    let h_prime_rad = h_prime * DEG_TO_RAD;
    let delta_prime_rad = delta_prime * DEG_TO_RAD;
    let lat_rad = lat * DEG_TO_RAD;
    let den =
        f64::cos(h_prime_rad) * f64::sin(lat_rad) - f64::tan(delta_prime_rad) * f64::cos(lat_rad);
    let mut gamma = f64::atan2(f64::sin(h_prime_rad), den);
    gamma = gamma * RAD_TO_DEG;
    gamma.rem_euclid(360.0)
}

fn topocentric_azimuth_angle(gamma: f64) -> f64 {
    let phi = gamma + 180.;
    ((phi % 360.) + 360.) % 360.
}

pub fn calculate_deltat(year: f64, month: f64) -> f64 {
    // https://eclipse.gsfc.nasa.gov/SEhelp/deltaT.html
    // https://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html
    // FIXME: Below should raise errors
    if (year < -1999.) | (year > 3000.) {
        panic!("`year` should be above -1999 and below 3000");
    }
    if (month < 1.) | (month > 12.) {
        panic!("`month` should be above 0 and below 13");
    }
    let y = year + (month - 0.5) / 12.;
    // TODO: Replace with match statement
    if year < -500. {
        return -20. + 32. * f64::powf((y - 1820.) / 100., 2.);
    } else if (year >= -500.) & (year < 500.) {
        return 10583.6 - 1014.41 * (y / 100.) + 33.78311 * f64::powf(y / 100., 2.)
            - 5.952053 * f64::powf(y / 100., 3.)
            - 0.1798452 * f64::powf(y / 100., 4.)
            + 0.022174192 * f64::powf(y / 100., 5.)
            + 0.0090316521 * f64::powf(y / 100., 6.);
    } else if (year >= 500.) & (year < 1600.) {
        return 1574.2 - 556.01 * (y - 1000.) / 100.
            + 71.23472 * f64::powf((y - 1000.) / 100., 2.)
            + 0.319781 * f64::powf((y - 1000.) / 100., 3.)
            - 0.8503463 * f64::powf((y - 1000.) / 100., 4.)
            - 0.005050998 * f64::powf((y - 1000.) / 100., 5.)
            + 0.0083572073 * f64::powf((y - 1000.) / 100., 6.);
    } else if (year >= 1600.) & (year < 1700.) {
        return 120. - 0.9808 * (y - 1600.) - 0.01532 * f64::powf(y - 1600., 2.)
            + f64::powf(y - 1600., 3.) / 7129.;
    } else if (year >= 1700.) & (year < 1800.) {
        return 8.83 + 0.1603 * (y - 1700.) - 0.0059285 * f64::powf(y - 1700., 2.)
            + 0.00013336 * f64::powf(y - 1700., 3.)
            - f64::powf(y - 1700., 4.) / 1174000.;
    } else if (year >= 1800.) & (year < 1860.) {
        return 13.72 - 0.332447 * (y - 1800.)
            + 0.0068612 * f64::powf(y - 1800., 2.)
            + 0.0041116 * f64::powf(y - 1800., 3.)
            - 0.00037436 * f64::powf(y - 1800., 4.)
            + 0.0000121272 * f64::powf(y - 1800., 5.)
            - 0.0000001699 * f64::powf(y - 1800., 6.)
            + 0.000000000875 * f64::powf(y - 1800., 7.);
    } else if (year >= 1860.) & (year < 1900.) {
        return 7.62 + 0.5737 * (y - 1860.) - 0.251754 * f64::powf(y - 1860., 2.)
            + 0.01680668 * f64::powf(y - 1860., 3.)
            - 0.0004473624 * f64::powf(y - 1860., 4.)
            + f64::powf(y - 1860., 5.) / 233174.;
    } else if (year >= 1900.) & (year < 1920.) {
        return -2.79 + 1.494119 * (y - 1900.) - 0.0598939 * f64::powf(y - 1900., 2.)
            + 0.0061966 * f64::powf(y - 1900., 3.)
            - 0.000197 * f64::powf(y - 1900., 4.);
    } else if (year >= 1920.) & (year < 1941.) {
        return 21.2 + 0.84493 * (y - 1920.) - 0.0761 * f64::powf(y - 1920., 2.)
            + 0.0020936 * f64::powf(y - 1920., 3.);
    } else if (year >= 1941.) & (year < 1961.) {
        return 29.07 + 0.407 * (y - 1950.) - f64::powf(y - 1950., 2.) / 233.
            + f64::powf(y - 1950., 3.) / 2547.;
    } else if (year >= 1961.) & (year < 1986.) {
        return 45.45 + 1.067 * (y - 1975.)
            - f64::powf(y - 1975., 2.) / 260.
            - f64::powf(y - 1975., 3.) / 718.;
    } else if (year >= 1986.) & (year < 2005.) {
        return 63.86 + 0.3345 * (y - 2000.) - 0.060374 * f64::powf(y - 2000., 2.)
            + 0.0017275 * f64::powf(y - 2000., 3.)
            + 0.000651814 * f64::powf(y - 2000., 4.)
            + 0.00002373599 * f64::powf(y - 2000., 5.);
    } else if (year >= 2005.) & (year < 2050.) {
        return 62.92 + 0.32217 * (y - 2000.) + 0.005589 * f64::powf(y - 2000., 2.);
    } else if (year >= 2050.) & (year < 2150.) {
        return -20. + 32. * f64::powf((y - 1820.) / 100., 2.) - 0.5628 * (2150. - y);
    }
    return -20. + 32. * f64::powf((y - 1820.) / 100., 2.);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_deltat() {
        let expected = [
            17217.912272222224,
            17211.72360555555,
            17184.831552348543,
            17177.310907656512,
            118.96186861631928,
            118.5375201220693,
            8.99069606382023,
            9.051572240380537,
            13.385367077090038,
            13.260910643977612,
            7.962908405417306,
            7.9713621281635225,
            -1.2918431282701737,
            -0.720127316638437,
            21.99992798755792,
            22.27683824508097,
            24.797268089496473,
            25.035728072816564,
            33.594798592945274,
            33.75427944537937,
            54.89627599023825,
            55.079355914075165,
            64.68633720312503,
            64.84502657812497,
            93.0847888888886,
            93.93328888888917,
            328.56800555555526,
            329.4486722222225,
        ];
        let years = [
            -501., -499., 1601., 1701., 1801., 1861., 1901., 1921., 1941., 1961., 1986., 2005.,
            2050., 2150.,
        ];
        for (i, year) in years.iter().enumerate() {
            for (j, month) in [1 as u8, 6 as u8].iter().enumerate() {
                assert_eq!(expected[i * 2 + j], calculate_deltat(*year, *month as f64));
            }
        }
    }

    #[test]
    #[should_panic]
    fn test_calculate_deltat_wrong_year_3001() {
        let year = 3001.;
        let month = 1.;
        calculate_deltat(year, month);
    }

    #[test]
    #[should_panic]
    fn test_calculate_deltat_wrong_year_min_2000() {
        let year = -2000.;
        let month = 1.;
        calculate_deltat(year, month);
    }

    #[test]
    #[should_panic]
    fn test_calculate_deltat_wrong_month_0() {
        let year = 2000.;
        let month = 0.;
        calculate_deltat(year, month);
    }

    #[test]
    #[should_panic]
    fn test_calculate_deltat_wrong_month_13() {
        let year = 2000.;
        let month = 13.;
        calculate_deltat(year, month);
    }

    #[test]
    fn test_julian_day() {
        let timestamp = 1357000200.;
        assert_eq!(2456293.5208333335, julian_day(timestamp));
        // taken from A4.1:
        let timestamp = 582724800.;
        assert_eq!(2447332.0, julian_day(timestamp));
    }

    #[test]
    fn test_julian_ephemeris_day() {
        let jd = 2456293.5208333335;
        let delta_t = 64.5;
        assert_eq!(2456293.521579861, julian_ephemeris_day(jd, delta_t));
    }

    #[test]
    fn test_julian_century() {
        let jd = 2456293.5208333335;
        assert_eq!(0.13000741501255272, julian_century(jd));
    }

    #[test]
    fn test_julian_ephemeris_century() {
        let jde = 2456293.521579861;
        assert_eq!(0.1300074354513669, julian_ephemeris_century(jde));
    }

    #[test]
    fn test_julian_ephemeris_millennium() {
        let jde = 0.1300074354513669;
        assert_eq!(0.01300074354513669, julian_ephemeris_millennium(jde));
    }

    #[test]
    fn test_heliocentric_longitude() {
        let jme = 0.01300074354513669;
        // pvlib:  100.75286502047948
        assert_eq!(100.75286504017458, heliocentric_longitude(jme));
        let jme = -0.01300074354513669;
        // pvlib   100.00464976366129
        assert_eq!(100.0046497683752, heliocentric_longitude(jme));
    }

    #[test]
    fn test_heliocentric_latitude() {
        let jme = 0.01300074354513669;
        // pvlib   0.00010663876509755685
        assert_eq!(0.00010663876509755686, heliocentric_latitude(jme));
        let jme = -0.01300074354513669;
        // pvlib   0.00010214440296148753
        assert_eq!(0.00010214440296148753, heliocentric_latitude(jme));
    }

    #[test]
    fn test_heliocentric_radius_vector() {
        let jme = 0.01300074354513669;
        // pvlib   0.9832940565275712
        assert_eq!(0.9832940565275712, heliocentric_radius_vector(jme));
        let jme = -0.01300074354513669;
        // pvlib   0.9833347325330166
        assert_eq!(0.9833347325330166, heliocentric_radius_vector(jme));
    }

    #[test]
    fn test_geocentric_longitude() {
        let l_deg = 90.;
        assert_eq!(270., geocentric_longitude(l_deg));
        let l_deg = 350.;
        assert_eq!(170., geocentric_longitude(l_deg));
    }

    #[test]
    fn test_geocentric_latitude() {
        let b_deg = 90.;
        assert_eq!(-90., geocentric_latitude(b_deg));
    }

    #[test]
    fn test_mean_elongation() {
        let jce = 0.13000741501255272;
        assert_eq!(58185.876481278865, mean_elongation(jce));
    }

    #[test]
    fn test_mean_anomaly_sun() {
        let jce = 0.13000741501255272;
        assert_eq!(5037.671194893454, mean_anomaly_sun(jce));
    }

    #[test]
    fn test_mean_anomaly_moon() {
        let jce = 0.13000741501255272;
        assert_eq!(62174.354324370404, mean_anomaly_moon(jce));
    }

    #[test]
    fn test_moon_argument_latitude() {
        let jce = 0.13000741501255272;
        assert_eq!(62913.117076730916, moon_argument_latitude(jce));
    }

    #[test]
    fn test_moon_ascending_longitude() {
        let jce = 0.13000741501255272;
        assert_eq!(-126.40750056925852, moon_ascending_longitude(jce));
    }

    #[test]
    fn test_longitude_obliquity_nutation() {
        let jce = 0.13000741501255272;
        let x0 = 58185.876481278865;
        let x1 = 5037.671194893454;
        let x2 = 62174.354324370404;
        let x3 = 62913.117076730916;
        let x4 = -126.40750056925852;
        let (n_longitude, n_obliquity) = longitude_obliquity_nutation(jce, x0, x1, x2, x3, x4);
        // pvlib 0.004069674875572919
        assert_eq!(0.0040696748755729446, n_longitude);
        // pvlib -0.0016509592741862002
        assert_eq!(-0.0016509592741862074, n_obliquity);
    }

    #[test]
    fn test_mean_ecliptic_obliquity() {
        let jme = 0.01300074354513669;
        assert_eq!(84375.3624447249, mean_ecliptic_obliquity(jme));
    }

    #[test]
    fn test_true_ecliptic_obliquity() {
        let epsilon0 = 84375.3624447249;
        let delta_epsilon = -0.0016509592741862074;
        assert_eq!(
            23.435949719816065,
            true_ecliptic_obliquity(epsilon0, delta_epsilon)
        );
    }

    #[test]
    fn test_aberration_correction() {
        let r = 0.983294088577116;
        assert_eq!(-0.005788310107047633, aberration_correction(r));
    }

    #[test]
    fn test_apparent_sun_longitude() {
        let theta = 1.;
        let delta_epsilon = 1.;
        let delta_tau = 1.;
        assert_eq!(3., apparent_sun_longitude(theta, delta_epsilon, delta_tau));
    }

    #[test]
    fn test_mean_sidereal_time() {
        let jd = 2456293.5208333335;
        let jc = 0.13000741501255272;
        // pvlib   108.3276781309396
        assert_eq!(108.3276781309396, mean_sidereal_time(jd, jc));
    }

    #[test]
    fn test_apparent_sidereal_time() {
        let v0 = 108.32767813109967;
        let delta_psi = 0.0040696748755729446;
        let epsilon = 23.435949719816065;
        assert_eq!(
            108.33141207919714,
            apparent_sidereal_time(v0, delta_psi, epsilon)
        );
    }

    #[test]
    fn test_geocentric_sun_right_ascension() {
        let lambda = 170.;
        let epsilon = 23.435949719816065;
        let beta = -0.0016509592741862074;
        assert_eq!(
            170.80960819870663,
            geocentric_sun_right_ascension(lambda, epsilon, beta)
        );
    }

    #[test]
    fn test_geocentric_sun_declination() {
        let lambda = 170.;
        let epsilon = 23.435949719816065;
        let beta = -0.0016509592741862074;
        assert_eq!(
            3.9587091134966355,
            geocentric_sun_declination(lambda, epsilon, beta)
        );
    }

    #[test]
    fn test_local_hour_angle() {
        let v = 170.;
        let lon = 100.;
        let alpha = 20.;
        assert_eq!(250., local_hour_angle(v, lon, alpha));
        let v = 170.;
        let lon = 250.;
        let alpha = 20.;
        assert_eq!(40., local_hour_angle(v, lon, alpha));
    }

    #[test]
    fn test_equatorial_horizontal_parallax() {
        let r = 0.983294088577116;
        assert_eq!(0.002484279938377968, equatorial_horizontal_parallax(r));
    }

    #[test]
    fn test_u_term() {
        let lat = 100.;
        assert_eq!(-1.3956881668836334, u_term(lat));
    }

    #[test]
    fn test_x_term() {
        let u = -1.3956881668836334;
        let elev = 200.;
        let lat = 100.;
        assert_eq!(0.17420919940606158, x_term(u, elev, lat));
    }

    #[test]
    fn test_y_term() {
        let u = -1.3956881668836334;
        let elev = 200.;
        let lat = 100.;
        assert_eq!(-0.9813752830759921, y_term(u, elev, lat));
    }

    #[test]
    fn test_parallax_sun_right_ascension() {
        let x = 0.17420919940606158;
        let xi = 0.002484279938377968;
        let h = 40.;
        let delta = 3.9587091134966355;
        // pvlib -0.0002788554074644773
        assert_eq!(
            -0.0002788554074644773,
            parallax_sun_right_ascension(x, xi, h, delta)
        )
    }

    #[test]
    fn test_topocentric_sun_declination() {
        let delta = 3.9587091134966355;
        let x = 0.17420919940606158;
        let y = -0.9813752830759921;
        let xi = 0.002484279938377968;
        let delta_alpha = 50.;
        let h = 40.;
        assert_eq!(
            2.548569550651887,
            topocentric_sun_declination(delta, x, y, xi, delta_alpha, h)
        );
    }

    #[test]
    fn test_topocentric_local_hour_angle() {
        let h = 50.;
        let delta_alpha = 30.;
        assert_eq!(20., topocentric_local_hour_angle(h, delta_alpha));
    }

    #[test]
    fn test_topocentric_elevation_angle_without_atmosphere() {
        let lat = 100.;
        let delta_prime = 2.548569550651887;
        let h_prime = 20.;
        assert_eq!(
            -6.847307469841585,
            topocentric_elevation_angle_without_atmosphere(lat, delta_prime, h_prime)
        );
    }

    #[test]
    fn test_atmospheric_refraction_correction() {
        let pressure = 1000.;
        let temp = 20.;
        let e0 = -6.847307469841585;
        let atmos_refract = 2.;
        assert_eq!(
            0.,
            atmospheric_refraction_correction(pressure, temp, e0, atmos_refract)
        );
        let e0: f64 = 6.847307469841585;
        assert_eq!(
            0.12010357702796845,
            atmospheric_refraction_correction(pressure, temp, e0, atmos_refract)
        );
    }

    #[test]
    fn test_topocentric_elevation_angle() {
        let e0 = 1.;
        let delta_e = 1.;
        assert_eq!(2., topocentric_elevation_angle(e0, delta_e));
    }

    #[test]
    fn test_topocentric_zenith_angle() {
        let e = 45.;
        assert_eq!(45., topocentric_zenith_angle(e));
    }

    #[test]
    fn test_topocentric_astronomers_azimuth() {
        let h_prime = 20.;
        let delta_prime = 2.548569550651887;
        let lat = 100.;
        assert_eq!(
            20.129089268987894,
            topocentric_astronomers_azimuth(h_prime, delta_prime, lat)
        );
    }

    #[test]
    fn test_topocentric_azimuth_angle() {
        let gamma = 20.;
        assert_eq!(200., topocentric_azimuth_angle(gamma));
        let gamma = 220.;
        assert_eq!(40., topocentric_azimuth_angle(gamma));
        let gamma = -200.;
        assert_eq!(340., topocentric_azimuth_angle(gamma));
    }
}
