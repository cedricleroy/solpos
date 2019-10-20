use std;
use dict_derive::{IntoPyObject};

mod tables;


#[derive(IntoPyObject)]
pub struct Result {
    pub theta: f64,
    pub theta0: f64,
    pub e: f64,
    pub e0: f64,
    pub phi: f64
}


impl std::fmt::Display for Result {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "theta: {}, theta0: {} e: {}, e0: {}, phi: {}",
               self.theta, self.theta0, self.e, self.e0, self.phi)
    }
}


pub fn solar_position(unixtime: &u64, lat: &f32, lon: &f32, elev: &f32,
                      pressure: &f32, temp: &f32, delta_t: &f64,
                      atmos_refract: &f32) -> Result {
    let jd = julian_day(&unixtime);
    let jde = julian_ephemeris_day(&jd, &delta_t); 
    let jc = julian_century(&jd);
    let jce = julian_ephemeris_century(&jde);
    let jme = julian_ephemeris_millennium(&jce);
    let l = heliocentric_longitude(&jme);
    let b = heliocentric_latitude(&jme);
    let r = heliocentric_radius_vector(&jme);
    let theta = geocentric_longitude(&l);
    let beta = geocentric_latitude(&b);
    let x0 = mean_elongation(&jce);
    let x1 = mean_anomaly_sun(&jce);
    let x2 = mean_anomaly_moon(&jce);
    let x3 = moon_argument_latitude(&jce);
    let x4 = moon_ascending_longitude(&jce);
    let (delta_psi, delta_epsilon) = longitude_obliquity_nutation(&jce, &x0, &x1, &x2, &x3, &x4);
    let epsilon0 = mean_ecliptic_obliquity(&jme);
    let epsilon = true_ecliptic_obliquity(&epsilon0, &delta_epsilon);
    let delta_tau = aberration_correction(&r);
    let lambda = apparent_sun_longitude(&theta, &delta_epsilon, &delta_tau);
    let v0 = mean_sidereal_time(&jd, &jc);
    let v = apparent_sidereal_time(&v0, &delta_psi, &epsilon);
    let alpha = geocentric_sun_right_ascension(&lambda, &epsilon, &beta);
    let delta = geocentric_sun_declination(&lambda, &epsilon, &beta);
    let h = local_hour_angle(&v, &lon, &alpha);
    let xi = equatorial_horizontal_parallax(&r);
    let u = u_term(&lat);
    let x = x_term(&u, &elev, &lat);
    let y = y_term(&u, &elev, &lat);
    let delta_alpha = parallax_sun_right_ascension(&x, &xi, &h, &delta);
    //let alpha_prime = topocentric_sun_right_ascension(&alpha, &delta_alpha);
    let delta_prime = topocentric_sun_declination(&delta, &x, &y, &xi, &delta_alpha, &h);
    let h_prime = topocentric_local_hour_angle(&h, &delta_alpha);
    let e0 = topocentric_elevation_angle_without_atmosphere(&lat, &delta_prime, &h_prime);
    let delta_e = atmospheric_refraction_correction(&pressure, &temp, &e0, &atmos_refract);
    let e = topocentric_elevation_angle(&e0, &delta_e);
    let theta = topocentric_zenith_angle(&e);
    let theta0 = topocentric_zenith_angle(&e0);
    let gamma = topocentric_astronomers_azimuth(&h_prime, &delta_prime, &lat);
    let phi = topocentric_azimuth_angle(&gamma);
    Result{theta, theta0, e, e0, phi}
}


fn julian_day(unixtime: &u64) -> f64 {
    // differ from SPA papper
    let jd = *unixtime as f64 * 1.0 / 86400. + 2440587.5;
    jd
}


fn julian_ephemeris_day(jd: &f64, delta_t: &f64) -> f64 {
    let jde = jd + delta_t * 1. / 86400.;
    jde
}


fn julian_century(jd: &f64) -> f64 {
    let jc = (jd - 2451545.) / 36525.;
    jc
}


fn julian_ephemeris_century(jde: &f64) -> f64 {
    let jce = (jde - 2451545.) / 36525.;
    jce
}


fn julian_ephemeris_millennium(jce: &f64) -> f64 {
    let jme = jce / 10.;
    jme
}


fn limit_angle_range(an: &mut f64) {
    let l_360 = *an / 360.;
    let f = l_360 - (l_360 as i64) as f64;
    if *an >= 0. {
        *an = 360. * f
    } else {
        *an = 360. - 360. * f.abs()
    }
}


fn heliocentric_longitude(jme: &f64) -> f64 {
    // FIXME: too verbose and not DRY
    let mut l0i: [f64; 65] = [0.; 65];
    for (i, row) in tables::L0.iter().enumerate() {
        l0i[i] = row[0] as f64 * f64::cos(row[1] as f64 + row[2] as f64 * jme);
    }
    let l0: f64 = l0i.iter().sum();
    let mut l1i: [f64; 34] = [0.; 34];
    for (i, row) in tables::L1.iter().enumerate() {
        l1i[i] = row[0] as f64 * f64::cos(row[1] as f64 + row[2] as f64 * jme);
    }
    let l1: f64 = l1i.iter().sum();
    let mut l2i: [f64; 20] = [0.; 20];
    for (i, row) in tables::L2.iter().enumerate() {
        l2i[i] = row[0] as f64 * f64::cos(row[1] as f64 + row[2] as f64 * jme);
    }
    let l2: f64 = l2i.iter().sum();
    let mut l3i: [f64; 7] = [0.; 7];
    for (i, row) in tables::L3.iter().enumerate() {
        l3i[i] = row[0] as f64 * f64::cos(row[1] as f64 + row[2] as f64 * jme);
    }
    let l3: f64 = l3i.iter().sum();
    let mut l4i: [f64; 3] = [0.; 3];
    for (i, row) in tables::L4.iter().enumerate() {
        l4i[i] = row[0] as f64 * f64::cos(row[1] as f64 + row[2] as f64 * jme);
    }
    let l4: f64 = l4i.iter().sum();
    let l5: f64 = tables::L5[0] as f64 * f64::cos(tables::L5[1] as f64 + tables::L5[2] as f64 * jme);
    let l_rad = (l0 + l1 * jme
                 + l2 * f64::powf(*jme, 2.) 
                 + l3 * f64::powf(*jme, 3.) 
                 + l4 * f64::powf(*jme, 4.) 
                 + l5 * f64::powf(*jme, 5.)) / f64::powf(10., 8.);
    let mut l_deg = l_rad * 180. / std::f64::consts::PI;
    limit_angle_range(&mut l_deg);
    l_deg
}


fn heliocentric_latitude(jme: &f64) -> f64 {
    let mut b0i: [f64; 5] = [0.; 5];
    for (i, row) in tables::B0.iter().enumerate() {
        b0i[i] = row[0] as f64 * f64::cos(row[1] as f64 + row[2] as f64 * jme);
    }
    let b0: f64 = b0i.iter().sum();
    let mut b1i: [f64; 2] = [0.; 2];
    for (i, row) in tables::B1.iter().enumerate() {
        b1i[i] = row[0] as f64 * f64::cos(row[1] as f64 + row[2] as f64 * jme);
    }
    let b1: f64 = b1i.iter().sum();
    let b_rad = (b0 + b1 * jme) / f64::powf(10., 8.);
    let b_deg = b_rad * 180. / std::f64::consts::PI;
    b_deg
}


fn heliocentric_radius_vector(jme: &f64) -> f64 {
    let mut r0i: [f64; 40] = [0.; 40];
    for (i, row) in tables::R0.iter().enumerate() {
        r0i[i] = row[0] as f64 * f64::cos(row[1] as f64 + row[2] as f64 * jme);
    }
    let r0: f64 = r0i.iter().sum();
    let mut r1i: [f64; 10] = [0.; 10];
    for (i, row) in tables::R1.iter().enumerate() {
        r1i[i] = row[0] as f64 * f64::cos(row[1] as f64 + row[2] as f64 * jme);
    }
    let r1: f64 = r1i.iter().sum();
    let mut r2i: [f64; 6] = [0.; 6];
    for (i, row) in tables::R2.iter().enumerate() {
        r2i[i] = row[0] as f64 * f64::cos(row[1] as f64 + row[2] as f64 * jme);
    }
    let r2: f64 = r2i.iter().sum();
    let mut r3i: [f64; 2] = [0.; 2];
    for (i, row) in tables::R3.iter().enumerate() {
        r3i[i] = row[0] as f64 * f64::cos(row[1] as f64 + row[2] as f64 * jme);
    }
    let r3: f64 = r3i.iter().sum();
    let r4: f64 = tables::R4[0] as f64 * f64::cos(tables::R4[1] as f64 + tables::R4[2] as f64 * jme);
    let r = (r0 + r1 * jme
             + r2 * f64::powf(*jme, 2.) 
             + r3 * f64::powf(*jme, 3.) 
             + r4 * f64::powf(*jme, 4.)) / f64::powf(10., 8.); 
    r
}


fn geocentric_longitude(l_deg: &f64) -> f64 {
    let mut theta = l_deg + 180.;
    limit_angle_range(&mut theta);
    theta
}


fn geocentric_latitude(b_deg: &f64) -> f64 {
    return -b_deg
}


fn mean_elongation(jce: &f64) -> f64 {
    let x0 = 297.85036 + 445267.11148 * jce
             - 0.0019142 * f64::powf(*jce, 2.)
             + f64::powf(*jce, 3.) / 189474.;
    x0
}


fn mean_anomaly_sun(jce: &f64) -> f64 {
    let x1 = 357.52772 + 35999.05034 * jce
             - 0.0001603 * f64::powf(*jce, 2.)
             - f64::powf(*jce, 3.) / 300000.;
    x1
}


fn mean_anomaly_moon(jce: &f64) -> f64 {
    let x2 = 134.96298 + 477198.867398 * jce
             + 0.0086972 * f64::powf(*jce, 2.)
             + f64::powf(*jce, 3.) / 56250.;
    x2
}


fn moon_argument_latitude(jce: &f64) -> f64 {
    let x3 = 93.27191 + 483202.017538 * jce
             - 0.0036825 * f64::powf(*jce, 2.)
             + f64::powf(*jce, 3.) / 327270.;
    x3
}


fn moon_ascending_longitude(jce: &f64) -> f64 {
    let x4 = 125.04452 - 1934.136261 * jce
             + 0.0020708 * f64::powf(*jce, 2.)
             + f64::powf(*jce, 3.) / 450000.;
    x4
}


fn longitude_obliquity_nutation(jce: &f64, x0: &f64, x1: &f64, x2: &f64,
                                x3: &f64, x4: &f64) -> (f64, f64) {
    let mut delta_psi: [f64; 63] = [0.; 63];
    let mut delta_eps: [f64; 63] = [0.; 63];
    let x = [x0, x1, x2, x3, x4];
    let mut sum_xy;
    for (i, y) in tables::NUT.iter().enumerate() {
        sum_xy = 0.;
        for j in 0..5 {
            sum_xy += x[j] * y[j] as f64;
        }
        sum_xy = sum_xy * std::f64::consts::PI / 180.;
        delta_psi[i] = (tables::ABCD[i][0] + tables::ABCD[i][1] * jce) * f64::sin(sum_xy);
        delta_eps[i] = (tables::ABCD[i][2] + tables::ABCD[i][3] * jce) * f64::cos(sum_xy);
    }
    let delta_psi_sum: f64 = delta_psi.iter().sum();
    let delta_eps_sum: f64 = delta_eps.iter().sum();
    return (delta_psi_sum / 36000000., delta_eps_sum / 36000000.)
}


fn mean_ecliptic_obliquity(jme: &f64) -> f64 {
    let u = jme / 10.;
    let e0 = 84381.448 - 4680.93 * u
             -1.55 * f64::powf(u, 2.)
             + 1999.25 * f64::powf(u, 3.)
             - 51.38 * f64::powf(u, 4.)
             - 249.67 * f64::powf(u, 5.)
             - 39.05 * f64::powf(u, 6.)
             + 7.12 * f64::powf(u, 7.)
             + 27.87 * f64::powf(u, 8.)
             + 5.79 * f64::powf(u, 9.)
             + 2.45 * f64::powf(u, 10.);
    e0
}


fn true_ecliptic_obliquity(epsilon0: &f64, delta_epsilon: &f64) -> f64 {
    let epsilon = epsilon0 / 3600. + delta_epsilon;
    epsilon
}


fn aberration_correction(r: &f64) -> f64 {
    let delta_tau = -20.4898 / (3600. * r);
    delta_tau
}


fn apparent_sun_longitude(theta: &f64, delta_eps: &f64, delta_tau: &f64) -> f64 {
    let lambda = theta + delta_eps + delta_tau;
    lambda
}


fn mean_sidereal_time(jd: &f64, jc: &f64) -> f64 {
    let mut v0 = 280.46061837 + 360.98564736629 * (jd - 2451545.)
                 + 0.000387933 * f64::powf(*jc, 2.)
                 - f64::powf(*jc, 3.) / 38710000.;
    limit_angle_range(&mut v0);
    v0
}


fn apparent_sidereal_time(v0: &f64, delta_psi: &f64, epsilon: &f64) -> f64 {
    let v = v0 + delta_psi * f64::cos(epsilon * std::f64::consts::PI / 180.);
    v
}


fn geocentric_sun_right_ascension(lambda: &f64, epsilon: &f64, beta: &f64) -> f64 {
    let lambda_rad = lambda * std::f64::consts::PI / 180.;
    let epsilon_rad = epsilon * std::f64::consts::PI / 180.;
    let beta_rad = beta * std::f64::consts::PI / 180.;
    let num = f64::sin(lambda_rad) * f64::cos(epsilon_rad)
              - f64::tan(beta_rad) * f64::sin(epsilon_rad);
    let mut alpha = f64::atan2(num, f64::cos(lambda_rad));
    alpha = alpha * 180. / std::f64::consts::PI;
    limit_angle_range(&mut alpha);
    alpha
}


fn geocentric_sun_declination(lambda: &f64, epsilon: &f64, beta: &f64) -> f64 {
    let lambda_rad = lambda * std::f64::consts::PI / 180.;
    let epsilon_rad = epsilon * std::f64::consts::PI / 180.;
    let beta_rad = beta * std::f64::consts::PI / 180.;
    let num = f64::sin(beta_rad) * f64::cos(epsilon_rad)
              + f64::cos(beta_rad) * f64::sin(epsilon_rad) * f64::sin(lambda_rad);
    let mut delta = f64::asin(num);
    delta = delta * 180. / std::f64::consts::PI;
    delta
}


fn local_hour_angle(v: &f64, lon: &f32, alpha: &f64) -> f64 {
    let mut h = v + *lon as f64 - alpha;
    limit_angle_range(&mut h);
    h
}


fn equatorial_horizontal_parallax(r: &f64) -> f64 {
    let xi = 8.794 / (3600. * r);
    xi
}


fn u_term(lat: &f32) -> f64 {
    let lat_rad = *lat as f64 * std::f64::consts::PI / 180.;
    let u = f64::atan(0.99664719 * f64::tan(lat_rad));
    u
}


fn x_term(u: &f64, elev: &f32, lat: &f32) -> f64 {
    let lat_rad = *lat as f64 * std::f64::consts::PI / 180.;
    let x = f64::cos(*u) + (*elev as f64 / 6378140.) * f64::cos(lat_rad);
    x
}


fn y_term(u: &f64, elev: &f32, lat: &f32) -> f64 {
    let lat_rad = *lat as f64 * std::f64::consts::PI / 180.;
    let y = 0.99664719 * f64::sin(*u) + (*elev as f64 / 6378140.) * f64::sin(lat_rad);
    y
}


fn parallax_sun_right_ascension(x: &f64, xi: &f64, h: &f64, delta: &f64) -> f64 {
    let xi_rad = xi * std::f64::consts::PI / 180.;
    let h_rad = h * std::f64::consts::PI / 180.;
    let delta_rad = delta * std::f64::consts::PI / 180.;
    let num = -x * f64::sin(xi_rad) * f64::sin(h_rad);
    let denom = f64::cos(delta_rad) - x * f64::sin(xi_rad) * f64::cos(h_rad);
    let mut delta_alpha = f64::atan2(num, denom);
    delta_alpha = delta_alpha * 180. / std::f64::consts::PI;
    delta_alpha
}


fn topocentric_sun_right_ascension(alpha: &f64, delta_alpha: &f64) -> f64 {
    let delta_prime = alpha + delta_alpha;
    delta_prime
}


fn topocentric_sun_declination(delta: &f64, x: &f64, y: &f64, xi: &f64,
                               delta_alpha: &f64, h: &f64) -> f64 {
    let xi_rad = xi * std::f64::consts::PI / 180.;
    let h_rad = h * std::f64::consts::PI / 180.;
    let delta_rad = delta * std::f64::consts::PI / 180.;
    let delta_alpha_rad = delta_alpha * std::f64::consts::PI / 180.;
    let num = (f64::sin(delta_rad) - y * f64::sin(xi_rad)) * f64::cos(delta_alpha_rad);
    let denom = f64::cos(delta_rad) - x * f64::sin(xi_rad) * f64::cos(h_rad);
    let mut delta_prime = f64::atan2(num, denom);
    delta_prime = delta_prime * 180. / std::f64::consts::PI;
    delta_prime
}


fn topocentric_local_hour_angle(h: &f64, delta_alpha: &f64) -> f64 {
    let h_prime = h - delta_alpha;
    h_prime
}


fn topocentric_elevation_angle_without_atmosphere(lat: &f32, delta_prime: &f64,
                                                  h_prime: &f64) -> f64 {
    let lat_rad = *lat as f64 * std::f64::consts::PI / 180.;
    let delta_prime_rad = delta_prime * std::f64::consts::PI / 180.;
    let h_prime_rad = h_prime * std::f64::consts::PI / 180.;
    let mut e0: f64 = 0.;
    e0 = f64::asin(f64::sin(lat_rad) * f64::sin(delta_prime_rad)
                   + f64::cos(lat_rad) * f64::cos(delta_prime_rad) * f64::cos(h_prime_rad));
    e0 = e0 * 180. / std::f64::consts::PI;
    e0
}


fn atmospheric_refraction_correction(pressure: &f32, temp: &f32, e0: &f64,
                                     atmos_refract: &f32) -> f64 {
    if e0 >= &(-0.26667 - *atmos_refract as f64) {
        let term = (e0 + 10.3 / (e0 + 5.11)) * std::f64::consts::PI / 180.;
        let delta_e = *pressure as f64 / 1010.
                    * 283. / (273. + *temp as f64)
                    * 1.02 / (60. * f64::tan(term));
        return delta_e;
    }
    let delta_e = 0.;
    delta_e
}


fn topocentric_elevation_angle(e0: &f64, delta_e: &f64) -> f64 {
    let e = e0 + delta_e;
    e
}


fn topocentric_zenith_angle(e: &f64) -> f64 {
    let theta = 90. - e;
    theta
}


fn topocentric_astronomers_azimuth(h_prime: &f64, delta_prime: &f64, lat: &f32) -> f64 {
    let h_prime_rad = h_prime * std::f64::consts::PI / 180.;
    let delta_prime_rad = delta_prime * std::f64::consts::PI / 180.;
    let lat_rad = *lat as f64 * std::f64::consts::PI / 180.;
    let den = f64::cos(h_prime_rad) * f64::sin(lat_rad) - f64::tan(delta_prime_rad) * f64::cos(lat_rad);
    let mut gamma = f64::atan2(f64::sin(h_prime_rad), den);
    gamma = gamma * 180. / std::f64::consts::PI;
    limit_angle_range(&mut gamma);
    gamma
}


fn topocentric_azimuth_angle(gamma: &f64) -> f64 {
    let phi = gamma + 180.;
    ((phi % 360.) + 360.) % 360.
}


pub fn calculate_deltat(year: &i16, month: &u8) -> f64 {
    // https://eclipse.gsfc.nasa.gov/SEhelp/deltaT.html
    // https://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html
    if (year < &-1999) | (year > &3000) {
        panic!("`year` should be above -1999 and below 3000");
    }
    if (month < &1) | (month > &12) {
        panic!("`month` should be above 0 and below 13");
    }
    let y = *year as f64 + (*month as f64 - 0.5) / 12.;
    if year < &-500 {
        return -20. + 32. * f64::powf((y - 1820.) / 100., 2.);
    }
    else if (year >= &-500) & (year < &500) {
        return 10583.6 - 1014.41 * (y /100.)
               + 33.78311 * f64::powf(y / 100., 2.)
               - 5.952053 * f64::powf(y / 100., 3.)
               - 0.1798452 * f64::powf(y / 100., 4.)
               + 0.022174192 * f64::powf(y / 100., 5.)
               + 0.0090316521 * f64::powf(y / 100., 6.)
    }
    else if (year >= &500) & (year < &1600) {
        return 1574.2 - 556.01 * (y - 1000.) / 100.
               + 71.23472 * f64::powf((y - 1000.) / 100., 2.)
               + 0.319781 * f64::powf((y - 1000.) / 100., 3.)
               - 0.8503463 * f64::powf((y - 1000.) / 100., 4.)
               - 0.005050998 * f64::powf((y - 1000.) / 100., 5.)
               + 0.0083572073 * f64::powf((y - 1000.) / 100., 6.)
    }
    else if (year >= &1600) & (year < &1700) {
        return 120. - 0.9808 * (y - 1600.)
               - 0.01532 * f64::powf(y - 1600., 2.) 
               + f64::powf(y - 1600., 3.) / 7129.
    }
    else if (year >= &1700) & (year < &1800) {
        return 8.83 + 0.1603 * (y - 1700.)
               - 0.0059285 * f64::powf(y - 1700., 2.) 
               + 0.00013336 * f64::powf(y - 1700., 3.) 
               - f64::powf(y - 1700., 4.) / 1174000.
    }
    else if (year >= &1800) & (year < &1860) {
        return 13.72 - 0.332447 * (y - 1800.)
               + 0.0068612 * f64::powf(y - 1800., 2.) 
               + 0.0041116 * f64::powf(y - 1800., 3.) 
               - 0.00037436 * f64::powf(y - 1800., 4.) 
               + 0.0000121272 * f64::powf(y - 1800., 5.) 
               - 0.0000001699 * f64::powf(y - 1800., 6.) 
               + 0.000000000875 * f64::powf(y - 1800., 7.) 
    }
    else if (year >= &1860) & (year < &1900) {
        return 7.62 + 0.5737 * (y - 1860.)
               - 0.251754 * f64::powf(y - 1860., 2.) 
               + 0.01680668 * f64::powf(y - 1860., 3.) 
               - 0.0004473624 * f64::powf(y - 1860., 4.) 
               + f64::powf(y - 1860., 5.) / 233174.
    }
    else if (year >= &1900) & (year < &1920) {
        return -2.79 + 1.494119 * (y - 1900.)
               - 0.0598939 * f64::powf(y - 1900., 2.) 
               + 0.0061966 * f64::powf(y - 1900., 3.) 
               - 0.000197 * f64::powf(y - 1900., 4.) 
    }
    else if (year >= &1920) & (year < &1941) {
        return 21.2 + 0.84493 * (y - 1920.)
               - 0.0761 * f64::powf(y - 1920., 2.) 
               + 0.0020936 * f64::powf(y - 1920., 3.) 
    }
    else if (year >= &1941) & (year < &1961) {
        return 29.07 + 0.407 * (y - 1950.)
               - f64::powf(y - 1950., 2.) / 233.
               + f64::powf(y - 1950., 3.) / 2547.
    }
    else if (year >= &1961) & (year < &1986) {
        return 45.45 + 1.067 * (y - 1975.)
               - f64::powf(y - 1975., 2.) / 260.
               - f64::powf(y - 1975., 3.) / 718.
    }
    else if (year >= &1986) & (year < &2005) {
        return 63.86 + 0.3345 * (y - 2000.)
               - 0.060374 * f64::powf(y - 2000., 2.) 
               + 0.0017275 * f64::powf(y - 2000., 3.) 
               + 0.000651814 * f64::powf(y - 2000., 4.) 
               + 0.00002373599 * f64::powf(y - 2000., 5.) 
    }
    else if (year >= &2005) & (year < &2050) {
        return 62.92 + 0.32217 * (y - 2000.)
               + 0.005589 * f64::powf(y - 2000., 2.) 
    }
    else if (year >= &2050) & (year < &2150) {
        return -20. + 32. * f64::powf((y - 1820.) / 100., 2.) - 0.5628 * (2150. - y)
    }
    return -20. + 32. * f64::powf((y - 1820.) / 100., 2.) 
}


#[cfg(test)]
mod tests {
    use spa;

    #[test]
    fn test_calculate_deltat() {
        let expected = [17217.912272222224, 17211.72360555555,
                        17184.831552348543, 17177.310907656512,
                        118.96186861631928, 118.5375201220693,
                        8.99069606382023, 9.051572240380537,
                        13.385367077090038, 13.260910643977612,
                        7.962908405417306, 7.9713621281635225,
                        -1.2918431282701737, -0.720127316638437,
                        21.99992798755792, 22.27683824508097,
                        24.797268089496473, 25.035728072816564,
                        33.594798592945274, 33.75427944537937,
                        54.89627599023825, 55.079355914075165,
                        64.68633720312503, 64.84502657812497,
                        93.0847888888886, 93.93328888888917,
                        328.56800555555526, 329.4486722222225];
        let years = [-501, -499, 1601, 1701, 1801, 1861, 1901,
                     1921, 1941, 1961, 1986, 2005, 2050, 2150];
        for (i, year) in years.iter().enumerate() {
            for (j, month) in [1 as u8, 6 as u8].iter().enumerate() {
                println!("{} {} {} {}", year, month, i, j);
                assert_eq!(expected[i * 2 + j], spa::calculate_deltat(&year, &month));
            }
        }
    }

    #[test]
    #[should_panic]
    fn test_calculate_deltat_wrong_year_3001() {
        let year = 3001;
        let month = 1;
        spa::calculate_deltat(&year, &month);
    }

    #[test]
    #[should_panic]
    fn test_calculate_deltat_wrong_year_min_2000() {
        let year = -2000;
        let month = 1;
        spa::calculate_deltat(&year, &month);
    }

    #[test]
    #[should_panic]
    fn test_calculate_deltat_wrong_month_0() {
        let year = 2000;
        let month = 0;
        spa::calculate_deltat(&year, &month);
    }

    #[test]
    #[should_panic]
    fn test_calculate_deltat_wrong_month_13() {
        let year = 2000;
        let month = 13;
        spa::calculate_deltat(&year, &month);
    }

    #[test]
    fn test_julian_day() {
        let timestamp = 1357000200;
        assert_eq!(2456293.5208333335, spa::julian_day(&timestamp));
        // taken from A4.1:
        let timestamp = 582724800;
        assert_eq!(2447332.0, spa::julian_day(&timestamp));
    }

    #[test]
    fn test_julian_ephemeris_day() {
        let jd = 2456293.5208333335;
        let delta_t = 64.5;
        assert_eq!(2456293.521579861, spa::julian_ephemeris_day(&jd, &delta_t));
    }

    #[test]
    fn test_julian_century() {
        let jd = 2456293.5208333335;
        assert_eq!(0.13000741501255272, spa::julian_century(&jd));
    }

    #[test]
    fn test_julian_ephemeris_century() {
        let jde = 2456293.521579861;
        assert_eq!(0.1300074354513669, spa::julian_ephemeris_century(&jde));
    }

    #[test]
    fn test_julian_ephemeris_millennium() {
        let jde = 0.1300074354513669;
        assert_eq!(0.01300074354513669, spa::julian_ephemeris_millennium(&jde));
    }

    #[test]
    fn test_heliocentric_longitude() {
        let jme = 0.01300074354513669;
        // pvlib: 100.75286502047948
        assert_eq!(100.75288554139583, spa::heliocentric_longitude(&jme));
        let jme = -0.01300074354513669;
        // pvlib 100.00464976366129
        assert_eq!(100.00462189209429, spa::heliocentric_longitude(&jme));
    }

    #[test]
    fn test_heliocentric_latitude() {
        let jme = 0.01300074354513669;
        // pvlib 0.00010663876509755685
        assert_eq!(0.00010663836325029436, spa::heliocentric_latitude(&jme));
        let jme = -0.01300074354513669;
        // 0.00010214440296148753
        assert_eq!(0.00010214469246789543, spa::heliocentric_latitude(&jme));
    }

    #[test]
    fn test_heliocentric_radius_vector() {
        let jme = 0.01300074354513669;
        // pvlib 0.9832940565275712
        assert_eq!(0.983294088577116, spa::heliocentric_radius_vector(&jme));
        let jme = -0.01300074354513669;
        // pvlib 0.9833347325330166
        assert_eq!(0.9833347608615821, spa::heliocentric_radius_vector(&jme));
    }

    #[test]
    fn test_geocentric_longitude() {
        let l_deg = 90.;
        assert_eq!(270., spa::geocentric_longitude(&l_deg));
        let l_deg = 350.;
        assert_eq!(170.00000000000003, spa::geocentric_longitude(&l_deg));
    }

    #[test]
    fn test_geocentric_latitude() {
        let b_deg = 90.;
        assert_eq!(-90., spa::geocentric_latitude(&b_deg));
    }

    #[test]
    fn test_mean_elongation() {
        let jce = 0.13000741501255272;
        assert_eq!(58185.876481278865, spa::mean_elongation(&jce));
    }

    #[test]
    fn test_mean_anomaly_sun() {
        let jce = 0.13000741501255272;
        assert_eq!(5037.671194893454, spa::mean_anomaly_sun(&jce));
    }

    #[test]
    fn test_mean_anomaly_moon() {
        let jce = 0.13000741501255272;
        assert_eq!(62174.354324370404, spa::mean_anomaly_moon(&jce));
    }

    #[test]
    fn test_moon_argument_latitude() {
        let jce = 0.13000741501255272;
        assert_eq!(62913.117076730916, spa::moon_argument_latitude(&jce));
    }

    #[test]
    fn test_moon_ascending_longitude() {
        let jce = 0.13000741501255272;
        assert_eq!(-126.40750056925852, spa::moon_ascending_longitude(&jce));
    }

    #[test]
    fn test_longitude_obliquity_nutation() {
        let jce = 0.13000741501255272;
        let x0 = 58185.876481278865;
        let x1 = 5037.671194893454;
        let x2 = 62174.354324370404;
        let x3 = 62913.117076730916;
        let x4 = -126.40750056925852;
        let (n_longitude, n_obliquity) = spa::longitude_obliquity_nutation(&jce, &x0, &x1, &x2, &x3, &x4);
        // pvlib 0.004069674875572919
        assert_eq!(0.0040696748755729446, n_longitude);
        // pvlib -0.0016509592741862002
        assert_eq!(-0.0016509592741862074, n_obliquity);
    }

    #[test]
    fn test_mean_ecliptic_obliquity() {
        let jme = 0.01300074354513669;
        assert_eq!(84375.3624447249, spa::mean_ecliptic_obliquity(&jme));
    }

    #[test]
    fn test_true_ecliptic_obliquity() {
        let epsilon0 = 84375.3624447249;
        let delta_epsilon = -0.0016509592741862074;
        assert_eq!(23.435949719816065, spa::true_ecliptic_obliquity(&epsilon0, &delta_epsilon));
    }

    #[test]
    fn test_aberration_correction() {
        let r = 0.983294088577116;
        assert_eq!(-0.005788310107047633, spa::aberration_correction(&r));
    }

    #[test]
    fn test_apparent_sun_longitude() {
        let theta = 1.;
        let delta_epsilon = 1.;
        let delta_tau = 1.;
        assert_eq!(3., spa::apparent_sun_longitude(&theta, &delta_epsilon, &delta_tau));
    }

    #[test]
    fn test_mean_sidereal_time() {
        let jd = 2456293.5208333335;
        let jc = 0.13000741501255272;
        // pvlib 108.3276781309396
        assert_eq!(108.32767813109967, spa::mean_sidereal_time(&jd, &jc));
    }

    #[test]
    fn test_apparent_sidereal_time() {
        let v0: f64 = 108.32767813109967;
        let delta_psi: f64 = 0.0040696748755729446;
        let epsilon: f64 = 23.435949719816065;
        assert_eq!(108.33141207919714, spa::apparent_sidereal_time(&v0, &delta_psi, &epsilon));
    }

    #[test]
    fn test_geocentric_sun_right_ascension() {
        let lambda: f64 = 170.;
        let epsilon: f64 = 23.435949719816065;
        let beta: f64 = -0.0016509592741862074;
        assert_eq!(170.80960819870663, spa::geocentric_sun_right_ascension(&lambda, &epsilon, &beta));
    }

    #[test]
    fn test_geocentric_sun_declination() {
        let lambda: f64 = 170.;
        let epsilon: f64 = 23.435949719816065;
        let beta: f64 = -0.0016509592741862074;
        assert_eq!(3.9587091134966355, spa::geocentric_sun_declination(&lambda, &epsilon, &beta));
    }

    #[test]
    fn test_local_hour_angle() {
        let v: f64 = 170.;
        let lon: f32 = 100.;
        let alpha: f64 = 20.;
        assert_eq!(250., spa::local_hour_angle(&v, &lon, &alpha));
        let v: f64 = 170.;
        let lon: f32 = 250.;
        let alpha: f64 = 20.;
        assert_eq!(40.000000000000014, spa::local_hour_angle(&v, &lon, &alpha));
    }

    #[test]
    fn test_equatorial_horizontal_parallax() {
        let r: f64 = 0.983294088577116;
        assert_eq!(0.002484279938377968, spa::equatorial_horizontal_parallax(&r));
    }

    #[test]
    fn test_u_term() {
        let lat: f32 = 100.;
        assert_eq!(-1.3956881668836334, spa::u_term(&lat));
    }

    #[test]
    fn test_x_term() {
        let u: f64 = -1.3956881668836334;
        let elev: f32 = 200.;
        let lat: f32 = 100.;
        assert_eq!(0.17420919940606158, spa::x_term(&u, &elev, &lat));
    }

    #[test]
    fn test_y_term() {
        let u: f64 = -1.3956881668836334;
        let elev: f32 = 200.;
        let lat: f32 = 100.;
        assert_eq!(-0.9813752830759921, spa::y_term(&u, &elev, &lat));
    }

    #[test]
    fn test_parallax_sun_right_ascension() {
        let x: f64 = 0.17420919940606158;
        let xi: f64 = 0.002484279938377968;
        let h: f64 = 40.;
        let delta: f64 = 3.9587091134966355;
        // pvlib -0.0002788554074644773
        assert_eq!(-0.0002788554074644773, spa::parallax_sun_right_ascension(&x, &xi, &h, &delta))
    }

    #[test]
    fn test_topocentric_sun_right_ascension() {
        let alpha: f64 = 50.;
        let delta_alpha: f64 = 50.;
        assert_eq!(100., spa::topocentric_sun_right_ascension(&alpha, &delta_alpha));
    }

    #[test]
    fn test_topocentric_sun_declination() {
        let delta: f64 = 3.9587091134966355;
        let x: f64 = 0.17420919940606158;
        let y: f64 = -0.9813752830759921;
        let xi: f64 = 0.002484279938377968;
        let delta_alpha: f64 = 50.;
        let h: f64 = 40.;
        assert_eq!(2.548569550651887, spa::topocentric_sun_declination(
            &delta, &x, &y, &xi, &delta_alpha, &h            
        ));
    }

    #[test]
    fn test_topocentric_local_hour_angle() {
        let h: f64 = 50.;
        let delta_alpha: f64 = 30.;
        assert_eq!(20., spa::topocentric_local_hour_angle(&h, &delta_alpha));
    }

    #[test]
    fn test_topocentric_elevation_angle_without_atmosphere() {
        let lat: f32 = 100.;
        let delta_prime = 2.548569550651887;
        let h_prime: f64 = 20.;
        assert_eq!(-6.847307469841585, spa::topocentric_elevation_angle_without_atmosphere(
            &lat, &delta_prime, &h_prime));
    }

    #[test]
    fn test_atmospheric_refraction_correction() {
        let pressure: f32 = 1000.;
        let temp: f32 = 20.;
        let e0: f64 = -6.847307469841585;
        let atmos_refract: f32 = 2.;
        assert_eq!(0., spa::atmospheric_refraction_correction(&pressure,
            &temp, &e0, &atmos_refract));
        let e0: f64 = 6.847307469841585;
        assert_eq!(0.12010357702796845, spa::atmospheric_refraction_correction(&pressure,
            &temp, &e0, &atmos_refract));
    }

    #[test]
    fn test_topocentric_elevation_angle() {
        let e0: f64 = 1.;
        let delta_e: f64 = 1.;
        assert_eq!(2., spa::topocentric_elevation_angle(&e0, &delta_e));
    }

    #[test]
    fn test_topocentric_zenith_angle() {
        let e: f64 = 45.;
        assert_eq!(45., spa::topocentric_zenith_angle(&e));
    }

    #[test]
    fn test_topocentric_astronomers_azimuth() {
        let h_prime: f64 = 20.;
        let delta_prime = 2.548569550651887;
        let lat: f32 = 100.;
        assert_eq!(20.129089268987894,
                   spa::topocentric_astronomers_azimuth(&h_prime, &delta_prime, &lat));
    }

    #[test]
    fn test_topocentric_azimuth_angle() {
        let gamma = 20.;
        assert_eq!(200., spa::topocentric_azimuth_angle(&gamma));
        let gamma = 220.;
        assert_eq!(40., spa::topocentric_azimuth_angle(&gamma));
        let gamma = -200.;
        assert_eq!(340., spa::topocentric_azimuth_angle(&gamma));
    }

} 
