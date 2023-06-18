use criterion::{black_box, criterion_group, criterion_main, Criterion};
use solpos::spa;
use std::error::Error;

pub fn criterion_benchmark(c: &mut Criterion) -> Result<(), Box<dyn Error>> {
    let mut g = c.benchmark_group("Solar Position");

    // Expansive functions
    g.bench_function("heliocentric_longitude", |b| {
        b.iter(|| black_box(spa::heliocentric_longitude(0.01300074354513669)))
    });

    g.bench_function("heliocentric_latitude", |b| {
        b.iter(|| black_box(spa::heliocentric_latitude(0.01300074354513669)))
    });

    g.bench_function("heliocentric_radius_vector", |b| {
        b.iter(|| black_box(spa::heliocentric_radius_vector(0.01300074354513669)))
    });

    g.bench_function("longitude_obliquity_nutation", |b| {
        b.iter(|| {
            black_box(spa::longitude_obliquity_nutation(
                0.13000741501255272,
                58185.876481278865,
                5037.671194893454,
                62174.354324370404,
                62913.117076730916,
                -126.40750056925852,
            ))
        })
    });

    g.bench_function("calculate_delta_t", |b| {
        b.iter(|| black_box(spa::calculate_deltat(2005., 6.)))
    });

    // Entrypoint
    let mut builder = spa::SolPosInputsBuilder::default();
    builder.unixtimes(vec![1538939951.; 8760 * 60]);
    builder.latitude(30.29);
    builder.longitude(-97.74);
    builder.elevation(213.);
    builder.pressure(101325.);
    builder.temperature(12.);
    builder.atmos_refract(0.5667);
    builder.parallel_calcs(true);
    let inputs = builder.build().unwrap();

    g.bench_function("solpos", |b| {
        b.iter(|| black_box(spa::calculate_solar_position(inputs.clone())))
    });

    g.finish();

    Ok(())
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
