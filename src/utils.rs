use chrono::{Datelike, NaiveDateTime};

pub fn get_year_and_month_from_unixtime(unixtime: i64) -> (f64, f64) {
    let datetime = NaiveDateTime::from_timestamp_opt(unixtime, 0).unwrap(); // FIXME
    let year = datetime.year();
    let month = datetime.month();
    (year as f64, month as f64)
}
