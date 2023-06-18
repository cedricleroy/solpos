import pandas as pd

import solpos


def test_simple():
    daterange = pd.date_range(
        "2018-01-01 00:00:00", "2018-12-31 23:59:00", inclusive="left", freq="T", tz="Etc/GMT+1"
    )
    solpos.solar_position(
        (daterange.astype("int64") // 1e9),
        {"lat": 30.29, "lon": -97.74, "elev": 213, "parallel_calcs": True},
    )
