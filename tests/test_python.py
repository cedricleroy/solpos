import numpy as np
import numpy.testing as npt

from solpos import spa_scalar, spa_vector


def test_spa():
    res = spa_scalar(
        1538939951,
        2018,
        10,
        30.29,
        -97.74,
        213,
        1008.6,
        23,
        0.5667
    )
    expected = {
        'theta': 38.788075669522726,
        'theta0': 38.80104099275841,
        'e': 51.211924330477274,
        'e0': 51.19895900724159,
        'phi': 204.46709671428397
    }
    npt.assert_almost_equal(
        list(expected.values()), list(res.values())
    )


def test_vector():
    inputs = []
    expected = {
        'theta': 38.788075669522726,
        'theta0': 38.80104099275841,
        'e': 51.211924330477274,
        'e0': 51.19895900724159,
        'phi': 204.46709671428397
    }
    for x in range(1000):
        inputs.append([
            1538939951,
            2018,
            10,
            30.29,
            -97.74,
            213,
            1008.6,
            23,
            0.5667
        ])
    inputs = np.asarray(inputs)
    inputs = inputs.astype(float)
    res = spa_vector(inputs)
    npt.assert_almost_equal(
        list(expected.values()), list(res[2])
    )
    assert res.shape == (1000, 5)
