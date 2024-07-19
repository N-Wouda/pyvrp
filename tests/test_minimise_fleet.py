from numpy.testing import assert_equal

from pyvrp import minimise_fleet
from pyvrp.stop import MaxRuntime
from tests.helpers import read


def test_OkSmall(ok_small):
    """
    TODO
    """
    assert_equal(ok_small.num_vehicles, 3)

    veh_types = minimise_fleet(ok_small, MaxRuntime(0.25))
    data = ok_small.replace(vehicle_types=veh_types)
    assert_equal(data.num_vehicles, 2)


def test_rc208(rc208):
    assert_equal(rc208.num_vehicles, 25)

    veh_types = minimise_fleet(rc208, MaxRuntime(1))
    data = rc208.replace(vehicle_types=veh_types)
    assert_equal(data.num_vehicles, 4)


def test_X_instance():
    data = read("data/X-n101-50-k13.vrp", round_func="round")
    assert_equal(data.num_vehicles, 100)

    veh_types = minimise_fleet(data, MaxRuntime(1))
    data = data.replace(vehicle_types=veh_types)
    assert_equal(data.num_vehicles, 13)