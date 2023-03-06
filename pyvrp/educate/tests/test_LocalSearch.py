from numpy.testing import assert_, assert_equal, assert_raises
from pytest import mark

from pyvrp import Individual, PenaltyManager, XorShift128
from pyvrp.educate import LocalSearch, NeighbourhoodParams, compute_neighbours
from pyvrp.tests.helpers import read


def test_empty_neighbourhood_search_returns_same_solution():
    data = read("data/OkSmall.txt")
    pm = PenaltyManager(data.vehicle_capacity)
    rng = XorShift128(seed=42)

    neighbourhood = [[] for _ in range(data.num_clients + 1)]  # is empty
    ls = LocalSearch(data, pm, rng, neighbourhood)

    for _ in range(100):  # repeat a few times to make sure
        individual = Individual.make_random(data, pm, rng)
        assert_equal(ls.search(individual), individual)


def test_no_node_operators_search_returns_same_solution():
    data = read("data/OkSmall.txt")
    pm = PenaltyManager(data.vehicle_capacity)
    rng = XorShift128(seed=42)
    ls = LocalSearch(data, pm, rng, compute_neighbours(data))

    for _ in range(100):  # repeat a few times to make sure
        individual = Individual.make_random(data, pm, rng)
        assert_equal(ls.search(individual), individual)


def test_no_route_operators_intensify_returns_same_solution():
    data = read("data/OkSmall.txt")
    pm = PenaltyManager(data.vehicle_capacity)
    rng = XorShift128(seed=42)
    ls = LocalSearch(data, pm, rng, compute_neighbours(data))

    for _ in range(100):  # repeat a few times to make sure
        individual = Individual.make_random(data, pm, rng)
        assert_equal(ls.intensify(individual), individual)


@mark.parametrize("size", [1, 2, 3, 4, 6, 7])  # num_clients + 1 == 5
def test_local_search_raises_when_neighbourhood_dimensions_do_not_match(size):
    data = read("data/OkSmall.txt")
    pm = PenaltyManager(data.vehicle_capacity)
    rng = XorShift128(seed=42)

    # Each of the given sizes is either smaller than or bigger than desired.
    neighbours = [[] for _ in range(size)]

    with assert_raises(RuntimeError):
        LocalSearch(data, pm, rng, neighbours)

    ls = LocalSearch(data, pm, rng, compute_neighbours(data))

    with assert_raises(RuntimeError):
        ls.set_neighbours(neighbours)


def test_local_search_raises_when_neighbourhood_contains_self_or_depot():
    data = read("data/OkSmall.txt")
    pm = PenaltyManager(data.vehicle_capacity)
    rng = XorShift128(seed=42)

    neighbours = [[client] for client in range(data.num_clients + 1)]

    with assert_raises(RuntimeError):
        LocalSearch(data, pm, rng, neighbours)


@mark.parametrize(
    "weight_wait_time,"
    "weight_time_warp,"
    "nb_granular,"
    "symmetric_proximity,"
    "symmetric_neighbours",
    [
        (20, 20, 10, True, False),
        (20, 20, 10, True, True),
        # From original c++ implementation
        # (18, 20, 34, False),
        (18, 20, 34, True, True),
    ],
)
def test_local_search_set_get_neighbours(
    weight_wait_time: int,
    weight_time_warp: int,
    nb_granular: int,
    symmetric_proximity: bool,
    symmetric_neighbours: bool,
):
    data = read("data/RC208.txt", "solomon", round_func="trunc")

    seed = 42
    rng = XorShift128(seed=seed)
    pen_manager = PenaltyManager(data.vehicle_capacity)

    params = NeighbourhoodParams(nb_granular=1)
    prev_neighbours = compute_neighbours(data, params)
    ls = LocalSearch(data, pen_manager, rng, prev_neighbours)

    params = NeighbourhoodParams(
        weight_wait_time,
        weight_time_warp,
        nb_granular,
        symmetric_proximity,
        symmetric_neighbours,
    )
    neighbours = compute_neighbours(data, params)

    # Test that before we set neighbours we don't have same
    assert_(ls.get_neighbours() != neighbours)

    # Test after we set we have the same
    ls.set_neighbours(neighbours)
    ls_neighbours = ls.get_neighbours()
    assert_equal(ls_neighbours, neighbours)

    # Check that the bindings make a copy (in both directions)
    assert_(ls_neighbours is not neighbours)
    ls_neighbours[1] = []
    assert_(ls.get_neighbours() != ls_neighbours)
    assert_equal(ls.get_neighbours(), neighbours)
    neighbours[1] = []
    assert_(ls.get_neighbours() != neighbours)
