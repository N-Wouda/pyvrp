from numpy.testing import assert_, assert_allclose, assert_equal, assert_raises
from pytest import mark

from pyvrp import (
    GeneticAlgorithm,
    GeneticAlgorithmParams,
    Individual,
    PenaltyManager,
    Population,
    PopulationParams,
    XorShift128,
)
from pyvrp.crossover import selective_route_exchange as srex
from pyvrp.diversity import broken_pairs_distance as bpd
from pyvrp.educate import Exchange10, LocalSearch, compute_neighbours
from pyvrp.stop import MaxIterations
from pyvrp.tests.helpers import make_random_solutions, read


@mark.parametrize(
    "repair_probability,"
    "collect_statistics,"
    "intensify_probability,"
    "intensify_on_best,"
    "nb_iter_no_improvement",
    [
        (-0.25, True, 0.5, True, 0),  # repair_probability < 0
        (1.25, True, 0.5, True, 0),  # repair_probability > 1
        (0.0, True, 0.5, True, -1),  # nb_iter_no_improvement < 0
    ],
)
def test_params_constructor_throws_when_arguments_invalid(
    repair_probability: float,
    collect_statistics: bool,
    intensify_probability: int,
    intensify_on_best: bool,
    nb_iter_no_improvement: int,
):
    """
    Tests that invalid configurations are not accepted.
    """
    with assert_raises(ValueError):
        GeneticAlgorithmParams(
            repair_probability,
            collect_statistics,
            intensify_probability,
            intensify_on_best,
            nb_iter_no_improvement,
        )


@mark.parametrize(
    "repair_probability,"
    "collect_statistics,"
    "intensify_probability,"
    "intensify_on_best,"
    "nb_iter_no_improvement",
    [
        (0.0, True, 0.5, True, 0),  # nb_iter_no_improvement == 0
        (0.0, True, 0.5, True, 1),  # repair_probability == 0
        (1.0, True, 0.5, True, 1),  # repair_probability == 1
        (0.5, False, 0.5, True, 1),  # collect_statistics is False
        (0.5, True, 0, True, 1),  # intensify_probability == 0
        (0.5, True, 1, True, 1),  # intensify_probability == 1
        (0.5, True, 0.5, False, 1),  # intensify_on_best is False
        (0.5, False, 0.5, False, 1),  # both False
    ],
)
def test_params_constructor_does_not_raise_when_arguments_valid(
    repair_probability: float,
    collect_statistics: bool,
    intensify_probability: float,
    intensify_on_best: bool,
    nb_iter_no_improvement: int,
):
    """
    Tests valid boundary cases.
    """
    params = GeneticAlgorithmParams(
        repair_probability,
        collect_statistics,
        intensify_probability,
        intensify_on_best,
        nb_iter_no_improvement,
    )

    assert_allclose(params.repair_probability, repair_probability)
    assert_equal(params.collect_statistics, collect_statistics)
    assert_equal(params.intensify_probability, intensify_probability)
    assert_equal(params.intensify_on_best, intensify_on_best)
    assert_equal(params.nb_iter_no_improvement, nb_iter_no_improvement)


def test_raises_when_too_small_population_and_no_initial_solutions():
    """
    Tests that GeneticAlgorithm rejects empty populations with no provided
    initial solutions, since that is insufficient to do crossover.
    """
    data = read("data/RC208.txt", "solomon", "dimacs")
    pen_manager = PenaltyManager(data.vehicle_capacity)
    rng = XorShift128(seed=42)
    ls = LocalSearch(data, pen_manager, rng, compute_neighbours(data))

    pop = Population(bpd)
    assert_equal(len(pop), 0)

    with assert_raises(ValueError):
        # No individuals should raise.
        GeneticAlgorithm(data, pen_manager, rng, pop, ls, srex)

    individual = Individual.make_random(data, pen_manager, rng)

    # We have provided an initial solution, so this should be OK.
    GeneticAlgorithm(data, pen_manager, rng, pop, ls, srex, [individual])

    pop.add(individual)
    assert_equal(len(pop), 1)

    # One individual in the population without providing initial solutions
    # should also be OK.
    GeneticAlgorithm(data, pen_manager, rng, pop, ls, srex)


def test_initial_solutions():
    """
    Tests that GeneticAlgorithm adds initial solutions to the population
    when running.
    """
    data = read("data/E-n22-k4.txt", round_func="round")
    pm = PenaltyManager(data.vehicle_capacity)
    rng = XorShift128(seed=42)
    pop = Population(bpd)
    ls = LocalSearch(data, pm, rng, compute_neighbours(data))
    init = [Individual.make_random(data, pm, rng) for _ in range(25)]
    algo = GeneticAlgorithm(data, pm, rng, pop, ls, srex, init)

    algo.run(MaxIterations(0))

    # Check that the initial population individuals have the same routes as the
    # initial solutions.
    current = {individual for individual in pop}
    assert_equal(len(current & set(init)), 25)


def test_restart():
    """
    Tests that GeneticAlgorithm upon restarting clears the population and
    adds the initial solutions.
    """
    data = read("data/E-n22-k4.txt", round_func="round")
    pm = PenaltyManager(data.vehicle_capacity)
    rng = XorShift128(seed=42)
    pop = Population(bpd)

    ls = LocalSearch(data, pm, rng, compute_neighbours(data))
    ls.add_node_operator(Exchange10(data, pm))

    init = [Individual.make_random(data, pm, rng) for _ in range(25)]
    params = GeneticAlgorithmParams(
        repair_probability=0,
        intensify_probability=0,
        intensify_on_best=False,
        nb_iter_no_improvement=100,
    )
    algo = GeneticAlgorithm(data, pm, rng, pop, ls, srex, init, params=params)

    algo.run(MaxIterations(100))

    # Check that the population contains the initial solutions, and one more
    # due to the education step without repair.
    current = {individual for individual in pop}
    assert_equal(len(current & set(init)), 25)
    assert_equal(len(pop), 26)


def test_best_solution_improves_with_more_iterations():
    data = read("data/RC208.txt", "solomon", "dimacs")
    rng = XorShift128(seed=42)
    pm = PenaltyManager(data.vehicle_capacity)

    pop_params = PopulationParams()
    pop = Population(bpd, params=pop_params)

    for indiv in make_random_solutions(pop_params.min_pop_size, data, pm, rng):
        pop.add(indiv)

    ls = LocalSearch(data, pm, rng, compute_neighbours(data))
    ls.add_node_operator(Exchange10(data, pm))

    ga_params = GeneticAlgorithmParams(
        intensify_probability=0, intensify_on_best=False
    )
    algo = GeneticAlgorithm(data, pm, rng, pop, ls, srex, params=ga_params)

    initial_best = algo.run(MaxIterations(0)).best
    new_best = algo.run(MaxIterations(25)).best

    assert_(new_best.cost() < initial_best.cost())
    assert_(new_best.is_feasible())  # best must be feasible


# TODO more functional tests

# TODO test statistics collection on Result.has_statistics
