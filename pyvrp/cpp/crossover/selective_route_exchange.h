#ifndef PYVRP_SELECTIVE_ROUTE_EXCHANGE_H
#define PYVRP_SELECTIVE_ROUTE_EXCHANGE_H

#include "CostEvaluator.h"
#include "DynamicBitset.h"
#include "ProblemData.h"
#include "Solution.h"

#include <functional>
#include <vector>

namespace pyvrp::crossover
{
/**
 * Performs two SREX crossovers of the given parents. SREX is a method that
 * selects a set of routes for each parent and replaces the selected routes of
 * the first parent with those of the second parent. The routes are selected by
 * minimizing the overlap between the two sets of routes. This is achieved
 * through a heuristic that iteratively shifts adjacent routes until no further
 * improvement in minimizing the overlap is observed. Then, two offspring are
 * generated by replacing the selected routes in two distinct ways, and the
 * offspring with the lowest cost is returned.
 *
 * @param parents          The parent solutions.
 * @param data             The problem data.
 * @param costEvaluator    The cost evaluator.
 * @param startIndices     Start indices of routes in parent solutions.
 * @param numMovedRoutes   Number of routes to move.
 * @return A new offspring.
 *
 * Note that this is an internal docstring: the SREX operator is wrapped on
 * the Python side.
 */
Solution selectiveRouteExchange(
    std::pair<Solution const *, Solution const *> const &parents,
    ProblemData const &data,
    CostEvaluator const &costEvaluator,
    std::pair<size_t, size_t> const startIndices,
    size_t const numMovedRoutes);
}  // namespace pyvrp::crossover

#endif  // PYVRP_SELECTIVE_ROUTE_EXCHANGE_H
