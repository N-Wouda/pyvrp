#ifndef CROSSOVER_H
#define CROSSOVER_H

#include "CostEvaluator.h"
#include "Individual.h"
#include "ProblemData.h"
#include "XorShift128.h"

#include <functional>
#include <vector>

namespace crossover
{
/**
 * Greedily inserts each unplanned client into the non-empty route that's
 * nearest to the client.
 */
void greedyRepair(std::vector<std::vector<int>> &routes,
                  std::vector<int> const &unplanned,
                  ProblemData const &data,
                  CostEvaluator const &costEvaluator);
}  // namespace crossover

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
 * @param parents          The parent individuals.
 * @param data             The problem data.
 * @param costEvaluator    The cost evaluator.
 * @param startIndices     Start indices of routes in parent individuals.
 * @param numMovedRoutes   Number of routes to move.
 * @return A new offspring.
 *
 * <br />
 * Yuichi Nagata and Shigenobu Kobayashi. "A memetic algorithm for
 * the pickup and delivery problem with time windows using selective route
 * exchange crossover". In: International Conference on Parallel Problem Solving
 * from Nature. Springer. 2010, pp. 536–545.
 */
Individual selectiveRouteExchange(
    std::pair<Individual const *, Individual const *> const &parents,
    ProblemData const &data,
    CostEvaluator const &costEvaluator,
    std::pair<size_t, size_t> const startIndices,
    size_t const numMovedRoutes);

#endif  // CROSSOVER_H
