#include "SwapRoutes.h"

#include <vector>

using pyvrp::Cost;
using pyvrp::search::Route;
using pyvrp::search::SwapRoutes;
using TWS = pyvrp::TimeWindowSegment;

Cost SwapRoutes::evaluate(Route *U,
                          Route *V,
                          CostEvaluator const &costEvaluator)
{
    if (U->vehicleType() == V->vehicleType() || U->empty() || V->empty())
        return 0;

    Cost deltaCost = 0;

    // Changes in load capacity violations.
    deltaCost += costEvaluator.loadPenalty(U->load(), V->capacity());
    deltaCost -= costEvaluator.loadPenalty(U->load(), U->capacity());

    deltaCost += costEvaluator.loadPenalty(V->load(), U->capacity());
    deltaCost -= costEvaluator.loadPenalty(V->load(), V->capacity());

    // Changes in distance from and to the depot.
    // TODO

    // Changes in time warp.
    // TODO

    return deltaCost;
}

void SwapRoutes::apply([[maybe_unused]] Route *U,
                       [[maybe_unused]] Route *V) const
{
    std::vector<Route::Node *> uNodes = {U->begin(), U->end()};
    std::vector<Route::Node *> vNodes = {V->begin(), V->end()};

    // We are swapping the routes completely, so we can just clear both.
    U->clear();
    V->clear();

    // Reinsert the nodes from route U in route V.
    for (auto *node : uNodes)
        V->push_back(node);

    // Reinsert the nodes from route V in route U.
    for (auto *node : vNodes)
        U->push_back(node);
}

SwapRoutes::SwapRoutes(ProblemData const &data)
    : LocalSearchOperator<Route>(data)
{
}
