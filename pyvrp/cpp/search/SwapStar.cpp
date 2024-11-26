#include "SwapStar.h"

#include <cassert>

using pyvrp::Cost;
using pyvrp::Load;
using pyvrp::search::Route;
using pyvrp::search::SwapStar;

void SwapStar::updateRemovalCosts(Route *R, CostEvaluator const &costEvaluator)
{
    isCached(R->idx(), 0) = true;

    for (size_t idx = 1; idx != R->size() + 1; ++idx)
    {
        auto const proposal
            = R->proposal(R->before(idx - 1), R->after(idx + 1));

        Cost deltaCost = 0;
        costEvaluator.deltaCost<true, true>(deltaCost, proposal);

        auto const *U = (*R)[idx];
        removalCosts(R->idx(), U->client()) = deltaCost;
    }

    for (size_t idx = data.numDepots(); idx != data.numLocations(); ++idx)
        isCached(R->idx(), idx) = false;
}

void SwapStar::updateInsertPoints(Route *R,
                                  Route::Node *U,
                                  CostEvaluator const &costEvaluator)
{
    auto &insertPositions = insertCache(R->idx(), U->client());
    insertPositions.fill({std::numeric_limits<Cost>::max(), nullptr});

    for (size_t idx = 0; idx != R->size() + 1; ++idx)
    {
        auto const proposal = R->proposal(
            R->before(idx), U->route()->at(U->idx()), R->after(idx + 1));

        Cost deltaCost = 0;
        costEvaluator.deltaCost<true, true>(deltaCost, proposal);

        auto *V = (*R)[idx];

        if (deltaCost >= insertPositions[2].first)
            continue;

        if (deltaCost >= insertPositions[1].first)
            insertPositions[2] = {deltaCost, V};
        else if (deltaCost >= insertPositions[0].first)
        {
            insertPositions[2] = insertPositions[1];
            insertPositions[1] = {deltaCost, V};
        }
        else
        {
            insertPositions[2] = insertPositions[1];
            insertPositions[1] = insertPositions[0];
            insertPositions[0] = {deltaCost, V};
        }
    }

    isCached(R->idx(), U->client()) = true;
}

Cost SwapStar::deltaLoadCost(Route::Node *U,
                             Route::Node *V,
                             CostEvaluator const &costEvaluator) const
{
    auto const *uRoute = U->route();
    auto const *vRoute = V->route();

    ProblemData::Client const &uClient = data.location(U->client());
    ProblemData::Client const &vClient = data.location(V->client());

    auto const &uLoad = uRoute->load();
    auto const &uCap = uRoute->capacity();

    auto const &vLoad = vRoute->load();
    auto const &vCap = vRoute->capacity();

    // Separating removal and insertion means that the effects on load are not
    // counted correctly: during insert, U is still in the route, and now V is
    // added as well. The following addresses this issue with an approximation,
    // which is inexact when there are both pickups and deliveries in the data.
    // So it's pretty rough but fast and seems to mostly work well enough.
    Cost delta = 0;
    for (size_t dim = 0; dim != data.numLoadDimensions(); ++dim)
    {
        auto const loadDiff
            = std::max(uClient.delivery[dim], uClient.pickup[dim])
              - std::max(vClient.delivery[dim], vClient.pickup[dim]);

        delta += costEvaluator.loadPenalty(uLoad[dim] - loadDiff, uCap[dim]);
        delta -= costEvaluator.loadPenalty(uLoad[dim], uCap[dim]);

        delta += costEvaluator.loadPenalty(vLoad[dim] + loadDiff, vCap[dim]);
        delta -= costEvaluator.loadPenalty(vLoad[dim], vCap[dim]);
    }

    return delta;
}

SwapStar::InsertPoint SwapStar::bestInsertPoint(
    Route::Node *U, Route::Node *V, CostEvaluator const &costEvaluator)
{
    auto *route = V->route();
    auto &best_ = insertCache(route->idx(), U->client());

    if (!isCached(route->idx(), U->client()))
        updateInsertPoints(route, U, costEvaluator);

    for (auto [cost, where] : best_)  // only OK if V is not adjacent
        if (where && where != V && n(where) != V)
            return std::make_pair(cost, where);

    // As a fallback option, we consider inserting in the place of V.
    Cost deltaCost = 0;
    costEvaluator.deltaCost<true, true>(
        deltaCost,
        route->proposal(route->before(V->idx() - 1),
                        U->route()->at(U->idx()),
                        route->after(V->idx() + 1)));

    return std::make_pair(deltaCost, p(V));
}

Cost SwapStar::evaluateMove(Route::Node const *U,
                            Route::Node const *V,
                            Route::Node const *remove,
                            CostEvaluator const &costEvaluator) const
{
    assert(V->route() == remove->route());
    assert(V != remove);

    auto const *route = V->route();

    Cost deltaCost = 0;

    if (V->idx() + 1 == remove->idx())  // then we insert U in place of remove
        costEvaluator.deltaCost<true>(
            deltaCost,
            route->proposal(route->before(V->idx()),
                            U->route()->at(U->idx()),
                            route->after(V->idx() + 2)));
    else if (V->idx() < remove->idx())
        costEvaluator.deltaCost<true>(
            deltaCost,
            route->proposal(route->before(V->idx()),
                            U->route()->at(U->idx()),
                            route->between(V->idx() + 1, remove->idx() - 1),
                            route->after(remove->idx() + 1)));
    else if (V->idx() > remove->idx())
        costEvaluator.deltaCost<true>(
            deltaCost,
            route->proposal(route->before(remove->idx() - 1),
                            route->between(remove->idx() + 1, V->idx()),
                            U->route()->at(U->idx()),
                            route->after(V->idx() + 1)));

    return deltaCost;
}

void SwapStar::init(Solution const &solution)
{
    LocalSearchOperator<Route>::init(solution);
    for (size_t row = 0; row != isCached.numRows(); ++row)
        isCached(row, 0) = false;
}

Cost SwapStar::evaluate(Route *routeU,
                        Route *routeV,
                        CostEvaluator const &costEvaluator)
{
    best = {};

    if (!isCached(routeU->idx(), 0))
        updateRemovalCosts(routeU, costEvaluator);

    if (!isCached(routeV->idx(), 0))
        updateRemovalCosts(routeV, costEvaluator);

    for (auto *U : *routeU)
        for (auto *V : *routeV)
        {
            // The following lines compute a delta cost of removing U and V from
            // their own routes and inserting them into the other's route in the
            // best place. This is approximate since removal and insertion are
            // evaluated separately, not taking into account that while U leaves
            // its route, V will be inserted (and vice versa).
            Cost deltaCost = 0;

            // Load is a bit tricky, so we compute that separately.
            deltaCost += deltaLoadCost(U, V, costEvaluator);

            deltaCost += removalCosts(routeU->idx(), U->client());
            deltaCost += removalCosts(routeV->idx(), V->client());

            auto [extraV, UAfter] = bestInsertPoint(U, V, costEvaluator);
            deltaCost += extraV;

            if (deltaCost >= 0)  // continuing here avoids evaluating another
                continue;        // costly insertion point below

            auto [extraU, VAfter] = bestInsertPoint(V, U, costEvaluator);
            deltaCost += extraU;

            if (deltaCost < best.cost)
            {
                best.cost = deltaCost;

                best.U = U;
                best.UAfter = UAfter;

                best.V = V;
                best.VAfter = VAfter;
            }
        }

    // It is possible for positive delta costs to turn negative when we do an
    // exact evaluation. But in practice that almost never happens, and is not
    // worth spending time on.
    if (best.cost >= 0)
        return best.cost;

    return evaluateMove(best.V, best.VAfter, best.U, costEvaluator)
           + evaluateMove(best.U, best.UAfter, best.V, costEvaluator);
}

void SwapStar::apply(Route *U, Route *V) const
{
    assert(best.U);
    assert(best.UAfter);
    assert(best.V);
    assert(best.VAfter);

    U->remove(best.U->idx());
    V->remove(best.V->idx());

    V->insert(best.UAfter->idx() + 1, best.U);
    U->insert(best.VAfter->idx() + 1, best.V);
}

void SwapStar::update(Route *U) { isCached(U->idx(), 0) = false; }

SwapStar::SwapStar(ProblemData const &data)
    : LocalSearchOperator<Route>(data),
      insertCache(data.numVehicles(), data.numLocations()),
      isCached(data.numVehicles(), data.numLocations()),
      removalCosts(data.numVehicles(), data.numLocations())
{
}
