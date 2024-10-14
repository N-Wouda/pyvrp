#include "bindings.h"
#include "Exchange.h"
#include "LocalSearch.h"
#include "Route.h"
#include "SwapRoutes.h"
#include "SwapStar.h"
#include "SwapTails.h"
#include "primitives.h"
#include "search_docs.h"

#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

#include <sstream>

namespace nb = nanobind;

using pyvrp::search::Exchange;
using pyvrp::search::inplaceCost;
using pyvrp::search::insertCost;
using pyvrp::search::LocalSearch;
using pyvrp::search::LocalSearchOperator;
using pyvrp::search::removeCost;
using pyvrp::search::Route;
using pyvrp::search::SwapRoutes;
using pyvrp::search::SwapStar;
using pyvrp::search::SwapTails;

NB_MODULE(_search, m)
{
    using NodeOp = LocalSearchOperator<pyvrp::search::Route::Node>;
    using RouteOp = LocalSearchOperator<pyvrp::search::Route>;

    nb::class_<NodeOp>(m, "NodeOperator");
    nb::class_<RouteOp>(m, "RouteOperator");

    nb::class_<Exchange<1, 0>, NodeOp>(
        m, "Exchange10", DOC(pyvrp, search, Exchange))
        .def(nb::init<pyvrp::ProblemData const &>(),
             nb::arg("data"),
             nb::keep_alive<1, 2>())  // keep data alive
        .def("evaluate",
             &Exchange<1, 0>::evaluate,
             nb::arg("U"),
             nb::arg("V"),
             nb::arg("cost_evaluator"))
        .def("apply", &Exchange<1, 0>::apply, nb::arg("U"), nb::arg("V"));

    nb::class_<Exchange<2, 0>, NodeOp>(
        m, "Exchange20", DOC(pyvrp, search, Exchange))
        .def(nb::init<pyvrp::ProblemData const &>(),
             nb::arg("data"),
             nb::keep_alive<1, 2>())  // keep data alive
        .def("evaluate",
             &Exchange<2, 0>::evaluate,
             nb::arg("U"),
             nb::arg("V"),
             nb::arg("cost_evaluator"))
        .def("apply", &Exchange<2, 0>::apply, nb::arg("U"), nb::arg("V"));

    nb::class_<Exchange<3, 0>, NodeOp>(
        m, "Exchange30", DOC(pyvrp, search, Exchange))
        .def(nb::init<pyvrp::ProblemData const &>(),
             nb::arg("data"),
             nb::keep_alive<1, 2>())  // keep data alive
        .def("evaluate",
             &Exchange<3, 0>::evaluate,
             nb::arg("U"),
             nb::arg("V"),
             nb::arg("cost_evaluator"))
        .def("apply", &Exchange<3, 0>::apply, nb::arg("U"), nb::arg("V"));

    nb::class_<Exchange<1, 1>, NodeOp>(
        m, "Exchange11", DOC(pyvrp, search, Exchange))
        .def(nb::init<pyvrp::ProblemData const &>(),
             nb::arg("data"),
             nb::keep_alive<1, 2>())  // keep data alive
        .def("evaluate",
             &Exchange<1, 1>::evaluate,
             nb::arg("U"),
             nb::arg("V"),
             nb::arg("cost_evaluator"))
        .def("apply", &Exchange<1, 1>::apply, nb::arg("U"), nb::arg("V"));

    nb::class_<Exchange<2, 1>, NodeOp>(
        m, "Exchange21", DOC(pyvrp, search, Exchange))
        .def(nb::init<pyvrp::ProblemData const &>(),
             nb::arg("data"),
             nb::keep_alive<1, 2>())  // keep data alive
        .def("evaluate",
             &Exchange<2, 1>::evaluate,
             nb::arg("U"),
             nb::arg("V"),
             nb::arg("cost_evaluator"))
        .def("apply", &Exchange<2, 1>::apply, nb::arg("U"), nb::arg("V"));

    nb::class_<Exchange<3, 1>, NodeOp>(
        m, "Exchange31", DOC(pyvrp, search, Exchange))
        .def(nb::init<pyvrp::ProblemData const &>(),
             nb::arg("data"),
             nb::keep_alive<1, 2>())  // keep data alive
        .def("evaluate",
             &Exchange<3, 1>::evaluate,
             nb::arg("U"),
             nb::arg("V"),
             nb::arg("cost_evaluator"))
        .def("apply", &Exchange<3, 1>::apply, nb::arg("U"), nb::arg("V"));

    nb::class_<Exchange<2, 2>, NodeOp>(
        m, "Exchange22", DOC(pyvrp, search, Exchange))
        .def(nb::init<pyvrp::ProblemData const &>(),
             nb::arg("data"),
             nb::keep_alive<1, 2>())  // keep data alive
        .def("evaluate",
             &Exchange<2, 2>::evaluate,
             nb::arg("U"),
             nb::arg("V"),
             nb::arg("cost_evaluator"))
        .def("apply", &Exchange<2, 2>::apply, nb::arg("U"), nb::arg("V"));

    nb::class_<Exchange<3, 2>, NodeOp>(
        m, "Exchange32", DOC(pyvrp, search, Exchange))
        .def(nb::init<pyvrp::ProblemData const &>(),
             nb::arg("data"),
             nb::keep_alive<1, 2>())  // keep data alive
        .def("evaluate",
             &Exchange<3, 2>::evaluate,
             nb::arg("U"),
             nb::arg("V"),
             nb::arg("cost_evaluator"))
        .def("apply", &Exchange<3, 2>::apply, nb::arg("U"), nb::arg("V"));

    nb::class_<Exchange<3, 3>, NodeOp>(
        m, "Exchange33", DOC(pyvrp, search, Exchange))
        .def(nb::init<pyvrp::ProblemData const &>(),
             nb::arg("data"),
             nb::keep_alive<1, 2>())  // keep data alive
        .def("evaluate",
             &Exchange<3, 3>::evaluate,
             nb::arg("U"),
             nb::arg("V"),
             nb::arg("cost_evaluator"))
        .def("apply", &Exchange<3, 3>::apply, nb::arg("U"), nb::arg("V"));

    nb::class_<SwapRoutes, RouteOp>(
        m, "SwapRoutes", DOC(pyvrp, search, SwapRoutes))
        .def(nb::init<pyvrp::ProblemData const &>(),
             nb::arg("data"),
             nb::keep_alive<1, 2>())  // keep data alive
        .def("evaluate",
             &SwapRoutes::evaluate,
             nb::arg("U"),
             nb::arg("V"),
             nb::arg("cost_evaluator"))
        .def("apply", &SwapRoutes::apply, nb::arg("U"), nb::arg("V"));

    nb::class_<SwapStar, RouteOp>(m, "SwapStar", DOC(pyvrp, search, SwapStar))
        .def(nb::init<pyvrp::ProblemData const &>(),
             nb::arg("data"),
             nb::keep_alive<1, 2>())  // keep data alive
        .def("evaluate",
             &SwapStar::evaluate,
             nb::arg("U"),
             nb::arg("V"),
             nb::arg("cost_evaluator"))
        .def("apply", &SwapStar::apply, nb::arg("U"), nb::arg("V"));

    nb::class_<SwapTails, NodeOp>(m, "SwapTails", DOC(pyvrp, search, SwapTails))
        .def(nb::init<pyvrp::ProblemData const &>(),
             nb::arg("data"),
             nb::keep_alive<1, 2>())  // keep data alive
        .def("evaluate",
             &SwapTails::evaluate,
             nb::arg("U"),
             nb::arg("V"),
             nb::arg("cost_evaluator"))
        .def("apply", &SwapTails::apply, nb::arg("U"), nb::arg("V"));

    nb::class_<LocalSearch>(m, "LocalSearch")
        .def(nb::init<pyvrp::ProblemData const &,
                      std::vector<std::vector<size_t>>>(),
             nb::arg("data"),
             nb::arg("neighbours"),
             nb::keep_alive<1, 2>())  // keep data alive until LS is freed
        .def("add_node_operator",
             &LocalSearch::addNodeOperator,
             nb::arg("op"),
             nb::keep_alive<1, 2>())
        .def("add_route_operator",
             &LocalSearch::addRouteOperator,
             nb::arg("op"),
             nb::keep_alive<1, 2>())
        .def("set_neighbours",
             &LocalSearch::setNeighbours,
             nb::arg("neighbours"))
        .def("neighbours",
             &LocalSearch::neighbours,
             nb::rv_policy::reference_internal)
        .def("__call__",
             &LocalSearch::operator(),
             nb::arg("solution"),
             nb::arg("cost_evaluator"))
        .def("search",
             nb::overload_cast<pyvrp::Solution const &,
                               pyvrp::CostEvaluator const &>(
                 &LocalSearch::search),
             nb::arg("solution"),
             nb::arg("cost_evaluator"))
        .def("intensify",
             nb::overload_cast<pyvrp::Solution const &,
                               pyvrp::CostEvaluator const &,
                               double const>(&LocalSearch::intensify),
             nb::arg("solution"),
             nb::arg("cost_evaluator"),
             nb::arg("overlap_tolerance") = 0.05)
        .def("shuffle", &LocalSearch::shuffle, nb::arg("rng"));

    nb::class_<Route>(m, "Route", DOC(pyvrp, search, Route))
        .def(nb::init<pyvrp::ProblemData const &, size_t, size_t>(),
             nb::arg("data"),
             nb::arg("idx"),
             nb::arg("vehicle_type"),
             nb::keep_alive<1, 2>())  // keep data alive
        .def_prop_ro("idx", &Route::idx)
        .def_prop_ro("vehicle_type", &Route::vehicleType)
        .def("__delitem__", &Route::remove, nb::arg("idx"))
        .def("__getitem__",
             &Route::operator[],
             nb::arg("idx"),
             nb::rv_policy::reference_internal)
        .def(
            "__iter__",
            [](Route const &route)
            {
                return nb::make_iterator(
                    nb::type<Route>(), "iterator", route.begin(), route.end());
            },
            nb::rv_policy::reference_internal)
        .def("__len__", &Route::size)
        .def("__str__",
             [](Route const &route)
             {
                 std::stringstream stream;
                 stream << route;
                 return stream.str();
             })
        .def("is_feasible", &Route::isFeasible)
        .def("has_excess_load", &Route::hasExcessLoad)
        .def("has_excess_distance", &Route::hasExcessDistance)
        .def("has_time_warp", &Route::hasTimeWarp)
        .def("capacity", &Route::capacity)
        .def("start_depot", &Route::startDepot)
        .def("end_depot", &Route::endDepot)
        .def("fixed_vehicle_cost", &Route::fixedVehicleCost)
        .def("load", &Route::load)
        .def("excess_load", &Route::excessLoad)
        .def("excess_distance", &Route::excessDistance)
        .def("distance", &Route::distance)
        .def("distance_cost", &Route::distanceCost)
        .def("unit_distance_cost", &Route::unitDistanceCost)
        .def("duration", &Route::duration)
        .def("duration_cost", &Route::durationCost)
        .def("unit_duration_cost", &Route::unitDurationCost)
        .def("max_duration", &Route::maxDuration)
        .def("max_distance", &Route::maxDistance)
        .def("time_warp", &Route::timeWarp)
        .def("profile", &Route::profile)
        .def(
            "dist_at",
            [](Route const &route, size_t idx, size_t profile)
            { return route.at(idx).distance(profile); },
            nb::arg("idx"),
            nb::arg("profile") = 0)
        .def(
            "dist_between",
            [](Route const &route, size_t start, size_t end, size_t profile)
            { return route.between(start, end).distance(profile); },
            nb::arg("start"),
            nb::arg("end"),
            nb::arg("profile") = 0)
        .def(
            "dist_after",
            [](Route const &route, size_t start, size_t profile)
            { return route.after(start).distance(profile); },
            nb::arg("start"),
            nb::arg("profile") = 0)
        .def(
            "dist_before",
            [](Route const &route, size_t end, size_t profile)
            { return route.before(end).distance(profile); },
            nb::arg("end"),
            nb::arg("profile") = 0)
        .def(
            "load_at",
            [](Route const &route, size_t idx) { return route.at(idx).load(); },
            nb::arg("idx"))
        .def(
            "load_between",
            [](Route const &route, size_t start, size_t end)
            { return route.between(start, end).load(); },
            nb::arg("start"),
            nb::arg("end"))
        .def(
            "load_after",
            [](Route const &route, size_t start)
            { return route.after(start).load(); },
            nb::arg("start"))
        .def(
            "load_before",
            [](Route const &route, size_t end)
            { return route.before(end).load(); },
            nb::arg("end"))
        .def(
            "duration_at",
            [](Route const &route, size_t idx, size_t profile)
            { return route.at(idx).duration(profile); },
            nb::arg("idx"),
            nb::arg("profile") = 0)
        .def(
            "duration_between",
            [](Route const &route, size_t start, size_t end, size_t profile)
            { return route.between(start, end).duration(profile); },
            nb::arg("start"),
            nb::arg("end"),
            nb::arg("profile") = 0)
        .def(
            "duration_after",
            [](Route const &route, size_t start, size_t profile)
            { return route.after(start).duration(profile); },
            nb::arg("start"),
            nb::arg("profile") = 0)
        .def(
            "duration_before",
            [](Route const &route, size_t end, size_t profile)
            { return route.before(end).duration(profile); },
            nb::arg("end"),
            nb::arg("profile") = 0)
        .def("centroid", &Route::centroid)
        .def("overlaps_with",
             &Route::overlapsWith,
             nb::arg("other"),
             nb::arg("tolerance"))
        .def("append",
             &Route::push_back,
             nb::arg("node"),
             nb::keep_alive<1, 2>(),  // keep node alive
             nb::keep_alive<2, 1>())  // keep route alive
        .def("clear", &Route::clear)
        .def("insert",
             &Route::insert,
             nb::arg("idx"),
             nb::arg("node"),
             nb::keep_alive<1, 3>(),  // keep node alive
             nb::keep_alive<3, 1>())  // keep route alive
        .def("update", &Route::update);

    nb::class_<Route::Node>(m, "Node", DOC(pyvrp, search, Route, Node))
        .def(nb::init<size_t>(), nb::arg("loc"))
        .def_prop_ro("client", &Route::Node::client)
        .def_prop_ro("idx", &Route::Node::idx)
        .def_prop_ro("route", &Route::Node::route)
        .def("is_depot", &Route::Node::isDepot);

    m.def("insert_cost",
          &insertCost,
          nb::arg("U"),
          nb::arg("V"),
          nb::arg("data"),
          nb::arg("cost_evaluator"),
          DOC(pyvrp, search, insertCost));

    m.def("inplace_cost",
          &inplaceCost,
          nb::arg("U"),
          nb::arg("V"),
          nb::arg("data"),
          nb::arg("cost_evaluator"),
          DOC(pyvrp, search, inplaceCost));

    m.def("remove_cost",
          &removeCost,
          nb::arg("U"),
          nb::arg("data"),
          nb::arg("cost_evaluator"),
          DOC(pyvrp, search, removeCost));
}
