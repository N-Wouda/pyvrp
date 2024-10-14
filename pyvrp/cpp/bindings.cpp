#include "bindings.h"
#include "CostEvaluator.h"
#include "DistanceSegment.h"
#include "DurationSegment.h"
#include "DynamicBitset.h"
#include "LoadSegment.h"
#include "Matrix.h"
#include "ProblemData.h"
#include "RandomNumberGenerator.h"
#include "Route.h"
#include "Solution.h"
#include "SubPopulation.h"
#include "pyvrp_docs.h"

#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <sstream>
#include <string>
#include <variant>

namespace nb = nanobind;

using pyvrp::CostEvaluator;
using pyvrp::DistanceSegment;
using pyvrp::DurationSegment;
using pyvrp::DynamicBitset;
using pyvrp::LoadSegment;
using pyvrp::Matrix;
using pyvrp::PopulationParams;
using pyvrp::ProblemData;
using pyvrp::RandomNumberGenerator;
using pyvrp::Route;
using pyvrp::Solution;
using pyvrp::SubPopulation;

NB_MODULE(_pyvrp, m)
{
    nb::class_<DynamicBitset>(m, "DynamicBitset", DOC(pyvrp, DynamicBitset))
        .def(nb::init<size_t>(), nb::arg("num_bits"))
        .def(nb::self == nb::self, nb::arg("other"))  // this is __eq__
        .def("all", &DynamicBitset::all)
        .def("any", &DynamicBitset::any)
        .def("none", &DynamicBitset::none)
        .def("count", &DynamicBitset::count)
        .def("__len__", &DynamicBitset::size)
        .def("reset", &DynamicBitset::reset)
        .def(
            "__getitem__",
            [](DynamicBitset const &bitset, size_t idx) { return bitset[idx]; },
            nb::arg("idx"))
        .def(
            "__setitem__",
            [](DynamicBitset &bitset, size_t idx, bool value)
            { bitset[idx] = value; },
            nb::arg("idx"),
            nb::arg("value"))
        .def("__or__", &DynamicBitset::operator|, nb::arg("other"))
        .def("__and__", &DynamicBitset::operator&, nb::arg("other"))
        .def("__xor__", &DynamicBitset::operator^, nb::arg("other"))
        .def("__invert__", &DynamicBitset::operator~);

    nb::class_<ProblemData::Client>(
        m, "Client", DOC(pyvrp, ProblemData, Client))
        .def(nb::init<pyvrp::Coordinate,
                      pyvrp::Coordinate,
                      pyvrp::Load,
                      pyvrp::Load,
                      pyvrp::Duration,
                      pyvrp::Duration,
                      pyvrp::Duration,
                      pyvrp::Duration,
                      pyvrp::Cost,
                      bool,
                      std::optional<size_t>,
                      char const *>(),
             nb::arg("x"),
             nb::arg("y"),
             nb::arg("delivery") = 0,
             nb::arg("pickup") = 0,
             nb::arg("service_duration") = 0,
             nb::arg("tw_early") = 0,
             nb::arg("tw_late") = std::numeric_limits<pyvrp::Duration>::max(),
             nb::arg("release_time") = 0,
             nb::arg("prize") = 0,
             nb::arg("required") = true,
             nb::arg("group") = nb::none(),
             nb::kw_only(),
             nb::arg("name") = "")
        .def_ro("x", &ProblemData::Client::x)
        .def_ro("y", &ProblemData::Client::y)
        .def_ro("delivery", &ProblemData::Client::delivery)
        .def_ro("pickup", &ProblemData::Client::pickup)
        .def_ro("service_duration", &ProblemData::Client::serviceDuration)
        .def_ro("tw_early", &ProblemData::Client::twEarly)
        .def_ro("tw_late", &ProblemData::Client::twLate)
        .def_ro("release_time", &ProblemData::Client::releaseTime)
        .def_ro("prize", &ProblemData::Client::prize)
        .def_ro("required", &ProblemData::Client::required)
        .def_ro("group", &ProblemData::Client::group)
        .def_ro("name",
                &ProblemData::Client::name,
                nb::rv_policy::reference_internal)
        .def(nb::self == nb::self)  // this is __eq__
        .def("__getstate__",
             [](ProblemData::Client const &client)
             {
                 return std::make_tuple(client.x,
                                        client.y,
                                        client.delivery,
                                        client.pickup,
                                        client.serviceDuration,
                                        client.twEarly,
                                        client.twLate,
                                        client.releaseTime,
                                        client.prize,
                                        client.required,
                                        client.group,
                                        client.name);
             })
        .def("__setstate__",
             [](ProblemData::Client &client,
                std::tuple<pyvrp::Coordinate,
                           pyvrp::Coordinate,
                           pyvrp::Load,
                           pyvrp::Load,
                           pyvrp::Duration,
                           pyvrp::Duration,
                           pyvrp::Duration,
                           pyvrp::Duration,
                           pyvrp::Cost,
                           bool,
                           std::optional<size_t>,
                           std::string> const &t)
             {
                 new (&client)
                     ProblemData::Client(std::get<0>(t),    // x
                                         std::get<1>(t),    // y
                                         std::get<2>(t),    // delivery
                                         std::get<3>(t),    // pickup
                                         std::get<4>(t),    // service duration
                                         std::get<5>(t),    // tw early
                                         std::get<6>(t),    // tw late
                                         std::get<7>(t),    // release time
                                         std::get<8>(t),    // prize
                                         std::get<9>(t),    // required
                                         std::get<10>(t),   // group
                                         std::get<11>(t));  // name
             })
        .def(
            "__str__",
            [](ProblemData::Client const &client) { return client.name; },
            nb::rv_policy::reference_internal);

    nb::class_<ProblemData::Depot>(m, "Depot", DOC(pyvrp, ProblemData, Depot))
        .def(nb::init<pyvrp::Coordinate, pyvrp::Coordinate, char const *>(),
             nb::arg("x"),
             nb::arg("y"),
             nb::kw_only(),
             nb::arg("name") = "")
        .def_ro("x", &ProblemData::Depot::x)
        .def_ro("y", &ProblemData::Depot::y)
        .def_ro("name",
                &ProblemData::Depot::name,
                nb::rv_policy::reference_internal)
        .def(nb::self == nb::self)  // this is __eq__
        .def("__getstate__",
             [](ProblemData::Depot const &depot)
             { return std::make_tuple(depot.x, depot.y, depot.name); })
        .def("__setstate__",
             [](ProblemData::Depot &depot,
                std::tuple<pyvrp::Coordinate,
                           pyvrp::Coordinate,
                           std::string> const &t)
             {
                 new (&depot) ProblemData::Depot(std::get<0>(t),   // x
                                                 std::get<1>(t),   // y
                                                 std::get<2>(t));  // name
             })
        .def(
            "__str__",
            [](ProblemData::Depot const &depot) { return depot.name; },
            nb::rv_policy::reference_internal);

    nb::class_<ProblemData::ClientGroup>(
        m, "ClientGroup", DOC(pyvrp, ProblemData, ClientGroup))
        .def(nb::init<std::vector<size_t>, bool>(),
             nb::arg("clients") = nb::list(),
             nb::arg("required") = true)
        .def("add_client",
             &ProblemData::ClientGroup::addClient,
             nb::arg("client"))
        .def("clear", &ProblemData::ClientGroup::clear)
        .def_prop_ro("clients", &ProblemData::ClientGroup::clients)
        .def_ro("required", &ProblemData::ClientGroup::required)
        .def_ro("mutually_exclusive",
                &ProblemData::ClientGroup::mutuallyExclusive)
        .def(nb::self == nb::self)  // this is __eq__
        .def("__getstate__",
             [](ProblemData::ClientGroup const &group)
             { return std::make_tuple(group.clients(), group.required); })
        .def("__setstate__",
             [](ProblemData::ClientGroup &group,
                std::tuple<std::vector<size_t>, bool> const &t)
             {
                 new (&group)
                     ProblemData::ClientGroup(std::get<0>(t),   // clients
                                              std::get<1>(t));  // required
             })
        .def("__len__", &ProblemData::ClientGroup::size)
        .def(
            "__iter__",
            [](ProblemData::ClientGroup const &group)
            {
                return nb::make_iterator(nb::type<ProblemData::ClientGroup>(),
                                         "iterator",
                                         group.begin(),
                                         group.end());
            },
            nb::rv_policy::reference_internal);

    nb::class_<ProblemData::VehicleType>(
        m, "VehicleType", DOC(pyvrp, ProblemData, VehicleType))
        .def(nb::init<size_t,
                      pyvrp::Load,
                      size_t,
                      size_t,
                      pyvrp::Cost,
                      pyvrp::Duration,
                      pyvrp::Duration,
                      pyvrp::Duration,
                      pyvrp::Distance,
                      pyvrp::Cost,
                      pyvrp::Cost,
                      size_t,
                      char const *>(),
             nb::arg("num_available") = 1,
             nb::arg("capacity") = 0,
             nb::arg("start_depot") = 0,
             nb::arg("end_depot") = 0,
             nb::arg("fixed_cost") = 0,
             nb::arg("tw_early") = 0,
             nb::arg("tw_late") = std::numeric_limits<pyvrp::Duration>::max(),
             nb::arg("max_duration")
             = std::numeric_limits<pyvrp::Duration>::max(),
             nb::arg("max_distance")
             = std::numeric_limits<pyvrp::Distance>::max(),
             nb::arg("unit_distance_cost") = 1,
             nb::arg("unit_duration_cost") = 0,
             nb::arg("profile") = 0,
             nb::kw_only(),
             nb::arg("name") = "")
        .def_ro("num_available", &ProblemData::VehicleType::numAvailable)
        .def_ro("capacity", &ProblemData::VehicleType::capacity)
        .def_ro("start_depot", &ProblemData::VehicleType::startDepot)
        .def_ro("end_depot", &ProblemData::VehicleType::endDepot)
        .def_ro("fixed_cost", &ProblemData::VehicleType::fixedCost)
        .def_ro("tw_early", &ProblemData::VehicleType::twEarly)
        .def_ro("tw_late", &ProblemData::VehicleType::twLate)
        .def_ro("max_duration", &ProblemData::VehicleType::maxDuration)
        .def_ro("max_distance", &ProblemData::VehicleType::maxDistance)
        .def_ro("unit_distance_cost",
                &ProblemData::VehicleType::unitDistanceCost)
        .def_ro("unit_duration_cost",
                &ProblemData::VehicleType::unitDurationCost)
        .def_ro("profile", &ProblemData::VehicleType::profile)
        .def_ro("name",
                &ProblemData::VehicleType::name,
                nb::rv_policy::reference_internal)
        .def("replace",
             &ProblemData::VehicleType::replace,
             nb::arg("num_available") = nb::none(),
             nb::arg("capacity") = nb::none(),
             nb::arg("start_depot") = nb::none(),
             nb::arg("end_depot") = nb::none(),
             nb::arg("fixed_cost") = nb::none(),
             nb::arg("tw_early") = nb::none(),
             nb::arg("tw_late") = nb::none(),
             nb::arg("max_duration") = nb::none(),
             nb::arg("max_distance") = nb::none(),
             nb::arg("unit_distance_cost") = nb::none(),
             nb::arg("unit_duration_cost") = nb::none(),
             nb::arg("profile") = nb::none(),
             nb::kw_only(),
             nb::arg("name") = nb::none(),
             DOC(pyvrp, ProblemData, VehicleType, replace))
        .def(nb::self == nb::self)  // this is __eq__
        .def("__getstate__",
             [](ProblemData::VehicleType const &vehicleType)
             {
                 return std::make_tuple(vehicleType.numAvailable,
                                        vehicleType.capacity,
                                        vehicleType.startDepot,
                                        vehicleType.endDepot,
                                        vehicleType.fixedCost,
                                        vehicleType.twEarly,
                                        vehicleType.twLate,
                                        vehicleType.maxDuration,
                                        vehicleType.maxDistance,
                                        vehicleType.unitDistanceCost,
                                        vehicleType.unitDurationCost,
                                        vehicleType.profile,
                                        vehicleType.name);
             })
        .def("__setstate__",
             [](ProblemData::VehicleType &vehicleType,
                std::tuple<size_t,
                           pyvrp::Load,
                           size_t,
                           size_t,
                           pyvrp::Cost,
                           pyvrp::Duration,
                           pyvrp::Duration,
                           pyvrp::Duration,
                           pyvrp::Distance,
                           pyvrp::Cost,
                           pyvrp::Cost,
                           size_t,
                           std::string> const &t)
             {
                 new (&vehicleType) ProblemData::VehicleType(
                     std::get<0>(t),    // num available
                     std::get<1>(t),    // capacity
                     std::get<2>(t),    // start depot
                     std::get<3>(t),    // end depot
                     std::get<4>(t),    // fixed cost
                     std::get<5>(t),    // tw early
                     std::get<6>(t),    // tw late
                     std::get<7>(t),    // max duration
                     std::get<8>(t),    // max distance
                     std::get<9>(t),    // unit distance cost
                     std::get<10>(t),   // unit duration cost
                     std::get<11>(t),   // profile
                     std::get<12>(t));  // name
             })
        .def(
            "__str__",
            [](ProblemData::VehicleType const &vehType)
            { return vehType.name; },
            nb::rv_policy::reference_internal);

    nb::class_<ProblemData>(m, "ProblemData", DOC(pyvrp, ProblemData))
        .def(nb::init<std::vector<ProblemData::Client>,
                      std::vector<ProblemData::Depot>,
                      std::vector<ProblemData::VehicleType>,
                      std::vector<Matrix<pyvrp::Distance>>,
                      std::vector<Matrix<pyvrp::Duration>>,
                      std::vector<ProblemData::ClientGroup>>(),
             nb::arg("clients"),
             nb::arg("depots"),
             nb::arg("vehicle_types"),
             nb::arg("distance_matrices"),
             nb::arg("duration_matrices"),
             nb::arg("groups") = nb::list())
        .def("replace",
             &ProblemData::replace,
             nb::arg("clients") = nb::none(),
             nb::arg("depots") = nb::none(),
             nb::arg("vehicle_types") = nb::none(),
             nb::arg("distance_matrices") = nb::none(),
             nb::arg("duration_matrices") = nb::none(),
             nb::arg("groups") = nb::none(),
             DOC(pyvrp, ProblemData, replace))
        .def_prop_ro("num_clients",
                     &ProblemData::numClients,
                     DOC(pyvrp, ProblemData, numClients))
        .def_prop_ro("num_depots",
                     &ProblemData::numDepots,
                     DOC(pyvrp, ProblemData, numDepots))
        .def_prop_ro("num_groups",
                     &ProblemData::numGroups,
                     DOC(pyvrp, ProblemData, numGroups))
        .def_prop_ro("num_locations",
                     &ProblemData::numLocations,
                     DOC(pyvrp, ProblemData, numLocations))
        .def_prop_ro("num_vehicle_types",
                     &ProblemData::numVehicleTypes,
                     DOC(pyvrp, ProblemData, numVehicleTypes))
        .def_prop_ro("num_vehicles",
                     &ProblemData::numVehicles,
                     DOC(pyvrp, ProblemData, numVehicles))
        .def_prop_ro("num_profiles",
                     &ProblemData::numProfiles,
                     DOC(pyvrp, ProblemData, numProfiles))
        .def(
            "location",
            [](ProblemData const &data,
               size_t idx) -> std::variant<ProblemData::Client const *,
                                           ProblemData::Depot const *>
            {
                if (idx >= data.numLocations())
                    throw nb::index_error();

                auto const proxy = data.location(idx);
                if (idx < data.numDepots())
                    return proxy.depot;
                else
                    return proxy.client;
            },
            nb::arg("idx"),
            nb::rv_policy::reference_internal,
            DOC(pyvrp, ProblemData, location))
        .def("clients",
             &ProblemData::clients,
             nb::rv_policy::reference_internal,
             DOC(pyvrp, ProblemData, clients))
        .def("depots",
             &ProblemData::depots,
             nb::rv_policy::reference_internal,
             DOC(pyvrp, ProblemData, depots))
        .def("groups",
             &ProblemData::groups,
             nb::rv_policy::reference_internal,
             DOC(pyvrp, ProblemData, groups))
        .def("vehicle_types",
             &ProblemData::vehicleTypes,
             nb::rv_policy::reference_internal,
             DOC(pyvrp, ProblemData, vehicleTypes))
        .def("distance_matrices",
             &ProblemData::distanceMatrices,
             nb::rv_policy::reference_internal,
             DOC(pyvrp, ProblemData, distanceMatrices))
        .def("duration_matrices",
             &ProblemData::durationMatrices,
             nb::rv_policy::reference_internal,
             DOC(pyvrp, ProblemData, durationMatrices))
        .def("centroid",
             &ProblemData::centroid,
             nb::rv_policy::reference_internal,
             DOC(pyvrp, ProblemData, centroid))
        .def("group",
             &ProblemData::group,
             nb::arg("group"),
             nb::rv_policy::reference_internal,
             DOC(pyvrp, ProblemData, group))
        .def("vehicle_type",
             &ProblemData::vehicleType,
             nb::arg("vehicle_type"),
             nb::rv_policy::reference_internal,
             DOC(pyvrp, ProblemData, vehicleType))
        .def("distance_matrix",
             &ProblemData::distanceMatrix,
             nb::arg("profile"),
             nb::rv_policy::reference_internal,
             DOC(pyvrp, ProblemData, distanceMatrix))
        .def("duration_matrix",
             &ProblemData::durationMatrix,
             nb::arg("profile"),
             nb::rv_policy::reference_internal,
             DOC(pyvrp, ProblemData, durationMatrix))
        .def(nb::self == nb::self)  // this is __eq__
        .def("__getstate__",
             [](ProblemData const &data)
             {
                 return std::make_tuple(data.clients(),
                                        data.depots(),
                                        data.vehicleTypes(),
                                        data.distanceMatrices(),
                                        data.durationMatrices(),
                                        data.groups());
             })
        .def("__setstate__",
             [](ProblemData &data,
                std::tuple<std::vector<ProblemData::Client>,
                           std::vector<ProblemData::Depot>,
                           std::vector<ProblemData::VehicleType>,
                           std::vector<pyvrp::Matrix<pyvrp::Distance>>,
                           std::vector<pyvrp::Matrix<pyvrp::Duration>>,
                           std::vector<ProblemData::ClientGroup>> const &t)
             {
                 new (&data) ProblemData(std::get<0>(t),   // clients
                                         std::get<1>(t),   // depots
                                         std::get<2>(t),   // vehicle types
                                         std::get<3>(t),   // distances
                                         std::get<4>(t),   // durations
                                         std::get<5>(t));  // groups
             });

    nb::class_<Route>(m, "Route", DOC(pyvrp, Route))
        .def(nb::init<ProblemData const &, std::vector<size_t>, size_t>(),
             nb::arg("data"),
             nb::arg("visits"),
             nb::arg("vehicle_type"))
        .def("visits",
             &Route::visits,
             nb::rv_policy::reference_internal,
             DOC(pyvrp, Route, visits))
        .def("distance", &Route::distance, DOC(pyvrp, Route, distance))
        .def("distance_cost",
             &Route::distanceCost,
             DOC(pyvrp, Route, distanceCost))
        .def("excess_distance",
             &Route::excessDistance,
             DOC(pyvrp, Route, excessDistance))
        .def("delivery", &Route::delivery, DOC(pyvrp, Route, delivery))
        .def("pickup", &Route::pickup, DOC(pyvrp, Route, pickup))
        .def("excess_load", &Route::excessLoad, DOC(pyvrp, Route, excessLoad))
        .def("duration", &Route::duration, DOC(pyvrp, Route, duration))
        .def("duration_cost",
             &Route::durationCost,
             DOC(pyvrp, Route, durationCost))
        .def("time_warp", &Route::timeWarp, DOC(pyvrp, Route, timeWarp))
        .def("start_time", &Route::startTime, DOC(pyvrp, Route, startTime))
        .def("end_time", &Route::endTime, DOC(pyvrp, Route, endTime))
        .def("slack", &Route::slack, DOC(pyvrp, Route, slack))
        .def("travel_duration",
             &Route::travelDuration,
             DOC(pyvrp, Route, travelDuration))
        .def("service_duration",
             &Route::serviceDuration,
             DOC(pyvrp, Route, serviceDuration))
        .def("wait_duration",
             &Route::waitDuration,
             DOC(pyvrp, Route, waitDuration))
        .def(
            "release_time", &Route::releaseTime, DOC(pyvrp, Route, releaseTime))
        .def("prizes", &Route::prizes, DOC(pyvrp, Route, prizes))
        .def("centroid", &Route::centroid, DOC(pyvrp, Route, centroid))
        .def(
            "vehicle_type", &Route::vehicleType, DOC(pyvrp, Route, vehicleType))
        .def("start_depot", &Route::startDepot, DOC(pyvrp, Route, startDepot))
        .def("end_depot", &Route::endDepot, DOC(pyvrp, Route, endDepot))
        .def("is_feasible", &Route::isFeasible, DOC(pyvrp, Route, isFeasible))
        .def("has_excess_load",
             &Route::hasExcessLoad,
             DOC(pyvrp, Route, hasExcessLoad))
        .def("has_excess_distance",
             &Route::hasExcessDistance,
             DOC(pyvrp, Route, hasExcessDistance))
        .def("has_time_warp",
             &Route::hasTimeWarp,
             DOC(pyvrp, Route, hasTimeWarp))
        .def("__len__", &Route::size, DOC(pyvrp, Route, size))
        .def(
            "__iter__",
            [](Route const &route)
            {
                return nb::make_iterator(
                    nb::type<Route>(), "iterator", route.begin(), route.end());
            },
            nb::rv_policy::reference_internal)
        .def(
            "__getitem__",
            [](Route const &route, int idx)
            {
                // int so we also support negative offsets from the end.
                idx = idx < 0 ? route.size() + idx : idx;
                if (idx < 0 || static_cast<size_t>(idx) >= route.size())
                    throw nb::index_error();
                return route[idx];
            },
            nb::arg("idx"))
        .def(nb::self == nb::self)  // this is __eq__
        .def("__getstate__",
             [](Route const &route)
             {
                 // Returns a tuple that completely encodes the route's state.
                 return std::make_tuple(route.visits(),
                                        route.distance(),
                                        route.distanceCost(),
                                        route.excessDistance(),
                                        route.delivery(),
                                        route.pickup(),
                                        route.excessLoad(),
                                        route.duration(),
                                        route.durationCost(),
                                        route.timeWarp(),
                                        route.travelDuration(),
                                        route.serviceDuration(),
                                        route.waitDuration(),
                                        route.releaseTime(),
                                        route.startTime(),
                                        route.slack(),
                                        route.prizes(),
                                        route.centroid(),
                                        route.vehicleType(),
                                        route.startDepot(),
                                        route.endDepot());
             })
        .def("__setstate__",
             [](Route &route,
                std::tuple<std::vector<size_t>,
                           pyvrp::Distance,
                           pyvrp::Cost,
                           pyvrp::Distance,
                           pyvrp::Load,
                           pyvrp::Load,
                           pyvrp::Load,
                           pyvrp::Duration,
                           pyvrp::Cost,
                           pyvrp::Duration,
                           pyvrp::Duration,
                           pyvrp::Duration,
                           pyvrp::Duration,
                           pyvrp::Duration,
                           pyvrp::Duration,
                           pyvrp::Duration,
                           pyvrp::Cost,
                           std::pair<double, double>,
                           size_t,
                           size_t,
                           size_t> const &t)
             {
                 new (&route) Route(std::get<0>(t),    // visits
                                    std::get<1>(t),    // distance
                                    std::get<2>(t),    // distance cost
                                    std::get<3>(t),    // excess distance
                                    std::get<4>(t),    // delivery
                                    std::get<5>(t),    // pickup
                                    std::get<6>(t),    // excess load
                                    std::get<7>(t),    // duration
                                    std::get<8>(t),    // duration cost
                                    std::get<9>(t),    // time warp
                                    std::get<10>(t),   // travel
                                    std::get<11>(t),   // service
                                    std::get<12>(t),   // wait
                                    std::get<13>(t),   // release
                                    std::get<14>(t),   // start time
                                    std::get<15>(t),   // slack
                                    std::get<16>(t),   // prizes
                                    std::get<17>(t),   // centroid
                                    std::get<18>(t),   // vehicle type
                                    std::get<19>(t),   // start depot
                                    std::get<20>(t));  // end depot
             })
        .def("__str__",
             [](Route const &route)
             {
                 std::stringstream stream;
                 stream << route;
                 return stream.str();
             });

    nb::class_<Solution>(m, "Solution", DOC(pyvrp, Solution))
        // Since Route implements __len__ and __getitem__, it is
        // convertible to std::vector<size_t> and thus a list of Routes
        // is a valid argument for both constructors. We want to avoid
        // using the second constructor since that would lose the
        // vehicle type associations. As pybind11 will use the first
        // matching constructor we put this one first.
        .def(nb::init<ProblemData const &, std::vector<Route> const &>(),
             nb::arg("data"),
             nb::arg("routes"))
        .def(nb::init<ProblemData const &,
                      std::vector<std::vector<size_t>> const &>(),
             nb::arg("data"),
             nb::arg("routes"))
        .def_prop_ro_static(
            "make_random",
            [](nb::object)
            {
                return nb::cpp_function(
                    [](ProblemData const &data, RandomNumberGenerator &rng)
                    { return Solution(data, rng); },
                    nb::arg("data"),
                    nb::arg("rng"),
                    DOC(pyvrp, Solution, Solution));
            })
        .def(
            "num_routes", &Solution::numRoutes, DOC(pyvrp, Solution, numRoutes))
        .def("num_clients",
             &Solution::numClients,
             DOC(pyvrp, Solution, numClients))
        .def("num_missing_clients",
             &Solution::numMissingClients,
             DOC(pyvrp, Solution, numMissingClients))
        .def("routes",
             &Solution::routes,
             nb::rv_policy::reference_internal,
             DOC(pyvrp, Solution, routes))
        .def("neighbours",
             &Solution::neighbours,
             nb::rv_policy::reference_internal,
             DOC(pyvrp, Solution, neighbours))
        .def("is_feasible",
             &Solution::isFeasible,
             DOC(pyvrp, Solution, isFeasible))
        .def("is_group_feasible",
             &Solution::isGroupFeasible,
             DOC(pyvrp, Solution, isGroupFeasible))
        .def("is_complete",
             &Solution::isComplete,
             DOC(pyvrp, Solution, isComplete))
        .def("has_excess_load",
             &Solution::hasExcessLoad,
             DOC(pyvrp, Solution, hasExcessLoad))
        .def("has_excess_distance",
             &Solution::hasExcessDistance,
             DOC(pyvrp, Solution, hasExcessDistance))
        .def("has_time_warp",
             &Solution::hasTimeWarp,
             DOC(pyvrp, Solution, hasTimeWarp))
        .def("distance", &Solution::distance, DOC(pyvrp, Solution, distance))
        .def("distance_cost",
             &Solution::distanceCost,
             DOC(pyvrp, Solution, distanceCost))
        .def("duration", &Solution::duration, DOC(pyvrp, Solution, duration))
        .def("duration_cost",
             &Solution::durationCost,
             DOC(pyvrp, Solution, durationCost))
        .def("excess_load",
             &Solution::excessLoad,
             DOC(pyvrp, Solution, excessLoad))
        .def("excess_distance",
             &Solution::excessDistance,
             DOC(pyvrp, Solution, excessDistance))
        .def("fixed_vehicle_cost",
             &Solution::fixedVehicleCost,
             DOC(pyvrp, Solution, fixedVehicleCost))
        .def("time_warp", &Solution::timeWarp, DOC(pyvrp, Solution, timeWarp))
        .def("prizes", &Solution::prizes, DOC(pyvrp, Solution, prizes))
        .def("uncollected_prizes",
             &Solution::uncollectedPrizes,
             DOC(pyvrp, Solution, uncollectedPrizes))
        .def("__copy__", [](Solution const &sol) { return Solution(sol); })
        .def(
            "__deepcopy__",
            [](Solution const &sol, nb::dict) { return Solution(sol); },
            nb::arg("memo"))
        .def("__hash__",
             [](Solution const &sol) { return std::hash<Solution>()(sol); })
        .def(nb::self == nb::self)  // this is __eq__
        .def("__getstate__",
             [](Solution const &sol)
             {
                 // Returns a tuple that completely encodes the
                 // solution's state.
                 return std::make_tuple(sol.numClients(),
                                        sol.numMissingClients(),
                                        sol.distance(),
                                        sol.distanceCost(),
                                        sol.duration(),
                                        sol.durationCost(),
                                        sol.excessDistance(),
                                        sol.excessLoad(),
                                        sol.fixedVehicleCost(),
                                        sol.prizes(),
                                        sol.uncollectedPrizes(),
                                        sol.timeWarp(),
                                        sol.isGroupFeasible(),
                                        sol.routes(),
                                        sol.neighbours());
             })
        .def("__setstate__",
             [](Solution &sol,
                std::tuple<
                    size_t,
                    size_t,
                    pyvrp::Distance,
                    pyvrp::Cost,
                    pyvrp::Duration,
                    pyvrp::Cost,
                    pyvrp::Distance,
                    pyvrp::Load,
                    pyvrp::Cost,
                    pyvrp::Cost,
                    pyvrp::Cost,
                    pyvrp::Duration,
                    bool,
                    std::vector<Route>,
                    std::vector<std::optional<std::pair<size_t, size_t>>>> const
                    &t)
             {
                 new (&sol) Solution(std::get<0>(t),    // num clients
                                     std::get<1>(t),    // num missing
                                     std::get<2>(t),    // distance
                                     std::get<3>(t),    // distance cost
                                     std::get<4>(t),    // duration
                                     std::get<5>(t),    // duration cost
                                     std::get<6>(t),    // excess distance
                                     std::get<7>(t),    // excess load
                                     std::get<8>(t),    // fixed veh cost
                                     std::get<9>(t),    // prizes
                                     std::get<10>(t),   // uncollected
                                     std::get<11>(t),   // time warp
                                     std::get<12>(t),   // is group feasible
                                     std::get<13>(t),   // routes
                                     std::get<14>(t));  // neighbours
             })
        .def("__str__",
             [](Solution const &sol)
             {
                 std::stringstream stream;
                 stream << sol;
                 return stream.str();
             });

    nb::class_<CostEvaluator>(m, "CostEvaluator", DOC(pyvrp, CostEvaluator))
        .def(nb::init<pyvrp::Cost, pyvrp::Cost, pyvrp::Cost>(),
             nb::arg("load_penalty"),
             nb::arg("tw_penalty"),
             nb::arg("dist_penalty"))
        .def("load_penalty",
             &CostEvaluator::loadPenalty,
             nb::arg("load"),
             nb::arg("capacity"),
             DOC(pyvrp, CostEvaluator, loadPenalty))
        .def("tw_penalty",
             &CostEvaluator::twPenalty,
             nb::arg("time_warp"),
             DOC(pyvrp, CostEvaluator, twPenalty))
        .def("dist_penalty",
             &CostEvaluator::distPenalty,
             nb::arg("distance"),
             nb::arg("max_distance"),
             DOC(pyvrp, CostEvaluator, twPenalty))
        .def("penalised_cost",
             &CostEvaluator::penalisedCost<Solution>,
             nb::arg("solution"),
             DOC(pyvrp, CostEvaluator, penalisedCost))
        .def("cost",
             &CostEvaluator::cost<Solution>,
             nb::arg("solution"),
             DOC(pyvrp, CostEvaluator, cost));

    nb::class_<PopulationParams>(
        m, "PopulationParams", DOC(pyvrp, PopulationParams))
        .def(nb::init<size_t, size_t, size_t, size_t, double, double>(),
             nb::arg("min_pop_size") = 25,
             nb::arg("generation_size") = 40,
             nb::arg("nb_elite") = 4,
             nb::arg("nb_close") = 5,
             nb::arg("lb_diversity") = 0.1,
             nb::arg("ub_diversity") = 0.5)
        .def(nb::self == nb::self, nb::arg("other"))  // this is __eq__
        .def_ro("min_pop_size", &PopulationParams::minPopSize)
        .def_ro("generation_size", &PopulationParams::generationSize)
        .def_prop_ro("max_pop_size",
                     &PopulationParams::maxPopSize,
                     DOC(pyvrp, PopulationParams, maxPopSize))
        .def_ro("nb_elite", &PopulationParams::nbElite)
        .def_ro("nb_close", &PopulationParams::nbClose)
        .def_ro("lb_diversity", &PopulationParams::lbDiversity)
        .def_ro("ub_diversity", &PopulationParams::ubDiversity);

    nb::class_<SubPopulation::Item>(m, "SubPopulationItem")
        .def_ro("solution",
                &SubPopulation::Item::solution,
                nb::rv_policy::reference_internal,
                R"doc(
                            Solution for this SubPopulationItem.

                            Returns
                            -------
                            Solution
                                Solution for this SubPopulationItem.
                      )doc")
        .def_ro("fitness",
                &SubPopulation::Item::fitness,
                R"doc(
                Fitness value for this SubPopulationItem.

                Returns
                -------
                float
                    Fitness value for this SubPopulationItem.

                .. warning::

                This is a cached property that is not automatically updated.
                Before accessing the property, 
                :meth:`~SubPopulation.update_fitness` should be called unless 
                the population has not changed since the last call.
            )doc")
        .def("avg_distance_closest",
             &SubPopulation::Item::avgDistanceClosest,
             R"doc(
                Determines the average distance of the solution wrapped by this
                item to a number of solutions that are most similar to it. This 
                provides a measure of the relative 'diversity' of the wrapped
                solution.

                Returns
                -------
                float
                    The average distance/diversity of the wrapped solution.
             )doc");

    nb::class_<SubPopulation>(m, "SubPopulation", DOC(pyvrp, SubPopulation))
        .def(nb::init<pyvrp::diversity::DiversityMeasure,
                      PopulationParams const &>(),
             nb::arg("diversity_op"),
             nb::arg("params"),
             nb::keep_alive<1, 3>())  // keep params alive
        .def("add",
             &SubPopulation::add,
             nb::arg("solution"),
             nb::arg("cost_evaluator"),
             DOC(pyvrp, SubPopulation, add))
        .def("__len__", &SubPopulation::size)
        .def(
            "__getitem__",
            [](SubPopulation const &subPop, int idx)
            {
                // int so we also support negative offsets from the
                // end.
                idx = idx < 0 ? subPop.size() + idx : idx;
                if (idx < 0 || static_cast<size_t>(idx) >= subPop.size())
                    throw nb::index_error();
                return subPop[idx];
            },
            nb::arg("idx"),
            nb::rv_policy::reference_internal)
        .def(
            "__iter__",
            [](SubPopulation const &subPop)
            {
                return nb::make_iterator(nb::type<SubPopulation>(),
                                         "iterator",
                                         subPop.cbegin(),
                                         subPop.cend());
            },
            nb::rv_policy::reference_internal)
        .def("purge",
             &SubPopulation::purge,
             nb::arg("cost_evaluator"),
             DOC(pyvrp, SubPopulation, purge))
        .def("update_fitness",
             &SubPopulation::updateFitness,
             nb::arg("cost_evaluator"),
             DOC(pyvrp, SubPopulation, updateFitness));

    nb::class_<DistanceSegment>(
        m, "DistanceSegment", DOC(pyvrp, DistanceSegment))
        .def(nb::init<size_t, size_t, pyvrp::Distance>(),
             nb::arg("idx_first"),
             nb::arg("idx_last"),
             nb::arg("distance"))
        .def("distance",
             &DistanceSegment::distance,
             DOC(pyvrp, DistanceSegment, distance))
        .def_static("merge",
                    &DistanceSegment::merge<>,
                    nb::arg("distance_matrix"),
                    nb::arg("first"),
                    nb::arg("second"))
        .def_static("merge",
                    &DistanceSegment::merge<DistanceSegment const &>,
                    nb::arg("distance_matrix"),
                    nb::arg("first"),
                    nb::arg("second"),
                    nb::arg("third"));

    nb::class_<LoadSegment>(m, "LoadSegment", DOC(pyvrp, LoadSegment))
        .def(nb::init<pyvrp::Load, pyvrp::Load, pyvrp::Load>(),
             nb::arg("delivery"),
             nb::arg("pickup"),
             nb::arg("load"))
        .def("delivery",
             &LoadSegment::delivery,
             DOC(pyvrp, LoadSegment, delivery))
        .def("pickup", &LoadSegment::pickup, DOC(pyvrp, LoadSegment, pickup))
        .def("load", &LoadSegment::load, DOC(pyvrp, LoadSegment, load))
        .def_static(
            "merge", &LoadSegment::merge<>, nb::arg("first"), nb::arg("second"))
        .def_static("merge",
                    &LoadSegment::merge<LoadSegment const &>,
                    nb::arg("first"),
                    nb::arg("second"),
                    nb::arg("third"));

    nb::class_<DurationSegment>(
        m, "DurationSegment", DOC(pyvrp, DurationSegment))
        .def(nb::init<size_t,
                      size_t,
                      pyvrp::Duration,
                      pyvrp::Duration,
                      pyvrp::Duration,
                      pyvrp::Duration,
                      pyvrp::Duration>(),
             nb::arg("idx_first"),
             nb::arg("idx_last"),
             nb::arg("duration"),
             nb::arg("time_warp"),
             nb::arg("tw_early"),
             nb::arg("tw_late"),
             nb::arg("release_time"))
        .def("duration",
             &DurationSegment::duration,
             DOC(pyvrp, DurationSegment, duration))
        .def("tw_early",
             &DurationSegment::twEarly,
             DOC(pyvrp, DurationSegment, twEarly))
        .def("tw_late",
             &DurationSegment::twLate,
             DOC(pyvrp, DurationSegment, twLate))
        .def("time_warp",
             &DurationSegment::timeWarp,
             nb::arg("max_duration")
             = std::numeric_limits<pyvrp::Duration>::max(),
             DOC(pyvrp, DurationSegment, timeWarp))
        .def_static("merge",
                    &DurationSegment::merge<>,
                    nb::arg("duration_matrix"),
                    nb::arg("first"),
                    nb::arg("second"))
        .def_static("merge",
                    &DurationSegment::merge<DurationSegment const &>,
                    nb::arg("duration_matrix"),
                    nb::arg("first"),
                    nb::arg("second"),
                    nb::arg("third"));

    nb::class_<RandomNumberGenerator>(
        m, "RandomNumberGenerator", DOC(pyvrp, RandomNumberGenerator))
        .def(nb::init<uint32_t>(), nb::arg("seed"))
        .def(nb::init<std::array<uint32_t, 4>>(), nb::arg("state"))
        .def("min", &RandomNumberGenerator::min)
        .def("max", &RandomNumberGenerator::max)
        .def("__call__", &RandomNumberGenerator::operator())
        .def("rand", &RandomNumberGenerator::rand)
        .def("randint", &RandomNumberGenerator::randint<int>, nb::arg("high"))
        .def("state", &RandomNumberGenerator::state);
}
