#include "greedy_repair.h"
#include "nearest_route_insert.h"
#include "repair_docs.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

NB_MODULE(_repair, m)
{
    m.def("greedy_repair",
          &pyvrp::repair::greedyRepair,
          nb::arg("routes"),
          nb::arg("unplanned"),
          nb::arg("data"),
          nb::arg("cost_evaluator"),
          DOC(pyvrp, repair, greedyRepair));

    m.def("nearest_route_insert",
          &pyvrp::repair::nearestRouteInsert,
          nb::arg("routes"),
          nb::arg("unplanned"),
          nb::arg("data"),
          nb::arg("cost_evaluator"),
          DOC(pyvrp, repair, nearestRouteInsert));
}
