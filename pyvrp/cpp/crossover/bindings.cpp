#include "crossover_docs.h"
#include "ordered_crossover.h"
#include "selective_route_exchange.h"

#include <nanobind/nanobind.h>

namespace nb = nanobind;

NB_MODULE(_crossover, m)
{
    m.def("ordered_crossover",
          &pyvrp::crossover::orderedCrossover,
          nb::arg("parents"),
          nb::arg("data"),
          nb::arg("indices"),
          DOC(pyvrp, crossover, orderedCrossover));

    m.def("selective_route_exchange",
          &pyvrp::crossover::selectiveRouteExchange,
          nb::arg("parents"),
          nb::arg("data"),
          nb::arg("cost_evaluator"),
          nb::arg("start_indices"),
          nb::arg("num_moved_routes"),
          DOC(pyvrp, crossover, selectiveRouteExchange));
}
