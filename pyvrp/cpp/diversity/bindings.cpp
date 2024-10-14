#include "diversity.h"
#include "diversity_docs.h"

#include <nanobind/nanobind.h>

namespace nb = nanobind;

NB_MODULE(_diversity, m)
{
    m.def("broken_pairs_distance",
          &pyvrp::diversity::brokenPairsDistance,
          nb::arg("first"),
          nb::arg("second"),
          DOC(pyvrp, diversity, brokenPairsDistance));
}
