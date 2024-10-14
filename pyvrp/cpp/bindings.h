#include "Matrix.h"
#include "Measure.h"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <cassert>
#include <type_traits>

namespace nanobind::detail
{
// This is not a fully general type caster for Matrix. Instead, it assumes
// we're casting elements that are, or have the same size as, pyvrp::Value,
// which is the case for e.g. the Measure types.
template <typename T> struct type_caster<pyvrp::Matrix<T>>
{
    static_assert(sizeof(T) == sizeof(pyvrp::Value)
                  && std::is_constructible_v<T, pyvrp::Value>);

    using Value = pyvrp::Matrix<T>;

    static constexpr auto Name = "numpy.ndarray[int]";

    bool from_python(handle src, uint8_t flags, cleanup_list *cleanup) noexcept
    {
        if (!convert && !array_t<pyvrp::Value>::check_(src))
            return false;

        auto const style = array::c_style | array::forcecast;
        auto const buf = array_t<pyvrp::Value, style>::ensure(src);

        if (!buf || buf.ndim() != 2)
            return false;

        if (buf.size() == 0)  // then the default constructed object is already
            return true;      // OK, and we have nothing to do.

        std::vector<T> data = {buf.data(), buf.data() + buf.size()};
        value = pyvrp::Matrix<T>(data, buf.shape(0), buf.shape(1));

        return true;
    }

    static handle cast(pyvrp::Matrix<T> const &src,  // C++ -> Python
                       [[maybe_unused]] return_value_policy policy,
                       handle parent)
    {
        auto constexpr elemSize = sizeof(pyvrp::Value);

        array_t<pyvrp::Value> array
            = {{src.numRows(), src.numCols()},                      // shape
               {elemSize * src.numCols(), elemSize},                // strides
               reinterpret_cast<pyvrp::Value const *>(src.data()),  // data
               parent};                                             // base

        // This is not pretty, but it makes the matrix non-writeable on the
        // Python side. That's needed because src is const, and we should
        // preserve that to avoid issues.
        detail::array_proxy(array.ptr())->flags
            &= ~detail::npy_api::NPY_ARRAY_WRITEABLE_;

        return array.release();
    }
};

// On the C++ side we have strong types for different measure values (for
// example distance, duration, etc.), but on the Python side those things
// are just ints. This type caster converts between the two.
template <pyvrp::MeasureType T> struct type_caster<pyvrp::Measure<T>>
{
    using Value = pyvrp::Measure<T>;

    static constexpr auto Name = "int";

    bool from_python(handle src, uint8_t flags, cleanup_list *cleanup) noexcept
    {
        static_assert(sizeof(long long) >= sizeof(pyvrp::Value));

        if (!convert && !PyLong_Check(src.ptr()))  // only int when conversion
            return false;                          // is not allowed.

        PyObject *tmp = PyNumber_Long(src.ptr());  // any argument for which
        if (!tmp)                                  // Python's int() succeeds.
            return false;

        auto const raw = PyLong_AsLongLong(tmp);
        Py_DECREF(tmp);

        // See https://docs.python.org/3/c-api/long.html#c.PyLong_AsLongLong:
        // -1 is returned on overflow, and OverflowError is set.
        if (raw == -1 && PyErr_Occurred())
        {
            PyErr_Clear();
            return false;
        }

        value = Value(raw);
        return true;
    }

    static handle from_cpp(Value const *value,
                           rv_policy policy,
                           cleanup_list *cleanup) noexcept
    {
        assert(value);
        return PyLong_FromLongLong(value->get());
    }
};
}  // namespace nanobind::detail
