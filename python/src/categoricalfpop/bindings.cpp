#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace ei = Eigen;
namespace py = pybind11;
using namespace pybind11::literals; // to have named arguments

ei::ArrayXi getPeltCP(const ei::ArrayXXi &signal, double penalty);

ei::ArrayXi getAFPCP(const ei::ArrayXXi &signal, double penalty);

ei::ArrayXi getPeltCP_py(const py::EigenDRef<ei::ArrayXXi> signal, double penalty)
{
      return getPeltCP(signal, penalty);
}

ei::ArrayXi getAFPCP_py(const py::EigenDRef<ei::ArrayXXi> signal, double penalty)
{
      return getAFPCP(signal, penalty);
}

PYBIND11_MODULE(categoricalfpop, m)
{
      m.doc() = "Change-point detection using functional pruning."; // optional
      m.def("get_peltcp", &getPeltCP_py,
            "Best change-points (with PELT).",
            "signal"_a.noconvert(), "penalty"_a);
      m.def("get_afpcp", &getAFPCP_py,
            "Best change-points (with AFP).",
            "signal"_a.noconvert(), "penalty"_a);
}
