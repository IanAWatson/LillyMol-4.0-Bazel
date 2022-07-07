#include <pybind11/pybind11.h>

#include "Molecule_Tools/unique_molecules_api.h"

namespace py = pybind11;
int add(int i, int j) {
    return i + j;
}

PYBIND11_MODULE(qunique_molecules, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &add, "A function that adds two numbers");
}

class FooBar {
  private:
    std::string _value;
  public:
    FooBar() {
      _value = "abc";
    }
    void set_value(const std::string& s) {
      _value = s;
    }
    const std::string& get_value() const {
      return _value;
    }
};

#ifdef THIS_WORKS
PYBIND11_MODULE(unique_molecules, m) {
  py::class_<FooBar>(m, "foobar")
    .def(py::init<>())
    .def("get_value", &FooBar::get_value)
    .def("set_value", &FooBar::set_value);
}
#endif

using unique_molecules::UniqueMoleculesImplementation;
PYBIND11_MODULE(unique_molecules, m) {
  py::class_<UniqueMoleculesImplementation>(m, "unique_molecules")
    .def(py::init<>())
    .def("set_reduce_to_largest_fragment", &UniqueMoleculesImplementation::set_reduce_to_largest_fragment)
    .def("set_exclude_chiral_info", &UniqueMoleculesImplementation::set_exclude_chiral_info)
    .def("set_exclude_cis_trans_bonding_info", &UniqueMoleculesImplementation::set_exclude_cis_trans_bonding_info)
    .def("set_ignore_isotopes", &UniqueMoleculesImplementation::set_ignore_isotopes)
    .def("unique_molecules", &UniqueMoleculesImplementation::unique_molecules)
    .def("activate_standardization", &UniqueMoleculesImplementation::ActivateChemicalStandardization)
    .def("is_unique", &UniqueMoleculesImplementation::SmilesIsUnique);
    
}
