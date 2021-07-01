#include <pybind11/pybind11.h>
#include <pybind11/stl.h>     

#include "plugins/Interstitial_mesh_plugins/interstitial_mesh.hpp"
#include "plugins/Interstitial_mesh_plugins/interstitial_mesh.cpp" 



PYBIND11_MODULE(utilsplugin, m) {

    m.def("make_orbit", &make_orbit);
}
