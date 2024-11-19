/***********************************************************************************************************************
 * @file    py11_iface.cpp
 * @author  Imanol Etxezarreta (ietxezarretam@gmail.com)
 * 
 * @brief   This file is an interface with Pybind11 to expose the BPOTF class and its methods to a python module.
 * 
 * @version 0.1
 * @date    21/08/2024
 * 
 * @copyright Copyright (c) 2024
 * 
 **********************************************************************************************************************/

// Pybind11 header libraries
#include <pybind11/pybind11.h>
// #include <pybind11/numpy.h>

// Custom headers
#include "BPOTF/OBPOTF.h"
#include "docstrings.h"

namespace py = pybind11;

// Bindings for the BPOTF module
// TODO: Add docstrings to the methods and attributes.
PYBIND11_MODULE(BPOTF, BPOTF) {
   
   // Set module's version (set at compile time)
   BPOTF.attr("__version__") = BPOTF_VERSION;

   // Export class
   auto py_BPOTF = py::class_<OBPOTF>(BPOTF, "OBPOTF");

   // Export enumeration typedef for different error sources supported
   py::enum_<OBPOTF::ENoiseType_t>(py_BPOTF, "NoiseType")
      .value("E_CC", OBPOTF::ENoiseType_t::E_CC, 
         R"pbdoc(
            Code Capacity kind of noise.
         )pbdoc")
      .value("E_PHEN", OBPOTF::ENoiseType_t::E_PHEN, 
         R"pbdoc(
            Phenomenological kind of noise. (No support yet)
         )pbdoc")
      .value("E_CLN", OBPOTF::ENoiseType_t::E_CLN, 
         R"pbdoc(
            Circuit-Level-Noise type. Select this to build object from DEMs.
         )pbdoc")
      .export_values();

   // Export class's public methods.
   py_BPOTF
      .def(py::init<py::object const &, float const &, OBPOTF::ENoiseType_t const, py::object const * const>(),
            py::arg("pcm"),   // Parity-check matrix parameter
            py::arg("p"),     // Physical error probability
            py::arg("noise_type") = OBPOTF::ENoiseType_t::E_CC,   // Noise model type. Default: E_CC
            py::arg("transfer_mat") = nullptr,  // Transfer matrix
            docstr_bpotf_constructor
         );
   py_BPOTF.def("otf_uf_probs", py::overload_cast<py::array_t<double, C_FMT> const &>(&OBPOTF::otf_uf_probs));
   py_BPOTF.def("decode", &OBPOTF::decode);
      
#if defined(DEBUG_OBPOTF)
   py_BPOTF.def("get_pcm", &OBPOTF::getPcm);
   py_BPOTF.def("get_phen_pcm", &OBPOTF::getPhenPcm);
   py_BPOTF.def("get_obs", &OBPOTF::getObs);
   py_BPOTF.def("get_transf", &OBPOTF::getTransfMat);
   py_BPOTF.def("get_priors", &OBPOTF::getPriors);
   py_BPOTF.def("get_cols", &OBPOTF::get_cols);
   py_BPOTF.def("get_rows", &OBPOTF::get_rows);
   py_BPOTF.def("get_obs_cols", &OBPOTF::get_cols_obs);
   py_BPOTF.def("get_obs_rows", &OBPOTF::get_rows_obs);
   py_BPOTF.def("get_phen_pcm_cols", &OBPOTF::get_cols_phen_pcm);
   py_BPOTF.def("get_phen_pcm_rows", &OBPOTF::get_rows_phen_pcm);
   py_BPOTF.def("get_transf_cols", &OBPOTF::get_cols_transf);
   py_BPOTF.def("get_transf_rows", &OBPOTF::get_rows_transf);
#endif
   py_BPOTF.def("print_object", &OBPOTF::print_object);

}