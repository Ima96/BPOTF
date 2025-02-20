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
#include "SDemData/SDemData.h"
#include "docstrings.h"

namespace py = pybind11;

// Bindings for the BPOTF module
// TODO: Add docstrings to the methods and attributes.
PYBIND11_MODULE(BPOTF, mBPOTF) {

   // Initialize dependencies. Done this way to ensure that python .so and .dll and dependencies are available when
   // the imports are called.
   initialize_BPOTF_dependencies();
   
   // Set module's version (set at compile time)
   mBPOTF.attr("__version__") = BPOTF_VERSION;

   // Export class
   auto py_BPOTF = py::class_<OBPOTF>(mBPOTF, "OBPOTF");

   // Export enumeration typedef for different error sources supported
   py::enum_<ENoiseType_t>(mBPOTF, "NoiseType")
      .value("E_CC", ENoiseType_t::E_CC, 
         R"pbdoc(
            Code Capacity kind of noise.
         )pbdoc")
      .value("E_PHEN", ENoiseType_t::E_PHEN, 
         R"pbdoc(
            Phenomenological kind of noise. (No support yet)
         )pbdoc")
      .value("E_CLN", ENoiseType_t::E_CLN, 
         R"pbdoc(
            Circuit-Level-Noise type. Select this to build object from DEMs.
         )pbdoc")
      .export_values();

   // External struct to receive Dem Data from user
   auto py_DemData = py::class_<SDemData_t>(mBPOTF, "DemData");
   py_DemData.def(py::init<>());
   py_DemData.def_readwrite("obs_matrix", &SDemData_t::po_obs_csr_mat);
   py_DemData.def_readwrite("phen_obs_matrix", &SDemData_t::po_phen_obs_csr);
   py_DemData.def_readwrite("phen_check_matrix", &SDemData_t::po_phen_pcm_csc);
   py_DemData.def_readwrite("transfer_matrix", &SDemData_t::po_transfer_csr_mat);
   py_DemData.def_readwrite("priors", &SDemData_t::af64_priors);

   // Export class's public methods.
   py_BPOTF
      .def(py::init<py::object const &, float const &, ENoiseType_t const, py::object const &, SDemData_t const *>(),
            py::arg("pcm"),   // Parity-check matrix parameter
            py::arg("p"),     // Physical error probability
            py::arg("noise_type") = ENoiseType_t::E_CC,   // Noise model type. Default: E_CC
            py::arg("po_ext_bp_iters") = py::none(),  // External BP max iterations (array of ints)
            py::arg("ps_ext_dem_data") = py::none(),  // External DEM data
            docstr_bpotf_constructor
         );
   py_BPOTF.def("decode", &OBPOTF::decode);
   py_BPOTF.def("has_converged", &OBPOTF::has_converged);
      
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