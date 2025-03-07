/***********************************************************************************************************************
 * @file    SDemData.h
 * @author  Imanol Etxezarreta (ietxezarretam@gmail.com)
 * 
 * @brief   Structure to expose to python for the user to fill and pass to the OBPOTF object in case a BPBPOTF decode
 *          scheme is wanted.
 * 
 * @version 0.1
 * @date    15/01/2025
 * 
 * @copyright Copyright (c) 2025
 * 
 **********************************************************************************************************************/
#ifndef SDEMDATA_H_
#define SDEMDATA_H_

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

typedef struct SDemData
{  
   //! Transfer matrix in case a DEM is provided.
   py::object po_transfer_csr_mat = py::none();

   //!< Observables matrix in case DEM is provided.
   py::object po_obs_csr_mat = py::none();

   //!< Phenomenological CSC matrix of the pcm.
   py::object po_phen_pcm_csc = py::none();

   //!< Phenomenological CSR matrix of the observables.
   py::object po_phen_obs_csr = py::none();

   //!< Array of prior probabilities
   py::object af64_priors = py::none();

} SDemData_t;

#endif // SDEMDATA_H_