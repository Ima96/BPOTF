/***********************************************************************************************************************
 * @file    OBPOTF.h
 * @author  Imanol Etxezarreta (ietxezarretam@gmail.com)
 * 
 * @brief   Object of a BPOTF decoder, with associated methods. BPOTF decoder uses Kruskal's algorithm and 
 *          Disjoint-Set Advanced Data Structure to offer a quick Low-Density Parity Check syndrome decodification. 
 *          (Revisar)
 * 
 * @version 0.1
 * @date    17/05/2024
 * 
 * @copyright Copyright (c) 2024
 * 
 **********************************************************************************************************************/
#ifndef OBPOTF_H_
#define OBPOTF_H_

// std includes
#include <vector>
#include <inttypes.h>

// Pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

// Custom includes
#include "CSC/OCSC.h"
#include "../ldpc_v2_src/bp.hpp"

namespace py = pybind11;

#define C_FMT py::array::c_style
#define F_FMT py::array::f_style

class OBPOTF
{
   /******************************************************************************************************************* 
    * PRIVATE MEMBER VARIABLE DECLARATION 
    *******************************************************************************************************************/
   private:

   //! Error probabilities.
   float const m_p;
   //! Parity check matrix number of rows.
   uint64_t m_u64_pcm_rows;
   //! Parity check matrix number of columns.
   uint64_t m_u64_pcm_cols;

   // TODO: Rename to specify explicitly that it is the pcm matrix
   //! Matrix in Compressed-Sparse-Column format.
   OCSC * m_po_csc_mat = nullptr;

   typedef struct SDemData
   {
      //! Transfer matrix in case a DEM is provided.
      OCSC * po_transfer_csc_mat = nullptr;

      //!< Observables CSC matrix in case DEM is provided.
      OCSC * po_obs_csc_mat = nullptr;

      //!< Phenomenological CSC matrix of the pcm.
      OCSC * po_phen_pcm_csc = nullptr;

      //!< Phenomenological CSC matrix of the observables.
      OCSC * po_phen_obs_csc = nullptr;

     //!< Array of prior probabilities
     std::vector<double> af64_priors; 

   } SDemData_t;

   SDemData_t * m_ps_dem_data = nullptr;

   //! Pointer to the pcm in format for BpDecoder object
   ldpc::bp::BpSparse * m_po_bpsparse_pcm = nullptr;
   //! Pointer to the phenomenological pcm in format for BpDecoder object
   ldpc::bp::BpSparse * m_po_bpsparse_phen = nullptr;
   //! Pointer to BP decoder object to use it against the pcm.
   ldpc::bp::BpDecoder * m_po_pcm_bp = nullptr;
   //! Pointer to BP decoder object to use it against the phenomenological pcm in case there is one.
   ldpc::bp::BpDecoder * m_po_phen_bp = nullptr;
   //! Pointer to BP decoder in case the OTF is performed and use BP against its result. 
   ldpc::bp::BpDecoder * m_po_otf_bp = nullptr;

   //! Array that holds indexes from 0 to m_u64_pcm_cols-1 to be sorted.
   std::vector<uint64_t> m_au64_index_array;

   //! Callback to the decoding method selected when constructing the object.
   py::array_t<uint8_t> (OBPOTF::*m_pf_decoding_func)(py::array_t<uint8_t, C_FMT> const &);

   /********************************************************************************************************************
    * PUBLIC MEMBER VARIABLE DECLARATION 
    *******************************************************************************************************************/
   public:

   /********************************************************************************************************************
    * @typedef ENoiseType_t
    * @brief   This typedef holds the different noise models that can be passed to the decoder. Depending on the type
    *          passed, the decoder could also use (by input argument or trying to obtaining it) a transference matrix
    *          to simplify the decoding procedure.
    *******************************************************************************************************************/
   typedef enum ENoiseType
   {
      E_CC     = 0,  //!< Code Capacity (Default mode).
      E_PHEN   = 1,  //!< Phenomenological.
      E_CLN    = 2   //!< Circuit-level noise.
   } ENoiseType_t;

   /********************************************************************************************************************
    * PRIVATE CLASS METHOD DECLARATION
    *******************************************************************************************************************/
   private:

   /********************************************************************************************************************
    * @brief Sub-routine that is called from the object constructor if it is called with a numpy array. It initialized 
    *        the object members from input parameters and executes necessary pre-processings.
    * 
    * @param pcm[in] Parity-check matrix from which to initialize the members.
    *******************************************************************************************************************/
   void OBPOTF_init_from_numpy(py::array_t<uint8_t, F_FMT> const & pcm);
   
   /********************************************************************************************************************
    * @brief Sub-routine that is called from the object constructor if it is called with a scipy_csc object. In this 
    *        case, the object is converted to a pyarray and the the OBPOTF_init_from_numpy is called with it. 
    * 
    * @param pcm[in] Parity-check matrix from which to initialize the members.
    *******************************************************************************************************************/
   void OBPOTF_init_from_scipy_csc(py::object const & pcm);

   /********************************************************************************************************************
    * @brief Sub-routine that is called from the object constructor if this is fed with a stim dem object. In this 
    *        case, the PCM is extracted from the DEM and a transfer matrix is passed as argument to simplify the 
    *        decoding process. In case there is no transfer matrix as input, it will try to generate one if possible.
    * 
    * @param po_dem[in]          Input DEM object. PCM is extracted from it.
    * @param po_transfer_mat[in] Input transfer matrix. If it is null, try to generate a new one.
    *******************************************************************************************************************/
   void OBPOTF_init_from_dem(py::object const & po_dem, py::object const * const po_transfer_mat);

   /********************************************************************************************************************
    * @brief This routine performs the OTF algorithm using the clasical Unified-Find method.
    * 
    * @param llrs[in]   The llrs is a vector containing the probabilities of the PCM columns.
    * @return std::vector<uint64_t> The return value is a vector containing the recovered error.
    *******************************************************************************************************************/
   std::vector<uint64_t> otf_classical_uf(std::vector<double> const & llrs);

   /********************************************************************************************************************
    * @brief This routine performs the OTF algorithm.
    * 
    * @param llrs[in]   The llrs is a vector containing the probabilities of the PCM columns.
    * @return std::vector<uint64_t> The return value is a vector containing the recovered error.
    *******************************************************************************************************************/
   std::vector<uint64_t> otf_uf(std::vector<double> const & llrs);

   /********************************************************************************************************************
    * @brief This routine returns a vector of the sorted indexes based on the probabilities of the llrs. It copies the 
    *        initial vector with the unsorted indexes from the member variable m_au64_index_array.
    * 
    * @param llrs[in]   The llrs is a vector containing the probabilities of the PCM columns.
    * @return std::vector<uint64_t> The return variable is a vector with the sorted indexes according the llrs.
    *******************************************************************************************************************/
   std::vector<uint64_t> sort_indexes(py::array_t<double> const & llrs);

   /********************************************************************************************************************
    * @brief This routine returns a sorted vector of pointers to the sorted indexes stored in the member variable
    *        m_au64_index_array, based on the probabilities in the llrs.
    * 
    * @param llrs[in]   The llrs is a vector containing the probabilities of the PCM columns.
    * @return std::vector<uint64_t *> The return variable is a sorted vector of pointers to the indexes based on the
    *                                 llrs.
    *******************************************************************************************************************/
   std::vector<uint64_t *> sort_indexes_nc(std::vector<double> const & llrs);

   /********************************************************************************************************************
    * @brief This routine returns a sorted vector of pointers to the sorted indexes stored in the member variable
    *        m_au64_index_array, based on the probabilities in the llrs.
    * 
    * @param llrs[in]   The llrs is a span containing the probabilities of the PCM columns.
    * @return std::vector<uint64_t *> The return variable is a sorted vector of pointers to the indexes based on the
    *                                 llrs.
    *******************************************************************************************************************/
   std::vector<uint64_t *> sort_indexes_nc(std::span<double> const & llrs);

   std::vector<double> propagate(std::vector<double> const & vec_f_probs);

   /********************************************************************************************************************
    * @brief This routine executes a generic decode procedure, which is done for surface-codes. It is registered as a 
    *        callback in the member variable m_pf_decoding_func when the object is created if the enumeration type is
    *        set to E_GENERIC.
    * 
    * @param syndrome[in]  A python array in c-style format that indicates the syndrome from which recover the error.
    * @return py::array_t<uint8_t> Output python array with the resulting recovered error.
    *******************************************************************************************************************/
   py::array_t<uint8_t> generic_decode(py::array_t<uint8_t, C_FMT> const & syndrome);

   /********************************************************************************************************************
    * @brief This routine executes the decoding process for circuit-level noise type of errors. It is registered as a 
    *        callback in the member variable m_pf_decoding_func when the object is created if the enumeration type is
    *        set to E_CLN.
    * 
    * @param syndrome[in]  A python array in c-style format that indicates the syndrome from which recover the error.
    * @return py::array_t<uint8_t> Output python array with the resulting recovered error.
    *******************************************************************************************************************/
   py::array_t<uint8_t> cln_decode(py::array_t<uint8_t, C_FMT> const & syndrome);

   /********************************************************************************************************************
    * PUBLIC CLASS METHOD DECLARATION
    *******************************************************************************************************************/
   public:

   /********************************************************************************************************************
    * @brief Construct a new OBPOTF object from the input values. It calls other sub-routines depending on the python 
    *        object that is passed as a parameter.
    * 
    * @param pcm[in]          Parity check matrix. It is passed as a py::object for speed and avoid copying the matrix.
    * @param p[in]            Phisical error to initialize the bp_decoder.
    * @param noise_type[in]   Type of the noise source.
    * @param transfer_mat[in] Transference matrix to try to simplify the decoding process.
    *******************************************************************************************************************/
   OBPOTF(py::object const & pcm, float const & p, ENoiseType_t const noise_type,
            py::object const * const transfer_mat);

   /********************************************************************************************************************
    * @brief Delete default constructor, to avoid empty objects.
    *******************************************************************************************************************/
   OBPOTF() = delete;

   /********************************************************************************************************************
    * @brief Destroy the OBPOTF object to free allocated memory.
    *******************************************************************************************************************/
   ~OBPOTF();

   /********************************************************************************************************************
    * @brief This routine performs the OTF algorithm.
    * 
    * (It is public right now because Ton needed only this part of the algorithm for some tests)
    * 
    * @param llrs[in]   The llrs is a vector containing the probabilities of the PCM columns.
    * @return py::array_t<uint64_t> The return value is a py::array_t containing the recovered error.
    *******************************************************************************************************************/
   py::array_t<uint64_t> otf_uf(py::array_t<double, C_FMT> const & llrs);

   /********************************************************************************************************************
    * @brief This is the main decoding routine. This routine calls the registered callback function in 
    *        m_pf_decoding_func with the syndrome to recover the error.
    * 
    * @param syndrome[in]  A python array in c-style format that indicates the syndrome from which to recover the error.
    * @return py::array_t<uint8_t> Returned value is a python array with the recovered error.
    *******************************************************************************************************************/
   py::array_t<uint8_t> decode(py::array_t<uint8_t, C_FMT> const & syndrome);

   /********************************************************************************************************************
    * @brief Prints the object's member. Developing purposes and testing.
    *******************************************************************************************************************/
   void print_object(void);

   py::array_t<uint8_t> getPcm(void);
   py::array_t<uint8_t> getPhenPcm(void);
   py::array_t<uint8_t> getObs(void);
   py::array_t<double> getPriors(void);

   inline uint64_t get_cols(void) { return m_u64_pcm_cols; }
   inline uint64_t get_rows(void) { return m_u64_pcm_rows; }
   inline uint64_t get_cols_obs(void) { return m_ps_dem_data->po_obs_csc_mat->get_col_num(); }
   inline uint64_t get_rows_obs(void) { return m_ps_dem_data->po_obs_csc_mat->get_row_num(); }
   inline uint64_t get_cols_phen_pcm(void) { return m_ps_dem_data->po_phen_pcm_csc->get_col_num(); }
   inline uint64_t get_rows_phen_pcm(void) { return m_ps_dem_data->po_phen_pcm_csc->get_row_num(); }
};

#endif // OBPOTF_H_