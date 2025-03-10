/***********************************************************************************************************************
 * @file    OBPOTF.cpp
 * @author  Imanol Etxezarreta (ietxezarretam@gmail.com)
 * 
 * @brief   Implementation of the object and its methods for OBPOTF.
 * 
 * @version 0.1
 * @date    17/05/2024
 * 
 * @copyright Copyright (c) 2024
 * 
 **********************************************************************************************************************/

// std includes
#include <numeric>
#include <span>
#include <iostream>

#if defined(TIMING_OBPOTF)
// Timing tests
#include <chrono>
#include <iomanip>
#endif

// Windows build includes
#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

// Custom includes
#include "OBPOTF.h"
#include "SDemData/SDemData.h"
#include "DisjointSet/DisjointSet.h"
#include "CSC/OCSC.h"

using namespace py::literals;

/***********************************************************************************************************************
 * DEFINES
 **********************************************************************************************************************/
#define CONV_2_U8_NDARR(X) X.attr("toarray")("F").attr("astype")("uint8")

#if defined(TIMING_OBPOTF)
std::chrono::_V2::system_clock::time_point vf_start_time;
std::chrono::_V2::system_clock::time_point vf_stop_time;
std::chrono::nanoseconds::rep vf_duration_ns;
#define  START_CHRONO   vf_start_time = std::chrono::high_resolution_clock::now();
#define  STOP_CHRONO(str)  {\
   vf_stop_time = std::chrono::high_resolution_clock::now(); \
   vf_duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(vf_stop_time-vf_start_time).count(); \
   std::cout << str << vf_duration_ns*1e-9 << std::setprecision(9) << " s" << std::endl; \
   }
#else
#define START_CHRONO
#define STOP_CHRONO(str)
#endif

/***********************************************************************************************************************
 * FILE GLOBAL VARIABLES
 **********************************************************************************************************************/
// Global variable that holds scipy.sparse.csc_matrix type 
py::object vf_scipy_csc_type;
// Global constant variable that holds the initial llr value for probability 1e-14
const double vfc_initial_llr_value = std::log((1.0 - 1e-14) / 1e-14);

/***********************************************************************************************************************
 * Helper functions
 **********************************************************************************************************************/
void initialize_BPOTF_dependencies(void)
{
   // Import scipy.sparse.csc_matrix from scipy
   vf_scipy_csc_type = py::module_::import("scipy.sparse").attr("csc_matrix");
}

/***********************************************************************************************************************
 * @brief helper function to avoid making a copy when returning a py::array_t
 * 
 *    author: https://github.com/YannickJadoul
 *    source: https://github.com/pybind/pybind11/issues/1042#issuecomment-642215028
 * 
 * @tparam Sequence::value_type 
 * @param seq[in]    input sequence to return as py::array.
 * @return py::array_t<typename Sequence::value_type> 
 **********************************************************************************************************************/
template <typename Sequence>
inline static py::array_t<typename Sequence::value_type> as_pyarray(Sequence &&seq) {
   auto size = seq.size();
   auto data = seq.data();
   std::unique_ptr<Sequence> seq_ptr = std::make_unique<Sequence>(std::move(seq));
   auto capsule = py::capsule(seq_ptr.get(), 
                              [](void *p) 
                              {
                                 std::unique_ptr<Sequence>(reinterpret_cast<Sequence *>(p));
                              }
                              );
   seq_ptr.release();

   return py::array(size, data, capsule);
}

/***********************************************************************************************************************
 * @brief Returns span<T> from py:array_T<T>. Efficient as zero-copy. Only works with py::array_t with 2 dimensions.
 * 
 * @tparam T Type
 * @tparam I Array format enumeration.
 * @param passthrough[in] Numpy array to get as span.
 * @return std::span<T> clean and safe reference to contents of Numpy array.
 **********************************************************************************************************************/
template<typename T, int I>
inline static std::span<T> toSpan2D(py::array_t<T, I> const & passthrough)
{
   py::buffer_info passthroughBuf = passthrough.request();
   if (passthroughBuf.ndim != 2) {
      throw std::runtime_error("Error. Number of dimensions must be two");
   }
   uint64_t length = passthroughBuf.shape[0] * passthroughBuf.shape[1];
   T* passthroughPtr = static_cast<T*>(passthroughBuf.ptr);
   std::span<T> passthroughSpan(passthroughPtr, length);
   return passthroughSpan;
}

/***********************************************************************************************************************
 * @brief Returns span<T> from py:array_T<T>. Efficient as zero-copy. Only works with py::array_t with 1 dimension.
 * 
 * @tparam T Type
 * @param passthrough[in] Numpy array to get as span.
 * @return std::span<T> clean and safe reference to contents of Numpy array.
 **********************************************************************************************************************/
template<typename T>
inline static std::span<T> toSpan1D(py::array_t<T, C_FMT> const & passthrough)
{
   py::buffer_info passthroughBuf = passthrough.request();
   if (passthroughBuf.ndim != 1) {
      throw std::runtime_error("Error. Number of dimensions must be 1");
   }
   uint64_t length = passthroughBuf.shape[0];
   T* passthroughPtr = static_cast<T*>(passthroughBuf.ptr);
   std::span<T> passthroughSpan(passthroughPtr, length);
   return passthroughSpan;
}

/***********************************************************************************************************************
 * @brief Returns vector<T> from py:array_T<T>. Only works with py::array_t with 1 dimension.
 * 
 * @tparam T Type
 * @param passthrough[in] Numpy array to get as vector.
 * @return std::vector<T> clean and safe reference to contents of Numpy array.
 **********************************************************************************************************************/
template<typename T>
inline static std::vector<T> toVect1D(py::array_t<T, C_FMT> const & passthrough)
{
   py::buffer_info passthroughBuf = passthrough.request();
   if (passthroughBuf.ndim != 1) {
      throw std::runtime_error("Error. Number of dimensions must be 1");
   }
   uint64_t length = passthroughBuf.shape[0];
   std::vector<T> passthroughVector(passthrough.data(), passthrough.data() + length);
   return passthroughVector;
}

inline OCSC * convert_u8_ndarray_to_csc(py::array_t<uint8_t, F_FMT> const & po_ndarr)
{
   py::buffer_info py_pcm_bufinfo = po_ndarr.request();

   if (py_pcm_bufinfo.ndim != 2)
   {
      throw std::runtime_error("[ERROR] The pcm must be of ndim = 2!");
   }

   uint64_t u64_ndarr_rows = po_ndarr.shape(0L);

   std::span<uint8_t> const au8_ndarr_sp = toSpan2D(po_ndarr);

   OCSC * po_csc_ndarr = new OCSC(au8_ndarr_sp, u64_ndarr_rows);

   return po_csc_ndarr;
}

inline OCSR * convert_u8_ndarray_to_csr(py::array_t<uint8_t, F_FMT> const & po_ndarr)
{
   py::buffer_info py_pcm_bufinfo = po_ndarr.request();

   if (py_pcm_bufinfo.ndim != 2)
   {
      throw std::runtime_error("[ERROR] The pcm must be of ndim = 2!");
   }

   uint64_t u64_ndarr_rows = po_ndarr.shape(0L);

   std::span<uint8_t> const au8_ndarr_sp = toSpan2D(po_ndarr);

   OCSR * po_csr_ndarr = new OCSR(au8_ndarr_sp, u64_ndarr_rows);

   return po_csr_ndarr;
}

inline std::vector<uint8_t> mat_vec_mul(OCSR const & csr_mat, std::vector<uint8_t> const & vec)
{
   // I am not going to make any checks for the moment to gain speed and I will assume that the 
   // lengths of the matrix and vector are the correct ones.
   std::vector<uint8_t> vec_u8_res(csr_mat.get_row_num(), 0U);

   uint64_t u64_num_rows = csr_mat.get_row_num();

   for (uint64_t u64_row_idx = 0U; u64_row_idx < u64_num_rows; ++u64_row_idx)
   {
      std::span<uint64_t> sp_u8_col_idxs = csr_mat.get_row_col_idxs_fast(u64_row_idx);
      uint64_t u64_row_nnz = sp_u8_col_idxs.size();
      uint64_t u64_sum = 0U;
      for (uint64_t u64_i = 0U; u64_i < u64_row_nnz; ++u64_i)
      {
         u64_sum += vec[sp_u8_col_idxs[u64_i]];
      }
      vec_u8_res[u64_row_idx] = u64_sum % 2;
   }

   return vec_u8_res;
}

/***********************************************************************************************************************
 * CLASS METHODS
 **********************************************************************************************************************/
OBPOTF::SIDemData::SIDemData(SDemData const * const ps_ext_dem_data)
{
   if (this->check_members(ps_ext_dem_data))
   {
      this->af64_priors = toVect1D<double>(ps_ext_dem_data->af64_priors);
      this->po_obs_csr_mat = convert_u8_ndarray_to_csr(ps_ext_dem_data->po_obs_csr_mat);

      if (true == this->m_bp_bp_otf_enable)
      {
         this->po_transfer_csr_mat = convert_u8_ndarray_to_csr(ps_ext_dem_data->po_transfer_csr_mat);
         this->po_phen_pcm_csc = convert_u8_ndarray_to_csc(ps_ext_dem_data->po_phen_pcm_csc);
         this->po_phen_obs_csr = convert_u8_ndarray_to_csr(ps_ext_dem_data->po_phen_obs_csr);
      }
   }
   else
   {
      throw std::runtime_error("Error! Wrong external DEM data passed... Please verify the following conditions:\n"
                               "  - Priors are of type numpy.array[double]\n"
                               "  - Rest of the parameters are of type numpy.array[uint8]");
   }
}

bool OBPOTF::SIDemData::check_members(SDemData const * const ps_ext_dem_data)
{
   bool ret_val = false;

   ret_val = this->check_member<double>(ps_ext_dem_data->af64_priors);
   ret_val = (ret_val == false) ? ret_val : this->check_member<uint8_t>(ps_ext_dem_data->po_obs_csr_mat);

   // Check rest of values to enable 2 stage bp
   bool & p_enable_ref = this->m_bp_bp_otf_enable;
   p_enable_ref = 
      (ret_val == false) ? ret_val : this->check_member<uint8_t>(ps_ext_dem_data->po_transfer_csr_mat);
   p_enable_ref = 
      (p_enable_ref == false) ? p_enable_ref : this->check_member<uint8_t>(ps_ext_dem_data->po_phen_pcm_csc);
    p_enable_ref = 
      (p_enable_ref == false) ? p_enable_ref : this->check_member<uint8_t>(ps_ext_dem_data->po_phen_obs_csr);

   return ret_val;
}

OBPOTF::OBPOTF(py::object const & au8_pcm, float const & p,
               ENoiseType_t const noise_type = E_CC,
               py::object const & ps_ext_bp_iters = py::none(),
               SDemData_t const * ps_ext_dem_data = nullptr)
               :m_p(p)
{
   // Initialize depending the py::object instance
   if (true == py::isinstance<py::array_t<uint8_t>>(au8_pcm))
   {
      this->OBPOTF_init_from_numpy(au8_pcm, noise_type, ps_ext_bp_iters, ps_ext_dem_data);
   }
   else if (true == py::isinstance(au8_pcm, vf_scipy_csc_type))
   {
     this->OBPOTF_init_from_scipy_csc(au8_pcm, noise_type, ps_ext_bp_iters, ps_ext_dem_data);
   }
   else
   {
      throw std::runtime_error("Input type not supported! Input type must be ndarray of uint8"
                               "or scipy.sparse.csc_matrix");
   }

}

void OBPOTF::OBPOTF_init_from_scipy_csc(py::object const & au8_pcm,
                                          ENoiseType_t const & noise_type,
                                          py::object const & po_ext_bp_iters,
                                          SDemData_t const * const ps_ext_dem_data)
{
   // Convert scipy.sparse.csc_matrix to ndarray of uint8_t
   py::object dense_mat = au8_pcm.attr("toarray")();
   py::array_t<uint8_t, F_FMT> au8_pcm_pyarr = dense_mat.attr("astype")("uint8");

   this->OBPOTF_init_from_numpy(au8_pcm_pyarr, noise_type, po_ext_bp_iters, ps_ext_dem_data);
}

void OBPOTF::OBPOTF_init_from_numpy(py::array_t<uint8_t, F_FMT> const & au8_pcm,
                                    ENoiseType_t const & noise_type,
                                    py::object const & po_ext_bp_iters,
                                    SDemData_t const * const ps_ext_dem_data)
{
   py::buffer_info py_pcm_bufinfo = au8_pcm.request();

   if (py_pcm_bufinfo.ndim != 2)
   {
      throw std::runtime_error("[ERROR] The pcm must be of ndim = 2!");
   }

   m_u64_pcm_rows = au8_pcm.shape(0L);
   m_u64_pcm_cols = au8_pcm.shape(1L);


   // Span pointer to the pcm python input to be faster
   std::span<uint8_t> const au8_pcm_sp = toSpan2D(au8_pcm);

   // Create CSC type object.
   m_po_csc_mat = new OCSC(au8_pcm_sp, m_u64_pcm_rows);

   if (nullptr != ps_ext_dem_data)
   {
      m_ps_dem_data = new SIDemData(ps_ext_dem_data);
   }

   // Create BpSparse object
   m_po_bpsparse_pcm = new ldpc::bp::BpSparse(m_u64_pcm_rows, m_u64_pcm_cols);

   // Iterate through the pcm to get the non-trivial row indexes and maximum number of nt per column.
   for (uint64_t u64_c_idx = 0UL; u64_c_idx < m_u64_pcm_cols; ++u64_c_idx)
   {
      uint16_t u16_col_nnz = 0U;
      std::span<uint64_t> p_cur_col_idxs_sp = m_po_csc_mat->get_col_row_idxs_fast(u64_c_idx);
      u16_col_nnz = m_po_csc_mat->get_col_nnz(u64_c_idx);
      for (uint64_t u64_idx = 0UL; u64_idx < u16_col_nnz; ++u64_idx)
      {
         uint64_t u64_r_idx = p_cur_col_idxs_sp[u64_idx];
         m_po_bpsparse_pcm->insert_entry(u64_r_idx, u64_c_idx);
      }
      
      // TODO: Specify to decode CSS codes separately for capacity noise.
      if (u16_col_nnz == 1U)
      {
         m_po_csc_mat->add_row_idx(m_u64_pcm_rows, u64_c_idx);
         // // Take advantage of the integer division truncation.
         // if (p_cur_col_idxs_sp[0] < m_u64_pcm_rows / 2) 
         // {
         //    m_po_csc_mat->add_row_idx(m_u64_pcm_rows, u64_c_idx);
         // }
         // else
         // {
         //    m_po_csc_mat->add_row_idx(m_u64_pcm_rows+1, u64_c_idx);
         // }
      }
   }

   // Register decoding callback and select appropriate variables.
   std::vector<double> channel_errors_1;
   std::vector<double> channel_errors_2;
   if (noise_type == E_CC)
   {
      // Form the vector of indexes to be sorted
      m_au64_index_array = std::vector<uint64_t>(m_u64_pcm_cols);
      std::iota(m_au64_index_array.begin(), m_au64_index_array.end(), 0UL);
      
      m_po_bpsparse_otf = m_po_bpsparse_pcm;
      
      channel_errors_1 = std::vector<double>(m_u64_pcm_cols, m_p);
      channel_errors_2 = std::vector<double>(m_u64_pcm_cols, 1e-14);

      this->m_pf_decoding_func = &OBPOTF::bp_otf_cc_decode;
   }
   else if (noise_type == E_CLN && m_ps_dem_data->m_bp_bp_otf_enable)
   {
      uint64_t u64_phen_row_num = m_ps_dem_data->po_phen_pcm_csc->get_row_num();
      uint64_t u64_phen_col_num = m_ps_dem_data->po_phen_pcm_csc->get_col_num();
      m_po_bpsparse_phen_pcm = new ldpc::bp::BpSparse(u64_phen_row_num,
                                                      u64_phen_col_num);

      for (uint64_t u64_c_idx = 0UL; u64_c_idx < u64_phen_col_num; ++u64_c_idx)
      {
         uint16_t u16_col_nnz = 0U;
         std::span<uint64_t> p_cur_col_idxs_sp = m_ps_dem_data->po_phen_pcm_csc->get_col_row_idxs_fast(u64_c_idx);
         u16_col_nnz = m_ps_dem_data->po_phen_pcm_csc->get_col_nnz(u64_c_idx);
         for (uint64_t u64_i = 0UL; u64_i < p_cur_col_idxs_sp.size(); ++u64_i)
         {
            uint64_t u64_r_idx = p_cur_col_idxs_sp[u64_i];
            m_po_bpsparse_phen_pcm->insert_entry(u64_r_idx, u64_c_idx);
         }

         if (u16_col_nnz == 1U)
         {
            m_ps_dem_data->po_phen_pcm_csc->add_row_idx(u64_phen_row_num, u64_c_idx);
         }
      }

      // Form the vector of indexes to be sorted
      m_au64_index_array = std::vector<uint64_t>(u64_phen_col_num);
      std::iota(m_au64_index_array.begin(), m_au64_index_array.end(), 0UL);

      m_po_bpsparse_otf = m_po_bpsparse_phen_pcm;
      // Create the channel error probabilities vector from the obtained priors.
      channel_errors_1 = std::vector<double>(m_ps_dem_data->af64_priors.data(), 
                                             m_ps_dem_data->af64_priors.data()+m_ps_dem_data->af64_priors.size());
      channel_errors_2 = std::vector<double>(u64_phen_col_num, 1e-14);

      this->m_pf_decoding_func = &OBPOTF::bp_bp_otf_cln_decode;
   }
   else if (noise_type == E_CLN)
   {
      // Form the vector of indexes to be sorted
      m_au64_index_array = std::vector<uint64_t>(m_u64_pcm_cols);
      std::iota(m_au64_index_array.begin(), m_au64_index_array.end(), 0UL);

      m_po_bpsparse_otf = m_po_bpsparse_pcm;

      // Create the channel error probabilities vector from the obtained priors.
      channel_errors_1 = std::vector<double>(m_ps_dem_data->af64_priors.data(), 
                                             m_ps_dem_data->af64_priors.data()+m_ps_dem_data->af64_priors.size());
      channel_errors_2 = std::vector<double>(m_u64_pcm_cols, 1e-14);

      this->m_pf_decoding_func = &OBPOTF::bp_otf_cln_decode;
   }
   else
   {
      throw std::runtime_error("ERROR! Introduced noise type is not supported!\n"
                                 "If you would like to include support submit a PR or open an "
                                 "issue on the repository. Thanks!");
   }

   m_ps_bp_max_iterations = new SIBpMaxIters();
   if (false == po_ext_bp_iters.is_none() && true == py::isinstance<py::array_t<int32_t>>(po_ext_bp_iters))
   {
      std::span<int32_t> sp_bp_iters = toSpan1D<int32_t>(py::cast<py::array_t<int32_t, 1>>(po_ext_bp_iters));
      if (sp_bp_iters.size() == 3U)
      {
         if (sp_bp_iters[0] > 0)
         {
            m_ps_bp_max_iterations->m_pcm_bp_iters = sp_bp_iters[0];
         }
         else
         {
            std::cout << "WARNING: Invalid max iteration number (" << sp_bp_iters[0] << 
               ") inserted for stage 1! Defaulting to " << ST1_DEFAULT_MAX_BP << std::endl;
         }

         if (sp_bp_iters[1] > 0)
         {
            m_ps_bp_max_iterations->m_phen_bp_iters = sp_bp_iters[1];
         }
         else
         {
            std::cout << "WARNING: Invalid max iteration number (" << sp_bp_iters[1] << 
               ") inserted for stage 2! Defaulting to " << ST2_DEFAULT_MAX_BP << std::endl;
         }
         if (sp_bp_iters[2] > 0)
         {
            m_ps_bp_max_iterations->m_otf_bp_iters = sp_bp_iters[2];
         }
         else
         {
            std::cout << "WARNING: Invalid max iteration number (" << sp_bp_iters[2] << 
               ") inserted for stage 3! Defaulting to " << ST3_DEFAULT_MAX_BP << std::endl;
         }
      }
      else
      {
         std::cout << "WARNING: Invalid po_ext_bp_iters size introduced! Default values will be used.\n" <<
            "\tPlease enter 3 values, and insert -1 to use default values." << std::endl;
      }
   }
   else
   {
      std::cout << "WARNING: Invalid po_ext_bp_iters format inserted! Default values will be used." << std::endl;
   }

   // Create BpDecoder to use against the pcm
   m_po_pcm_bp = new ldpc::bp::BpDecoder(*m_po_bpsparse_pcm,
                                          channel_errors_1,
                                          m_ps_bp_max_iterations->m_pcm_bp_iters,
                                          // ldpc::bp::PRODUCT_SUM,
                                          ldpc::bp::MINIMUM_SUM,
                                          ldpc::bp::PARALLEL,
                                          1.0, 1e-14, 1,
                                          ldpc::bp::NULL_INT_VECTOR,
                                          0, true, ldpc::bp::SYNDROME);
   
   if (noise_type == E_CLN && m_ps_dem_data->m_bp_bp_otf_enable)
   {
      m_po_phen_bp = new ldpc::bp::BpDecoder(*m_po_bpsparse_phen_pcm,
                                                channel_errors_2,
                                                m_ps_bp_max_iterations->m_phen_bp_iters,
                                                // ldpc::bp::PRODUCT_SUM,
                                                ldpc::bp::MINIMUM_SUM,
                                                ldpc::bp::PARALLEL,
                                                1.0, 1e-14, 1,
                                                ldpc::bp::NULL_INT_VECTOR,
                                                0, true, ldpc::bp::SYNDROME);
   }

   // Create BpDecoder to use against the pcm after OTF.
   m_po_otf_bp = new ldpc::bp::BpDecoder(*m_po_bpsparse_otf,
                                          channel_errors_2,
                                          m_ps_bp_max_iterations->m_otf_bp_iters,
                                          ldpc::bp::PRODUCT_SUM,
                                          // ldpc::bp::MINIMUM_SUM,
                                          ldpc::bp::PARALLEL,
                                          1.0, 1e-14, 1,
                                          ldpc::bp::NULL_INT_VECTOR,
                                          0, true, ldpc::bp::SYNDROME);

}

/*
 * Debug purposes, eventually could be removed or improved...
 */
void OBPOTF::print_object(void)
{
   std::cout << "m_u64_pcm_rows: " << m_u64_pcm_rows << std::endl;
   std::cout << "m_u64_pcm_cols: " << m_u64_pcm_cols << std::endl;

   // std::cout << "m_au64_index_array: " << std::endl;
   // for (uint64_t u64_idx = 0U; u64_idx < m_au64_index_array.size(); ++u64_idx)
   //    std::cout << m_au64_index_array[u64_idx] << " ";
   // std::cout << std::endl;

   std::cout << "m_ps_bp_max_iterations->m_pcm_bp_iters: " << m_ps_bp_max_iterations->m_pcm_bp_iters << std::endl;
   std::cout << "m_ps_bp_max_iterations->m_pehn_bp_iters: " << m_ps_bp_max_iterations->m_phen_bp_iters<< std::endl;
   std::cout << "m_ps_bp_max_iterations->m_otf_bp_iters: " << m_ps_bp_max_iterations->m_otf_bp_iters << std::endl;

}

py::array_t<uint8_t> OBPOTF::decode(py::array_t<uint8_t, C_FMT> const & syndrome)
{
   return (this->*m_pf_decoding_func)(syndrome);
   // Below method is more "safe" due to readability... both work the same.
   // return std::invoke(this->m_pf_decoding_func, this, syndrome);
}

py::array_t<uint8_t> OBPOTF::bp_otf_cc_decode(py::array_t<uint8_t, C_FMT> const & syndrome)
{
   std::vector<uint8_t> u8_syndrome(syndrome.data(), syndrome.data() + syndrome.size());
   std::vector<uint8_t> u8_recovered_err = m_po_pcm_bp->decode(u8_syndrome);
   m_b_converged = m_po_pcm_bp->converge;

   if (false == m_po_pcm_bp->converge)
   {
      std::vector<double> llrs = m_po_pcm_bp->log_prob_ratios;
      std::vector<double> vec_f_probs = this->get_probs_from_llrs(llrs);

      std::vector<uint64_t> columns_chosen = this->otf_uf_probs(m_po_csc_mat, vec_f_probs);

      std::vector<double> updated_probs(m_u64_pcm_cols, 1e-14);
      uint64_t u64_col_chosen_sz = columns_chosen.size();
      for (uint64_t u64_idx = 0U; u64_idx < u64_col_chosen_sz; ++u64_idx)
         updated_probs[columns_chosen[u64_idx]] = m_p; // TODO: Think about this value

      // TODO: change channel probs to initial_logs and test...
      m_po_otf_bp->update_channel_probs(updated_probs);
      u8_recovered_err = m_po_otf_bp->decode(u8_syndrome);
      m_b_converged = m_po_otf_bp->converge;
   }

   return as_pyarray(std::move(u8_recovered_err));
}

py::array_t<uint8_t> OBPOTF::bp_otf_cln_decode(py::array_t<uint8_t, C_FMT> const & syndrome)
{
   std::vector<uint8_t> vec_u8_res;
   // First BP stage
   std::vector<uint8_t> u8_syndrome(syndrome.data(), syndrome.data() + syndrome.size());

   START_CHRONO
   std::vector<uint8_t> u8_recovered_err = m_po_pcm_bp->decode(u8_syndrome);
   STOP_CHRONO("First BP Stage decode: ")
   m_b_converged = m_po_pcm_bp->converge;

   if (false == m_po_pcm_bp->converge)
   {
      // OTF stage
      std::vector<double> vec_f_llrs = m_po_pcm_bp->log_prob_ratios;
      std::vector<double> vec_f_probs = this->get_probs_from_llrs(vec_f_llrs);
      START_CHRONO
      // std::vector<uint64_t> columns_chosen = this->otf_uf(m_ps_dem_data->po_phen_pcm_csc, vec_f_llrs);
      std::vector<uint64_t> columns_chosen = this->otf_uf_probs(m_po_csc_mat, vec_f_probs);
      STOP_CHRONO("OTF: ")

      // std::vector<double> updated_llrs(m_ps_dem_data->po_phen_pcm_csc->get_col_num() , vfc_initial_llr_value);
      std::vector<double> updated_llrs(m_po_csc_mat->get_col_num() , 1e-9);
      uint64_t u64_col_chosen_sz = columns_chosen.size();
      // std::cout << "CPP OTF column chosen num: " << u64_col_chosen_sz << std::endl;
      for (uint64_t u64_idx = 0U; u64_idx < u64_col_chosen_sz; ++u64_idx)
      {
         uint64_t u64_col_idx = columns_chosen[u64_idx];
         // updated_llrs[u64_col_idx] = vec_f_llrs[u64_col_idx];
         updated_llrs[u64_col_idx] = vec_f_probs[u64_col_idx];
      }
      // m_po_otf_bp->initial_log_prob_ratios = updated_llrs;
      m_po_otf_bp->update_channel_probs(updated_llrs);

      START_CHRONO
      u8_recovered_err = m_po_otf_bp->decode(u8_syndrome);
      STOP_CHRONO("OTF stage decode: ")
      m_b_converged = m_po_otf_bp->converge;
   }

   START_CHRONO
   vec_u8_res = mat_vec_mul(*(m_ps_dem_data->po_obs_csr_mat), u8_recovered_err);
   STOP_CHRONO("Observables check multiplication: ")

   return as_pyarray(std::move(vec_u8_res));
}

py::array_t<uint8_t> OBPOTF::bp_bp_otf_cln_decode(py::array_t<uint8_t, C_FMT> const & syndrome)
{
   std::vector<uint8_t> vec_u8_res;
   // First BP stage
   std::vector<uint8_t> u8_syndrome(syndrome.data(), syndrome.data() + syndrome.size());

   START_CHRONO
   std::vector<uint8_t> u8_recovered_err = m_po_pcm_bp->decode(u8_syndrome);
   STOP_CHRONO("First Stage decode: ")
   m_b_converged = m_po_pcm_bp->converge;

   if (false == m_po_pcm_bp->converge)
   {
      // Second BP stage in case it does not converge.
      std::vector<double> vec_f_llrs = m_po_pcm_bp->log_prob_ratios;

      START_CHRONO
      std::vector<double> vec_f_mapped_probs = this->propagate(vec_f_llrs);
      STOP_CHRONO("Propagate: ")
      for (auto &prob : vec_f_mapped_probs) {
         prob = std::clamp(prob, 1e-40, 1 - 1e-40);
      }
      std::cout << "Dimensions of vec_f_mapped_probs: " << vec_f_mapped_probs.size() << std::endl;
      auto minmax = std::minmax_element(vec_f_mapped_probs.begin(), vec_f_mapped_probs.end());
      std::cout << "Minimum value in vec_f_mapped_probs: " << *minmax.first << std::endl;
      std::cout << "Maximum value in vec_f_mapped_probs: " << *minmax.second << std::endl;

      // m_po_phen_bp->channel_probabilities = vec_f_mapped_probs;
      m_po_phen_bp->update_channel_probs(vec_f_mapped_probs);

      START_CHRONO
      u8_recovered_err = m_po_phen_bp->decode(u8_syndrome);
      STOP_CHRONO("Second Stage decode: ")
      m_b_converged = m_po_phen_bp->converge;
      
      if (false == m_po_phen_bp->converge)
      {
         // Third stage with OTF
         vec_f_llrs = m_po_phen_bp->log_prob_ratios;
         std::vector<double> vec_f_probs = this->get_probs_from_llrs(vec_f_llrs);
         START_CHRONO
         // std::vector<uint64_t> columns_chosen = this->otf_uf(m_ps_dem_data->po_phen_pcm_csc, vec_f_llrs);
         std::vector<uint64_t> columns_chosen = this->otf_uf_probs(m_ps_dem_data->po_phen_pcm_csc, vec_f_probs);
         STOP_CHRONO("OTF: ")

         // std::vector<double> updated_llrs(m_ps_dem_data->po_phen_pcm_csc->get_col_num() , vfc_initial_llr_value);
         std::vector<double> updated_llrs(m_ps_dem_data->po_phen_pcm_csc->get_col_num() , 1e-9);
         uint64_t u64_col_chosen_sz = columns_chosen.size();
         // std::cout << "CPP OTF column chosen num: " << u64_col_chosen_sz << std::endl;
         for (uint64_t u64_idx = 0U; u64_idx < u64_col_chosen_sz; ++u64_idx)
         {
            uint64_t u64_col_idx = columns_chosen[u64_idx];
            // updated_llrs[u64_col_idx] = vec_f_llrs[u64_col_idx];
            updated_llrs[u64_col_idx] = vec_f_probs[u64_col_idx];
         }
         // m_po_otf_bp->initial_log_prob_ratios = updated_llrs;
         m_po_otf_bp->update_channel_probs(updated_llrs);

         START_CHRONO
         u8_recovered_err = m_po_otf_bp->decode(u8_syndrome);
         STOP_CHRONO("Third Stage decode: ")
         m_b_converged = m_po_otf_bp->converge;
         std::cout << "OTF converged: " << m_b_converged << std::endl;
      }
      
      START_CHRONO
      vec_u8_res = mat_vec_mul(*(m_ps_dem_data->po_phen_obs_csr), u8_recovered_err);
      STOP_CHRONO("Observables check multiplication: ")
   }
   else
   {
      START_CHRONO
      vec_u8_res = mat_vec_mul(*(m_ps_dem_data->po_obs_csr_mat), u8_recovered_err);
      STOP_CHRONO("Observables check multiplication: ")
   }

   return as_pyarray(std::move(vec_u8_res));
}

std::vector<uint64_t> OBPOTF::sort_probs_indexes(py::array_t<double> const & probs) 
{
   // COPY! Not happy...
   std::vector<uint64_t> idx = m_au64_index_array;

   std::sort(idx.begin(), idx.end(), 
               [&probs](uint64_t i1, uint64_t i2) 
               {
                  return std::greater<double>{}(probs.data()[i1], probs.data()[i2]);
               }
            );

   return idx;
}

std::vector<uint64_t *> OBPOTF::sort_probs_indexes_nc(std::vector<double> const & probs) 
{
   std::vector<uint64_t *> idx(m_au64_index_array.size());

   uint64_t u64_idx_sz = idx.size();
   for (uint64_t u64_idx = 0U; u64_idx < u64_idx_sz; ++u64_idx)
      idx[u64_idx] = m_au64_index_array.data() + u64_idx;

   // Using std::sort assures a worst-case scenario of time complexity O(n*log(n))
   std::sort(idx.begin(), idx.end(), 
               [&probs](uint64_t * i1, uint64_t * i2) 
               {
                  // return std::less<double>{}(llrs[*i1], llrs[*i2]);
                  return std::greater<double>{}(probs[*i1], probs[*i2]);
               }
            );

   return idx;
}

std::vector<uint64_t *> OBPOTF::sort_probs_indexes_nc(std::span<double> const & probs) 
{
   std::vector<uint64_t *> idx(m_au64_index_array.size());

   uint64_t u64_idx_sz = idx.size();
   for (uint64_t u64_idx = 0U; u64_idx < u64_idx_sz; ++u64_idx)
      idx[u64_idx] = m_au64_index_array.data() + u64_idx;

   // Using std::sort assures a worst-case scenario of time complexity O(n*log(n))
   std::sort(idx.begin(), idx.end(), 
               [&probs](uint64_t * i1, uint64_t * i2) 
               {
                  return std::greater<double>{}(probs[*i1], probs[*i2]);
               }
            );

   return idx;
}

std::vector<uint64_t> OBPOTF::otf_classical_uf_probs(OCSC const * const po_csc_mat, std::vector<double> const & probs)
{
   uint16_t const & hog_rows = po_csc_mat->get_row_num();
   uint64_t const & hog_cols = po_csc_mat->get_col_num();
   
   std::vector<uint64_t> columns_chosen;
   columns_chosen.reserve(hog_cols);

   std::vector<uint64_t *> sorted_idxs = this->sort_probs_indexes_nc(probs);

   DisjSet clstr_set = DisjSet(hog_rows);

   uint64_t u64_sorted_idxs_sz = sorted_idxs.size();
   for (uint64_t col_idx = 0UL; col_idx < u64_sorted_idxs_sz; ++col_idx)
   {
      uint64_t effective_col_idx = *(sorted_idxs[col_idx]);
      //uint64_t col_offset = effective_col_idx * hog_rows;
      std::span<uint64_t> column_sp = m_po_csc_mat->get_col_row_idxs_fast(effective_col_idx);

      long int retcode = -1;
      long int root_set = clstr_set.find(column_sp[0]);
      uint64_t u64_col_sp_sz = column_sp.size();
      for (uint64_t nt_elem_idx = 1UL; nt_elem_idx < u64_col_sp_sz; ++nt_elem_idx)
      {
         // Could be of interest to add some kind of flag to indicate end of the column when 
         // column elements are different across the pcm
         // if (column_sp[nt_elem_idx] == -1L)
         //    break;

         retcode = clstr_set.set_union_opt(root_set, column_sp[nt_elem_idx]);
         if (retcode == -1L)
            break;
         root_set = retcode;
      }
      if (retcode != -1L)
      {
         columns_chosen.push_back(effective_col_idx);
         // if (columns_chosen.size() == rank)
         //    break;
      }
   }

   return columns_chosen;
}

std::vector<uint64_t> OBPOTF::otf_uf_probs(OCSC const * const po_csc_mat, std::vector<double> const & probs)
{
   uint16_t const & csc_rows = po_csc_mat->get_row_num(); // Already accounts for virtual checks
   uint64_t const & csc_cols = po_csc_mat->get_col_num();
   
   std::vector<uint64_t> columns_chosen;
   columns_chosen.reserve(csc_cols);

   std::vector<uint64_t *> sorted_idxs = this->sort_probs_indexes_nc(probs);

   DisjSet clstr_set = DisjSet(csc_rows);

   uint64_t u64_sorted_idxs_sz = sorted_idxs.size();
   for (uint64_t col_idx = 0UL; col_idx < u64_sorted_idxs_sz; ++col_idx)
   {
      uint64_t effective_col_idx = *(sorted_idxs[col_idx]);
      std::span<uint64_t> column_sp = po_csc_mat->get_col_row_idxs_fast(effective_col_idx);
      
      std::vector<uint8_t> checker(csc_rows, 0);
      std::vector<int> depths = {0, -1, 0};
      bool boolean_condition = true;

      uint64_t u64_col_sp_sz = column_sp.size();
      for (uint64_t nt_elem_idx = 0UL; nt_elem_idx < u64_col_sp_sz; ++nt_elem_idx)
      {
         int elem_root = clstr_set.find(column_sp[nt_elem_idx]);
         int elem_depth = clstr_set.get_rank(elem_root);

         if (checker[elem_root] == 1U)
         {
            boolean_condition = false;
            break;
         }
         checker[elem_root] = 1U;

         if (elem_depth > depths[1])
         {
            depths = {elem_root, elem_depth, 0};
         }
         if (elem_depth == depths[1])
         {
            depths[2] = 1;
         }
      }

      if (boolean_condition == true)
      {
         uint64_t u64_checker_sz = checker.size();
         for(uint64_t elem = 0UL; elem < u64_checker_sz; ++elem)
         {
            if (checker[elem] == 1U)
               clstr_set.set_parent(elem, depths[0]);
         }
         columns_chosen.push_back(effective_col_idx);
         if (depths[2] == 1)
         {
            clstr_set.increase_rank(depths[0]);
         }
         // if (columns_chosen.size() == rank)
         //    break;
      }
   }

   return columns_chosen;
}

py::array_t<uint64_t> OBPOTF::otf_uf_probs(py::array_t<double, C_FMT> const & probs)
{
   uint16_t const & csc_rows = m_po_csc_mat->get_row_num(); // Already accounts for virtual checks
   uint64_t const & csc_cols = m_po_csc_mat->get_col_num();
   
   std::vector<uint64_t> columns_chosen;
   columns_chosen.reserve(csc_cols);

   std::span<double> llrs_sp = toSpan1D<double>(probs);
   std::vector<uint64_t *> sorted_idxs = this->sort_probs_indexes_nc(llrs_sp);

   DisjSet clstr_set = DisjSet(csc_rows);

   uint64_t u64_sorted_idxs_sz = sorted_idxs.size();
   for (uint64_t col_idx = 0UL; col_idx < u64_sorted_idxs_sz; ++col_idx)
   {
      uint64_t effective_col_idx = *(sorted_idxs[col_idx]);
      std::span<uint64_t> column_sp = m_po_csc_mat->get_col_row_idxs_fast(effective_col_idx);
      
      std::vector<uint8_t> checker(csc_rows, 0);
      std::vector<int> depths = {0, -1, 0};
      bool boolean_condition = true;

      uint64_t u64_col_sp_sz = column_sp.size();
      for (uint64_t nt_elem_idx = 0UL; nt_elem_idx < u64_col_sp_sz; ++nt_elem_idx)
      {
         int elem_root = clstr_set.find(column_sp[nt_elem_idx]);
         int elem_depth = clstr_set.get_rank(elem_root);

         if (checker[elem_root] == 1U)
         {
            boolean_condition = false;
            break;
         }
         checker[elem_root] = 1U;

         if (elem_depth > depths[1])
         {
            depths = {elem_root, elem_depth, 0};
         }
         if (elem_depth == depths[1])
         {
            depths[2] = 1;
         }
      }

      if (boolean_condition == true)
      {
         uint64_t u64_checker_sz = checker.size();
         for(uint64_t elem = 0UL; elem < u64_checker_sz; ++elem)
         {
            if (checker[elem] == 1U)
               clstr_set.set_parent(elem, depths[0]);
         }
         columns_chosen.push_back(effective_col_idx);
         if (depths[2] == 1)
         {
            clstr_set.increase_rank(depths[0]);
         }
         // if (columns_chosen.size() == rank)
         //    break;
      }
   }

   return as_pyarray(std::move(columns_chosen));
}

double OBPOTF::compute_probability_from_log(double const & f64_log_val)
{
   double f64_res_prob = 1e-40;

   f64_res_prob = 1.0 / (1.0 + std::exp(f64_log_val));
   if (f64_res_prob > (1.0 - 1e-40))
   {
      f64_res_prob = 1.0 - 1e-40;
   }
   else if (f64_res_prob < 1e-40)
   {
      f64_res_prob = 1e-40;
   }

   return f64_res_prob;
}

std::vector<double> OBPOTF::get_probs_from_llrs(std::vector<double> const & vec_f_llrs)
{
   std::vector<double> vec_f_probs;
   vec_f_probs.resize(vec_f_llrs.size());
   uint64_t u64_probs_len = vec_f_probs.size();
   // To get actual probabilities from the log_prob_ratios
   for (uint64_t u64_idx = 0U; u64_idx < u64_probs_len; ++u64_idx)
   {
      vec_f_probs[u64_idx] = this->compute_probability_from_log(vec_f_llrs[u64_idx]);
   }

   return vec_f_probs;
}

std::vector<double> OBPOTF::propagate(std::vector<double> const & vec_f_llrs)
{
   std::vector<double> vec_f_probs = this->get_probs_from_llrs(vec_f_llrs);
   std::vector<double> vec_f_mapped_probs;
   uint64_t u64_res_len = m_ps_dem_data->po_transfer_csr_mat->get_row_num();
   vec_f_mapped_probs.resize(u64_res_len);

   for (uint64_t u64_row_idx = 0U; u64_row_idx < u64_res_len; ++u64_row_idx)
   {
      double f_prod_sum = 1.0;
      std::span<uint64_t> sp_u64_curr_row = m_ps_dem_data->po_transfer_csr_mat->get_row_col_idxs_fast(u64_row_idx);
      uint64_t u64_curr_row_len = sp_u64_curr_row.size();
      for (uint64_t u64_idx = 0U; u64_idx < u64_curr_row_len; ++u64_idx)
      {
         f_prod_sum *= (1.0 - (2.0 * vec_f_probs[sp_u64_curr_row[u64_idx]]));
      }

      vec_f_mapped_probs[u64_row_idx] = 0.5 * (1.0 - f_prod_sum);
   }

   return vec_f_mapped_probs;
}

#if defined(DEBUG_OBPOTF)
py::array_t<uint8_t> OBPOTF::getPcm(void)
{
   return as_pyarray(std::move(m_po_csc_mat->expand_to_row_major()));
}

py::array_t<uint8_t> OBPOTF::getPhenPcm(void)
{
   return as_pyarray(std::move(m_ps_dem_data->po_phen_pcm_csc->expand_to_row_major()));
}

py::array_t<uint8_t> OBPOTF::getObs(void)
{
   return as_pyarray(std::move(m_ps_dem_data->po_obs_csr_mat->expand_to_row_major()));
}

py::array_t<uint8_t> OBPOTF::getTransfMat(void)
{
   return as_pyarray(std::move(m_ps_dem_data->po_transfer_csr_mat->expand_to_row_major()));
}

py::array_t<double> OBPOTF::getPriors(void)
{
   return as_pyarray(std::move(m_ps_dem_data->af64_priors));
}
#endif

OBPOTF::~OBPOTF(void)
{
   if (nullptr != m_po_csc_mat)
      delete m_po_csc_mat;

   if (nullptr != m_ps_bp_max_iterations)
      delete m_ps_bp_max_iterations;

   if (nullptr != m_ps_dem_data)
   {
      if (nullptr != m_ps_dem_data->po_transfer_csr_mat)
         delete m_ps_dem_data->po_transfer_csr_mat;
      
      if (nullptr != m_ps_dem_data->po_obs_csr_mat)
         delete m_ps_dem_data->po_obs_csr_mat;
      
      if (nullptr != m_ps_dem_data->po_phen_pcm_csc)
         delete m_ps_dem_data->po_phen_pcm_csc;
      
      if (nullptr != m_ps_dem_data->po_phen_obs_csr)
         delete m_ps_dem_data->po_phen_obs_csr;
      
      delete m_ps_dem_data;
   }

   if (nullptr != m_po_pcm_bp)
      delete m_po_pcm_bp;

   if (nullptr != m_po_phen_bp)
      delete m_po_phen_bp;

   if (nullptr != m_po_otf_bp)
      delete m_po_otf_bp;

   if (nullptr != m_po_bpsparse_pcm)
      delete m_po_bpsparse_pcm;
}