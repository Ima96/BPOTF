/***********************************************************************************************************************
 * @file    OCSR.h
 * @author  Imanol Etxezarreta (ietxezarretam@gmail.com)
 * 
 * @brief   Object declaration to handle Compressed Sparse Row matrices. This object is a particularization used
 *          for binary low density parity check matrices.
 * 
 * @version 0.1
 * @date    05/11/2024
 * 
 * @copyright Copyright (c) 2024
 **********************************************************************************************************************/
#ifndef OCSR_H_
#define OCSR_H_

#include <inttypes.h>
#include <vector>
#include <span> // C++20

class OCSR
{
   private:
      /* data */
      //! 
      uint64_t * m_pu64_indptr;
      //! Column indexes
      uint64_t * m_pu64_c_indices;
      //! Number of non-zero values
      uint64_t m_u64_nnz;
      //! Row number
      uint64_t m_u64_m;
      //! Column number
      uint64_t m_u64_n;

   private:

      void add_col_idx_entry(uint64_t const & u64_row_idx, uint64_t const & u64_col_idx);

   public:
      OCSR() = delete;

      OCSR(OCSR const & csr_mat);

      OCSR(std::vector<uint64_t> const & u64_col_indices,
            std::vector<uint64_t> const & u64_ext_indptr,
            uint64_t const & u64_nnz);

      OCSR(std::vector<uint8_t> const & pcm, uint64_t const & u64_row_num);

      OCSR(std::span<uint8_t> const & pcm, uint64_t const & u64_row_num);

      OCSR(std::vector<std::vector<uint8_t>> const & pcm);

      void print_csr(void);

      inline uint64_t get_col_num(void) const
      {
         return m_u64_n;
      };

      inline uint64_t get_row_num(void) const
      {
         return m_u64_m;
      };

      inline uint64_t get_nnz(void) const
      {
         return m_u64_nnz;
      };

      std::vector<std::vector<uint8_t>> expand_to_mat(void) const;

      std::vector<uint8_t> expand_to_column_major(void);

      std::vector<uint8_t> expand_to_row_major(void);

      uint64_t get_row_nnz(uint64_t const & u64_row) const;

      std::vector<uint64_t> get_row_col_idxs(uint64_t const & u64_row);

      std::span<uint64_t> get_row_col_idxs_fast(uint64_t const & u64_row) const;

      void add_col_idx(uint64_t const & u64_row_idx, uint64_t const & u64_col_idx);

      ~OCSR();
};

#endif // OCSR_H_