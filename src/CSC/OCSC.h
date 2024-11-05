/***********************************************************************************************************************
 * @file    OCSC.h
 * @author  Imanol Etxezarreta (ietxezarretam@gmail.com)
 * 
 * @brief   Object declaration to handle Compressed Sparse Column matrices. This object is a particularization used
 *          for binary low density parity check matrices.
 * 
 * @version 0.1
 * @date    28/05/2024
 * 
 * @copyright Copyright (c) 2024
 **********************************************************************************************************************/
#ifndef OCSC_H_
#define OCSC_H_

#include <inttypes.h>
#include <vector>
#include <span> // C++20

class OCSC
{
   private:
      /* data */
      uint64_t * m_pu64_indptr;
      uint64_t * m_pu64_r_indices;
      uint64_t m_u64_nnz;
      uint64_t m_u64_m;
      uint64_t m_u64_n;

   private:

      void add_row_idx_entry(uint64_t const & u64_row_idx, uint64_t const & u64_col_idx);

   public:
      OCSC() = delete;

      OCSC(OCSC const & csc_mat);

      OCSC(std::vector<uint64_t> const & u64_nnz_cols_vec,
            std::vector<uint64_t> const & u64_row_idxs_vec,
            uint64_t const & u64_nnz);

      OCSC(std::vector<uint8_t> const & pcm, uint64_t const & u64_row_num);

      OCSC(std::span<uint8_t> const & pcm, uint64_t const & u64_row_num);

      OCSC(std::vector<std::vector<uint8_t>> const & pcm);

      void print_csc(void);

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

      uint64_t get_col_nnz(uint64_t const & u64_col);

      std::vector<uint64_t> get_col_row_idxs(uint64_t const & u64_col);

      std::span<uint64_t> get_col_row_idxs_fast(uint64_t const & u64_col);

      void add_row_idx(uint64_t const & u64_row_idx, uint64_t const & u64_col_idx);

      ~OCSC();
};




#endif // OCSC_H_