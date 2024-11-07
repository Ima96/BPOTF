/***********************************************************************************************************************
 * @file    OCSR.cpp
 * @author  Imanol Etxezarreta (ietxezarretam@gmail.com)
 * 
 * @brief   
 * 
 * @version 0.1
 * @date    05/11/2024
 * 
 * @copyright Copyright (c) 2024
 **********************************************************************************************************************/

#include <iostream>
#include <cstring>
#include <algorithm>

#include "OCSR.h"

/***********************************************************************************************************************
 * PRIVATE FUNCTIONS
 **********************************************************************************************************************/
void OCSR::add_col_idx_entry(uint64_t const & u64_row_idx, uint64_t const & u64_col_idx)
{
   std::vector<uint64_t> au64_row_checks = this->get_row_col_idxs(u64_row_idx);
   // Check if already exist
   if (std::find(au64_row_checks.begin(), au64_row_checks.end(), u64_col_idx) != au64_row_checks.end())
   {
      return;
   }

   uint64_t * pu64_temp_buffer = new uint64_t[m_u64_nnz+1];

   uint16_t u16_count = 0;
   for (uint16_t u16_idx = 0, u16_end_cond = au64_row_checks.size(); u16_idx < u16_end_cond; ++u16_idx)
   {
      if (au64_row_checks[u16_idx] < u64_col_idx)
      {
         ++u16_count;
      }
      else
      {
         break;
      }
   }

   uint64_t u64_num_entries_until_insert = m_pu64_indptr[u64_row_idx] + u16_count;
   std::memcpy(pu64_temp_buffer, m_pu64_c_indices, u64_num_entries_until_insert * sizeof(uint64_t));
   pu64_temp_buffer[u64_num_entries_until_insert] = u64_col_idx;
   uint64_t u64_rest = m_u64_nnz-u64_num_entries_until_insert;
   std::memcpy(&pu64_temp_buffer[u64_num_entries_until_insert]+1, 
                  &m_pu64_c_indices[u64_num_entries_until_insert], 
                  u64_rest * sizeof(uint64_t));
   
   delete [] m_pu64_c_indices;
   m_pu64_c_indices = pu64_temp_buffer;

   for (uint64_t u64_idx = u64_row_idx + 1; u64_idx < m_u64_m + 1; ++u64_idx)
   {
      ++m_pu64_indptr[u64_idx];
   }

   ++m_u64_nnz;
   if (u64_col_idx+1 > m_u64_n)
   {
      m_u64_n = u64_col_idx+1;
   }
}


/***********************************************************************************************************************
 * PUBLIC FUNCTIONS
 **********************************************************************************************************************/
OCSR::OCSR(OCSR const & csr_mat):
   m_u64_nnz(csr_mat.m_u64_nnz),
   m_u64_m(csr_mat.m_u64_m),
   m_u64_n(csr_mat.m_u64_n)
{
   m_pu64_indptr = new uint64_t[m_u64_m+1];
   m_pu64_c_indices = new uint64_t[m_u64_nnz];

   std::memcpy(m_pu64_indptr, csr_mat.m_pu64_indptr, (m_u64_m+1) * sizeof(uint64_t));
   std::memcpy(m_pu64_c_indices, csr_mat.m_pu64_c_indices, m_u64_nnz * sizeof(uint64_t));
}

OCSR::OCSR(std::vector<uint64_t> const & u64_col_indices,
            std::vector<uint64_t> const & u64_ext_indptr,
            uint64_t const & u64_nnz):
            m_u64_nnz(u64_nnz)
{
   m_u64_m = u64_ext_indptr.size() - 1;
   // I think it is not real if the last column is all 0s but al least is something
   m_u64_n = *std::max_element(u64_col_indices.begin(), u64_col_indices.end()) + 1;

   m_pu64_indptr= new uint64_t[m_u64_m+1];
   m_pu64_c_indices = new uint64_t[m_u64_nnz];

   std::memcpy(m_pu64_indptr, u64_ext_indptr.data(), (m_u64_m+1) * sizeof(uint64_t));
   std::memcpy(m_pu64_c_indices, u64_col_indices.data(), m_u64_nnz * sizeof(uint64_t));
}

OCSR::OCSR(std::vector<std::vector<uint8_t>> const & pcm)
{
   uint64_t const u64_m = pcm.size(); // rows
   uint64_t const u64_n = (u64_m > 0) ? pcm[0].size() : 0U; // cols

   m_u64_m = u64_m;
   m_u64_n = u64_n;
   m_u64_nnz = 0U;

   m_pu64_indptr = new uint64_t[u64_m+1];
   m_pu64_indptr[0] = 0;

   // count number of non-zero elements
   for (uint64_t i = 0U; i < u64_m; i++)
   {
      for (uint64_t j = 0U; j < u64_n; j++)
      {
         if (pcm[i][j] != 0U)
         {
            m_u64_nnz++;
         }
      }
      m_pu64_indptr[i+1] = m_u64_nnz;
   }

   // allocate memory
   m_pu64_c_indices = new uint64_t[m_u64_nnz];

   // fill array
   uint64_t u64_cont = 0U;
   for (uint64_t i = 0U; i < u64_m; i++)
   {
      for (uint64_t j = 0U; j < u64_n; j++)
      {
         if (pcm[i][j] != 0U)
         {
            m_pu64_c_indices[u64_cont] = j;
            u64_cont++;
         }
      }
   }
}

// This would be a constructor from a Column-Major vector
OCSR::OCSR(std::vector<uint8_t> const & pcm, uint64_t const & u64_row_num)
{
   if (pcm.size() % u64_row_num != 0)
   {
      throw std::runtime_error("Error. Number of rows must be equal in all columns.");
   }

   m_u64_m = u64_row_num;
   m_u64_n = pcm.size() / u64_row_num;
   m_u64_nnz = 0U;

   m_pu64_indptr = new uint64_t[m_u64_m+1];
   m_pu64_indptr[0] = 0;

   // count number of non-zero elements
   for (uint64_t i = 0U; i < m_u64_m; i++)
   {
      for (uint64_t j = 0U; j < m_u64_n; j++)
      {
         if (pcm[(j*m_u64_m)+i] != 0U)
         {
            m_u64_nnz++;
         }
      }
      m_pu64_indptr[i+1] = m_u64_nnz;
   }

   // allocate memory
   m_pu64_c_indices = new uint64_t[m_u64_nnz];

   // fill array
   uint64_t u64_cont = 0U;
   for (uint64_t i = 0U; i < m_u64_m; ++i)
   {
      for (uint64_t j = 0U; j < m_u64_n; ++j)
      {
         if (pcm[(j*m_u64_m)+i] != 0U)
         {
            m_pu64_c_indices[u64_cont] = j;
            u64_cont++;
         }
      }
   }
}

OCSR::OCSR(std::span<uint8_t> const & pcm, uint64_t const & u64_row_num)
{
   if (pcm.size() % u64_row_num != 0)
   {
      throw std::runtime_error("Error. Number of rows must be equal in all columns.");
   }

   m_u64_m = u64_row_num;
   m_u64_n = pcm.size() / u64_row_num;
   m_u64_nnz = 0U;

   m_pu64_indptr = new uint64_t[m_u64_m+1];
   m_pu64_indptr[0] = 0;

   // count number of non-zero elements
   for (uint64_t i = 0U; i < m_u64_m; i++)
   {
      for (uint64_t j = 0U; j < m_u64_n; j++)
      {
         if (pcm[(j*m_u64_m)+i] != 0U)
         {
            m_u64_nnz++;
         }
      }
      m_pu64_indptr[i+1] = m_u64_nnz;
   }

   // allocate memory
   m_pu64_c_indices = new uint64_t[m_u64_nnz];

   // fill array
   uint64_t u64_cont = 0U;
   for (uint64_t i = 0U; i < m_u64_m; ++i)
   {
      for (uint64_t j = 0U; j < m_u64_n; ++j)
      {
         if (pcm[(j*m_u64_m)+i] != 0U)
         {
            m_pu64_c_indices[u64_cont] = j;
            u64_cont++;
         }
      }
   }
}

OCSR::~OCSR()
{
   delete [] m_pu64_c_indices;
   delete [] m_pu64_indptr;
}

void OCSR::print_csr(void)
{
   std::cout << "Shape (MxN): " << m_u64_m << "x" << m_u64_n << std::endl;
   std::cout << "m_u64_nnz: " << m_u64_nnz << std::endl;
   std::cout << "m_pu64_c_indices: [";
   for (uint64_t i = 0; i < m_u64_nnz; i++)
   {
      std::string str = ", ";
      if (i == m_u64_nnz-1)
         str = "";
      std::cout << m_pu64_c_indices[i] << str;
   }
   std::cout << "]\n";

   std::cout << "m_pu64_indptr: [";
   for (uint64_t i = 0; i < m_u64_m+1; i++)
   {
      std::string str = ", ";
      if (i == m_u64_m)
         str = "";
      std::cout << m_pu64_indptr[i] << str;
   }
   std::cout << "]\n";

}

std::vector<std::vector<uint8_t>> OCSR::expand_to_mat(void) const
{
   std::vector<std::vector<uint8_t>> res_mat(m_u64_m, std::vector<uint8_t>(m_u64_n, 0U));

   for (uint64_t i = 1U; i < m_u64_m+1; i++)
   {
      for (uint64_t j = m_pu64_indptr[i-1]; j < m_pu64_indptr[i]; j++)
      {
         uint64_t col_idx = m_pu64_c_indices[j];
         res_mat[i-1][col_idx] = 1U;
      }
   }

   return res_mat;
}

std::vector<uint8_t> OCSR::expand_to_row_major(void)
{
   std::vector<uint8_t> res_vec(m_u64_m * m_u64_n, 0U);

   for (uint64_t i = 1U; i < m_u64_m+1; i++)
   {
      for (uint64_t j = m_pu64_indptr[i-1]; j < m_pu64_indptr[i]; j++)
      {
         uint64_t col_idx = m_pu64_c_indices[j];
         res_vec[((i-1)*m_u64_n)+col_idx] = 1U;
      }
   }

   return res_vec;
}

std::vector<uint8_t> OCSR::expand_to_column_major(void)
{
   std::vector<uint8_t> res_vec(m_u64_m * m_u64_n, 0U);

   std::vector<std::vector<uint8_t>> ppu8_exp_mat = this->expand_to_mat();

   for (uint64_t u64_i = 0U; u64_i < m_u64_m; ++u64_i)
   {
      std::vector<uint64_t> vec_u8_col_idxs = this->get_row_col_idxs(u64_i);
      uint64_t u64_row_nnz = vec_u8_col_idxs.size();
      for (uint64_t u64_idx = 0U; u64_idx < u64_row_nnz; ++u64_idx)
      {
         res_vec[(vec_u8_col_idxs[u64_idx]*m_u64_m)+u64_i] = 1U;
      }
   }

   return res_vec;
}

uint64_t OCSR::get_row_nnz(uint64_t const & u64_row) const
{
   if (u64_row > m_u64_m)
   {
      throw std::runtime_error("Error. Invalid row number!\n");
   }

   return m_pu64_indptr[u64_row+1] - m_pu64_indptr[u64_row];
}

std::vector<uint64_t> OCSR::get_row_col_idxs(uint64_t const & u64_row)
{
   std::vector<uint64_t> u64_res_vec;
   uint64_t u64_row_nnz = this->get_row_nnz(u64_row);
   if (u64_row_nnz != 0)
   {
      u64_res_vec.resize(u64_row_nnz);
      for (uint64_t i = 0U; i < u64_row_nnz; i++)
      {
         u64_res_vec[i] = m_pu64_c_indices[m_pu64_indptr[u64_row]+i];
      }
   }

   return u64_res_vec;
}

std::span<uint64_t> OCSR::get_row_col_idxs_fast(uint64_t const & u64_row) const
{
   std::span<uint64_t> u64_res_sp;
   uint64_t u64_row_nnz = this->get_row_nnz(u64_row);
   if (u64_row_nnz != 0)
   {
      u64_res_sp = std::span<uint64_t>(&m_pu64_c_indices[m_pu64_indptr[u64_row]], u64_row_nnz);
   }

   return u64_res_sp;
}

void OCSR::add_col_idx(uint64_t const & u64_row_idx, uint64_t const & u64_col_idx)
{
   if (u64_row_idx >= this->m_u64_m)
   {
      throw std::runtime_error("Error. Row index out-of-bounds...");
   }

   this->add_col_idx_entry(u64_row_idx, u64_col_idx);
}