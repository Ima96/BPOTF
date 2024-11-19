/***********************************************************************************************************************
 * @file    OCSC.cpp
 * @author  Imanol Etxezarreta (ietxezarretam@gmail.com)
 * 
 * @brief   
 * 
 * @version 0.1
 * @date    28/05/2024
 * 
 * @copyright Copyright (c) 2024
 **********************************************************************************************************************/

#include <iostream>
#include <cstring>
#include <algorithm>

#include "OCSC.h"

/***********************************************************************************************************************
 * PRIVATE FUNCTIONS
 **********************************************************************************************************************/
void OCSC::add_row_idx_entry(uint64_t const & u64_row_idx, uint64_t const & u64_col_idx)
{
   std::vector<uint64_t> au64_col_checks = this->get_col_row_idxs(u64_col_idx);
   // Check if already exist
   if (std::find(au64_col_checks.begin(), au64_col_checks.end(), u64_row_idx) != au64_col_checks.end())
   {
      return;
   }

   uint64_t * pu64_temp_buffer = new uint64_t[m_u64_nnz+1];

   uint16_t u16_count = 0;
   for (uint16_t u16_idx = 0, u16_end_cond = au64_col_checks.size(); u16_idx < u16_end_cond; ++u16_idx)
   {
      if (au64_col_checks[u16_idx] < u64_row_idx)
      {
         ++u16_count;
      }
      else
      {
         break;
      }
   }

   uint64_t u64_num_entries_until_insert = m_pu64_indptr[u64_col_idx] + u16_count;
   std::memcpy(pu64_temp_buffer, m_pu64_r_indices, u64_num_entries_until_insert * sizeof(uint64_t));
   pu64_temp_buffer[u64_num_entries_until_insert] = u64_row_idx;
   uint64_t u64_rest = m_u64_nnz-u64_num_entries_until_insert;
   std::memcpy(&pu64_temp_buffer[u64_num_entries_until_insert]+1, 
                  &m_pu64_r_indices[u64_num_entries_until_insert], 
                  u64_rest * sizeof(uint64_t));
   
   delete [] m_pu64_r_indices;
   m_pu64_r_indices = pu64_temp_buffer;

   for (uint64_t u64_idx = u64_col_idx + 1; u64_idx < m_u64_n + 1; ++u64_idx)
   {
      ++m_pu64_indptr[u64_idx];
   }

   ++m_u64_nnz;
   if (u64_row_idx+1 > m_u64_m)
   {
      m_u64_m = u64_row_idx+1;
   }
}


/***********************************************************************************************************************
 * PUBLIC FUNCTIONS
 **********************************************************************************************************************/
OCSC::OCSC(OCSC const & csc_mat):
   m_u64_nnz(csc_mat.m_u64_nnz),
   m_u64_m(csc_mat.m_u64_m),
   m_u64_n(csc_mat.m_u64_n)
{
   m_pu64_indptr = new uint64_t[m_u64_n+1];
   m_pu64_r_indices = new uint64_t[m_u64_nnz];

   std::memcpy(m_pu64_indptr, csc_mat.m_pu64_indptr, (m_u64_n+1) * sizeof(uint64_t));
   std::memcpy(m_pu64_r_indices, csc_mat.m_pu64_r_indices, m_u64_nnz * sizeof(uint64_t));
}

OCSC::OCSC(std::vector<uint64_t> const & u64_nnz_cols_vec,
            std::vector<uint64_t> const & u64_row_idxs_vec,
            uint64_t const & u64_nnz):
            m_u64_nnz(u64_nnz)
{
   m_u64_n = u64_nnz_cols_vec.size() - 1;
   // I think it is not real if the last column is all 0s but al least is something
   m_u64_m = *std::max_element(u64_row_idxs_vec.begin(), u64_row_idxs_vec.end()) + 1;

   m_pu64_indptr= new uint64_t[m_u64_n+1];
   m_pu64_r_indices = new uint64_t[m_u64_nnz];

   std::memcpy(m_pu64_indptr, u64_nnz_cols_vec.data(), (m_u64_n+1) * sizeof(uint64_t));
   std::memcpy(m_pu64_r_indices, u64_row_idxs_vec.data(), m_u64_nnz * sizeof(uint64_t));
}

OCSC::OCSC(std::vector<std::vector<uint8_t>> const & pcm)
{
   uint64_t const u64_m = pcm.size(); // rows
   uint64_t const u64_n = (u64_m > 0) ? pcm[0].size() : 0U; // cols

   m_u64_m = u64_m;
   m_u64_n = u64_n;
   m_u64_nnz = 0U;

   m_pu64_indptr = new uint64_t[u64_n+1];
   m_pu64_indptr[0] = 0;

   // count number of non-zero elements
   for (uint64_t j = 0U; j < u64_n; j++)
   {
      for (uint64_t i = 0U; i < u64_m; i++)
      {
         if (pcm[i][j] != 0U)
         {
            m_u64_nnz++;
         }
      }
      m_pu64_indptr[j+1] = m_u64_nnz;
   }

   // allocate memory
   m_pu64_r_indices = new uint64_t[m_u64_nnz];

   // fill array
   uint64_t u64_cont = 0U;
   for (uint64_t j = 0U; j < u64_n; j++)
   {
      for (uint64_t i = 0U; i < u64_m; i++)
      {
         if (pcm[i][j] != 0U)
         {
            m_pu64_r_indices[u64_cont] = i;
            u64_cont++;
         }
      }
   }
}

OCSC::OCSC(std::vector<uint8_t> const & pcm, uint64_t const & u64_row_num)
{
   if (pcm.size() % u64_row_num != 0)
   {
      throw std::runtime_error("Error. Number of rows must be equal in all columns.");
   }

   m_u64_m = u64_row_num;
   m_u64_n = pcm.size() / u64_row_num;
   m_u64_nnz = 0U;

   m_pu64_indptr = new uint64_t[m_u64_n+1];
   m_pu64_indptr[0] = 0;

   // count number of non-zero elements
   for (uint64_t j = 0U; j < m_u64_n; j++)
   {
      for (uint64_t i = 0U; i < m_u64_m; i++)
      {
         if (pcm[(j*m_u64_m)+i] != 0U)
         {
            m_u64_nnz++;
         }
      }
      m_pu64_indptr[j+1] = m_u64_nnz;
   }

   // allocate memory
   m_pu64_r_indices = new uint64_t[m_u64_nnz];

   // fill array
   uint64_t u64_cont = 0U;
   for (uint64_t j = 0U; j < m_u64_n; j++)
   {
      for (uint64_t i = 0U; i < m_u64_m; i++)
      {
         if (pcm[(j*m_u64_m)+i] != 0U)
         {
            m_pu64_r_indices[u64_cont] = i;
            u64_cont++;
         }
      }
   }
}

OCSC::OCSC(std::span<uint8_t> const & pcm, uint64_t const & u64_row_num)
{
   if (pcm.size() % u64_row_num != 0)
   {
      throw std::runtime_error("Error. Number of rows must be equal in all columns.");
   }

   m_u64_m = u64_row_num;
   m_u64_n = pcm.size() / u64_row_num;
   m_u64_nnz = 0U;

   m_pu64_indptr = new uint64_t[m_u64_n+1];
   m_pu64_indptr[0] = 0;

   // count number of non-zero elements
   for (uint64_t j = 0U; j < m_u64_n; j++)
   {
      for (uint64_t i = 0U; i < m_u64_m; i++)
      {
         if (pcm[(j*m_u64_m)+i] != 0U)
         {
            m_u64_nnz++;
         }
      }
      m_pu64_indptr[j+1] = m_u64_nnz;
   }

   // allocate memory
   m_pu64_r_indices = new uint64_t[m_u64_nnz];

   // fill array
   uint64_t u64_cont = 0U;
   for (uint64_t j = 0U; j < m_u64_n; j++)
   {
      for (uint64_t i = 0U; i < m_u64_m; i++)
      {
         if (pcm[(j*m_u64_m)+i] != 0U)
         {
            m_pu64_r_indices[u64_cont] = i;
            u64_cont++;
         }
      }
   }
}

OCSC::~OCSC()
{
   delete [] m_pu64_r_indices;
   delete [] m_pu64_indptr;
}

void OCSC::print_csc(void)
{
   std::cout << "Shape (MxN): " << m_u64_m << "x" << m_u64_n << std::endl;
   std::cout << "m_u64_nnz: " << m_u64_nnz << std::endl;
   std::cout << "m_pu64_r_indices: [";
   for (uint64_t i = 0; i < m_u64_nnz; i++)
   {
      std::string str = ", ";
      if (i == m_u64_nnz-1)
         str = "";
      std::cout << m_pu64_r_indices[i] << str;
   }
   std::cout << "]\n";

   std::cout << "m_pu64_indptr: [";
   for (uint64_t i = 0; i < m_u64_n+1; i++)
   {
      std::string str = ", ";
      if (i == m_u64_n)
         str = "";
      std::cout << m_pu64_indptr[i] << str;
   }
   std::cout << "]\n";

}

std::vector<std::vector<uint8_t>> OCSC::expand_to_mat(void) const
{
   std::vector<std::vector<uint8_t>> res_mat(m_u64_m, std::vector<uint8_t>(m_u64_n, 0U));

   for (uint64_t i = 1U; i < m_u64_n+1; i++)
   {
      for (uint64_t j = m_pu64_indptr[i-1]; j < m_pu64_indptr[i]; j++)
      {
         uint64_t row_idx = m_pu64_r_indices[j];
         res_mat[row_idx][i-1] = 1U;
      }
   }

   return res_mat;
}

std::vector<uint8_t> OCSC::expand_to_column_major(void)
{
   std::vector<uint8_t> res_vec(m_u64_m * m_u64_n, 0U);

   for (uint64_t i = 1U; i < m_u64_n+1; i++)
   {
      for (uint64_t j = m_pu64_indptr[i-1]; j < m_pu64_indptr[i]; j++)
      {
         uint64_t row_idx = m_pu64_r_indices[j];
         res_vec[((i-1)*m_u64_m)+row_idx] = 1U;
      }
   }

   return res_vec;
}

std::vector<uint8_t> OCSC::expand_to_row_major(void)
{
   std::vector<uint8_t> res_vec;

   std::vector<std::vector<uint8_t>> ppu8_exp_mat = this->expand_to_mat();

   for (uint64_t u64_i = 0U; u64_i < m_u64_m; ++u64_i)
   {
      res_vec.insert(res_vec.end(), ppu8_exp_mat[u64_i].begin(), ppu8_exp_mat[u64_i].end());
   }

   return res_vec;
}

uint64_t OCSC::get_col_nnz(uint64_t const & u64_col) const
{
   if (u64_col > m_u64_n)
   {
      throw std::runtime_error("Error. Invalid column number!\n");
   }

   return m_pu64_indptr[u64_col+1] - m_pu64_indptr[u64_col];
}

std::vector<uint64_t> OCSC::get_col_row_idxs(uint64_t const & u64_col) const
{
   std::vector<uint64_t> u64_res_vec;
   uint64_t u64_col_nnz = this->get_col_nnz(u64_col);
   if (u64_col_nnz != 0)
   {
      u64_res_vec.resize(u64_col_nnz);
      for (uint64_t i = 0U; i < u64_col_nnz; i++)
      {
         u64_res_vec[i] = m_pu64_r_indices[m_pu64_indptr[u64_col]+i];
      }
   }

   return u64_res_vec;
}

std::span<uint64_t> OCSC::get_col_row_idxs_fast(uint64_t const & u64_col) const
{
   std::span<uint64_t> u64_res_sp;
   uint64_t u64_col_nnz = this->get_col_nnz(u64_col);
   if (u64_col_nnz != 0)
   {
      u64_res_sp = std::span<uint64_t>(&m_pu64_r_indices[m_pu64_indptr[u64_col]], u64_col_nnz);
   }

   return u64_res_sp;
}

void OCSC::add_row_idx(uint64_t const & u64_row_idx, uint64_t const & u64_col_idx)
{
   if (u64_col_idx >= this->m_u64_n)
   {
      throw std::runtime_error("Error. Column index out-of-bounds...");
   }

   this->add_row_idx_entry(u64_row_idx, u64_col_idx);
}