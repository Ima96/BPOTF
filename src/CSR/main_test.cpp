
#include <vector>
#include <iostream>
#include "OCSR.h"

using namespace std;

void printMat(vector<vector<uint8_t>> const & mat)
{
   for (uint64_t i =0U; i < mat.size(); i++)
   {
      for (uint64_t j=0U; j < mat[0].size(); j++)
         cout << int(mat[i][j]) << " ";
      cout << endl;
   }
   cout << endl;
}

int main(void)
{
   vector<vector<uint8_t>> matrix = {
                                       { 0, 0, 0, 0, 1 },
                                       { 5, 8, 0, 0, 0 },
                                       { 0, 0, 3, 0, 0 },
                                       { 0, 6, 0, 0, 1 },
                                       { 0, 0, 0, 7, 0 },
                                       { 4, 0, 1, 0, 0 }
                                    };

   vector<uint8_t> u8_cm_mat_vec = {
                                       0, 5, 0, 0, 0, 4,
                                       0, 8, 0, 6, 0, 0,
                                       0, 0, 3, 0, 0, 1,
                                       0, 0, 0, 0, 7, 0,
                                       1, 0, 0, 1, 0, 0
                                    };

   vector<vector<uint8_t>> u8_expanded_mat;

   OCSR o_csr_mat(matrix);
   cout << "Constructed from Matrix:\n";
   o_csr_mat.print_csr();
   //u8_expanded_mat = o_csr_mat.expand_to_mat();
   printMat(o_csr_mat.expand_to_mat());
   cout << "Number of nnz in row 1: " << o_csr_mat.get_row_nnz(1) << endl;
   vector<uint64_t> col1_nnz_idxs = o_csr_mat.get_row_col_idxs(1);
   cout << "Indeces of nnz in row 1: ";
   for (int i = 0, end = col1_nnz_idxs.size(); i < end; i++)
      cout << int(col1_nnz_idxs[i]) << " ";
   cout << endl << endl;

   cout << "Inserting tests:" << endl;
   o_csr_mat.add_col_idx(2, 5);
   cout << " - Added an entry\n";
   printMat(o_csr_mat.expand_to_mat());
   o_csr_mat.add_col_idx(3, 5);
   cout << " - Added an entry\n";
   printMat(o_csr_mat.expand_to_mat());
   o_csr_mat.add_col_idx(3, 2);
   cout << " - Added an entry\n";
   printMat(o_csr_mat.expand_to_mat());
   o_csr_mat.add_col_idx(3, 0);
   cout << " - Added an entry\n";
   printMat(o_csr_mat.expand_to_mat());
   o_csr_mat.add_col_idx(1, 1);
   cout << " - Added an entry\n";
   printMat(o_csr_mat.expand_to_mat());
   o_csr_mat.print_csr();

   OCSR o_csr_cm_mat(u8_cm_mat_vec, 5);
   cout << "Constructed from Column-Major:\n";
   o_csr_cm_mat.print_csr();

   std::span<uint8_t> sp_u8_cm_mat_vec(u8_cm_mat_vec.data(), u8_cm_mat_vec.size());
   OCSR o_sp_csr_cm_mat(sp_u8_cm_mat_vec, 6);
   cout << "Constructed from Column-Major SPAN:\n";
   o_sp_csr_cm_mat.print_csr();
   std::vector<uint8_t> vec_u8_exp_rm = o_sp_csr_cm_mat.expand_to_row_major();
   std::cout << "Expanded to row major: ";
   for (int i = 0; i < vec_u8_exp_rm.size(); ++i)
   {
      std::cout << int(vec_u8_exp_rm[i]) << " ";
   }
   std::cout << std::endl;
   
   std::vector<uint8_t> vec_u8_exp_cm = o_sp_csr_cm_mat.expand_to_column_major();
   std::cout << "Expanded to col major: ";
   for (int i = 0; i < vec_u8_exp_cm.size(); ++i)
   {
      std::cout << int(vec_u8_exp_cm[i]) << " ";
   }
   std::cout << std::endl;

   OCSR o_csr_mat_cp(o_csr_mat);
   cout << "Constructed from Copy constructor:\n";
   o_csr_mat_cp.print_csr();

   vector<uint64_t> indptr = {0, 1, 3, 4, 5, 7};
   vector<uint64_t> indeces = {1, 1, 3, 2, 4, 0, 3};
   uint64_t nnz = 7;

   OCSR o_csr_data(indeces, indptr, nnz);
   cout << "Constructed from csr data:\n";
   o_csr_data.print_csr();
   printMat(o_csr_data.expand_to_mat());

   return 0;
}