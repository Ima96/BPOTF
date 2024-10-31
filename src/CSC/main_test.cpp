
#include <vector>
#include <iostream>
#include "OCSC.h"

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
                                       { 0, 0, 0, 7, 0 }
                                    };

   vector<uint8_t> u8_cm_mat_vec = {
                                       0, 5, 0, 0, 0,
                                       0, 8, 0, 6, 0,
                                       0, 0, 3, 0, 0,
                                       0, 0, 0, 0, 7,
                                       1, 0, 0, 1, 0
                                    };

   vector<vector<uint8_t>> u8_expanded_mat;

   OCSC o_csc_mat(matrix);
   cout << "Constructed from Matrix:\n";
   o_csc_mat.print_csc();
   //u8_expanded_mat = o_csc_mat.expand_to_mat();
   printMat(o_csc_mat.expand_to_mat());
   cout << "Number of nnz in column 1: " << o_csc_mat.get_col_nnz(1) << endl;
   vector<uint64_t> col1_nnz_idxs = o_csc_mat.get_col_row_idxs(1);
   cout << "Indeces of nnz in column 1: ";
   for (int i = 0, end = col1_nnz_idxs.size(); i < end; i++)
      cout << int(col1_nnz_idxs[i]) << " ";
   cout << endl << endl;

   cout << "Inserting tests:" << endl;
   o_csc_mat.add_row_idx(5, 2);
   cout << " - Added an entry\n";
   printMat(o_csc_mat.expand_to_mat());
   o_csc_mat.add_row_idx(5, 3);
   cout << " - Added an entry\n";
   printMat(o_csc_mat.expand_to_mat());
   o_csc_mat.add_row_idx(2, 3);
   cout << " - Added an entry\n";
   printMat(o_csc_mat.expand_to_mat());
   o_csc_mat.add_row_idx(0, 3);
   cout << " - Added an entry\n";
   printMat(o_csc_mat.expand_to_mat());
   o_csc_mat.add_row_idx(1, 1);
   cout << " - Added an entry\n";
   printMat(o_csc_mat.expand_to_mat());
   o_csc_mat.print_csc();

   OCSC o_csc_cm_mat(u8_cm_mat_vec, 5);
   cout << "Constructed from Column-Major:\n";
   o_csc_cm_mat.print_csc();

   std::span<uint8_t> sp_u8_cm_mat_vec(u8_cm_mat_vec);
   OCSC o_sp_csc_cm_mat(sp_u8_cm_mat_vec, 5);
   cout << "Constructed from Column-Major SPAN:\n";
   o_sp_csc_cm_mat.print_csc();

   OCSC o_csc_mat_cp(o_csc_mat);
   cout << "Constructed from Copy constructor:\n";
   o_csc_mat_cp.print_csc();

   vector<uint64_t> indeces = {0, 1, 3, 4, 5, 7};
   vector<uint64_t> indptr = {1, 1, 3, 2, 4, 0, 3};
   uint64_t nnz = 7;

   OCSC o_csc_data(indeces, indptr, nnz);
   cout << "Constructed from CSC data:\n";
   o_csc_data.print_csc();

   return 0;
}