/***********************************************************************************************************************
 * @file    DisjointSet.cpp
 * @author  Imanol Etxezarreta (ietxezarretam@gmail.com)
 * 
 * @brief   C++ implementation of the disjoint-set DSA. Taken from 
 *          https://www.geeksforgeeks.org/introduction-to-disjoint-set-data-structure-or-union-find-algorithm/ and may
 *          be subject to some changes depending the necessities of bpotf.
 * 
 * @version 0.2
 * @date    08/05/2024
 * 
 * @copyright Copyright (c) 2024
 * 
 **********************************************************************************************************************/

#include "DisjointSet.h"

DisjSet::DisjSet(long int const & n) 
{
   rank = new int[n];
   parent = new long int[n];
   this->n = n;
   make_set();
}

// Creates n single item sets
void DisjSet::make_set()
{
   for (long int i = 0; i < n; i++) {
      parent[i] = i;
      rank[i] = 0;
   }
}

DisjSet::~DisjSet(void)
{
   delete [] rank;
   delete [] parent;
}

// Finds set of given item x
long int DisjSet::find(long int const & x)
{
   // Finds the representative of the set
   // that x is an element of
   if (parent[x] != x) {

      // if x is not the parent of itself
      // Then x is not the representative of
      // his set,
      parent[x] = find(parent[x]);

      // so we recursively call Find on its parent
      // and move i's node directly under the
      // representative of this set
   }

   return parent[x];
}

void DisjSet::set_parent(long int const & x, long int const & new_parent)
{
   parent[x] = new_parent;
}

int DisjSet::get_rank(long int const & x)
{
   return rank[x];
}

void DisjSet::increase_rank(long int const & x)
{
   rank[x] += 1;
}

// Do union of two sets by rank represented
// by x and y.
int DisjSet::set_union(long int const & x, long int const & y)
{
   // Find current sets of x and y
   long int xset = find(x);
   long int yset = find(y);

   // If they are already in same set
   if (xset == yset)
      return -1;

   // Put smaller ranked item under
   // bigger ranked item if ranks are
   // different
   if (rank[xset] < rank[yset]) 
   {
      parent[xset] = yset;
   }
   else if (rank[xset] > rank[yset]) 
   {
      parent[yset] = xset;
   }
   else 
   {
      // If ranks are same, then increment
      // rank.
      parent[yset] = xset;
      rank[xset] = rank[xset] + 1;
   }

   return 0;
}

long int DisjSet::set_union_opt(long int const & root_set, 
                                 long int const & candidate)
{
   long int ret_val = -1;
   long int candidate_root = find(candidate);

   if (root_set == candidate_root)
      return ret_val;

   if (rank[root_set] < rank[candidate_root]) 
   {
      parent[root_set] = candidate_root;
      ret_val = candidate_root;
   }
   else if (rank[root_set] > rank[candidate_root])
   {
      parent[candidate_root] = root_set;
      ret_val = root_set;
   }
   else 
   {
      parent[candidate_root] = root_set;
      rank[root_set] = rank[root_set] + 1;
      ret_val = root_set;
   }

   return ret_val;
}