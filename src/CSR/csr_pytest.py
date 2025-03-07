import scipy
import scipy.sparse
import numpy as np

matrix = [
            [0, 0, 0, 0, 1],
            [5, 8, 0, 0, 0],
            [0, 0, 3, 0, 0],
            [0, 6, 0, 0, 1],
            [0, 0, 0, 7, 0]
        ]

csr_mat = scipy.sparse.csr_matrix(matrix)

print(csr_mat)
print(csr_mat.count_nonzero())
print(csr_mat.indices)
print(csr_mat.indptr)

data = np.array([1, 1, 1, 1, 1, 1, 1])
indptr = np.array([0, 1, 3, 4, 5, 7])
indices = np.array([1, 1, 3, 2, 4, 0, 3])

csr_2 = scipy.sparse.csr_matrix((data,indices,indptr))
# print(csr_2)
print(csr_2.toarray())