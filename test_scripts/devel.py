import BPOTF
import numpy as np
from scipy import sparse

def generate_random_matrix():
    N = np.random.randint(3, 10)  # Random number of rows between 1 and 10
    M = np.random.randint(3, 10)  # Random number of columns between 1 and 10

    # Generate a random NxM matrix with 0s and 1s
    matrix = np.random.randint(2, size=(N, M), dtype=np.uint8)
    return matrix

# Example usage:
mat = generate_random_matrix()
scipy_mat = sparse.csc_matrix(mat)

n, m = np.shape(mat)
print("Matrix shape in python is: {}x{}".format(n, m))
print(mat)
n, m = scipy_mat.get_shape()
print("Matrix shape in SCIPY python is: {}x{}".format(n, m))
print(scipy_mat)

#bpbp.koh_v2(mat)

print(type(mat))
bpotf_obj = BPOTF.OBPOTF(mat, 0.01)
print("first created")
print(type(scipy_mat))
bpotf_obj_scipy = BPOTF.OBPOTF(scipy_mat, 0.01)

bpotf_obj.print_object()
bpotf_obj_scipy.print_object()

fake_syndrome = np.random.randint(2, size=(n, 1), dtype=np.uint8)
print(fake_syndrome)
time = 0.0
error = bpotf_obj.decode(fake_syndrome)
error_scipy = bpotf_obj_scipy.decode(fake_syndrome)

assert np.allclose(error, error_scipy), "Errors not equal!"

print("TEST OK")


bpotf_obj_cln = BPOTF.OBPOTF(mat, 0.01, BPOTF.ECodeType.E_CLN)
bpotf_obj_cln.decode(fake_syndrome)