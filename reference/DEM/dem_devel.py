import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

import scipy.io as sio
conts = sio.loadmat('testBBCLNmap144_12_12_12rounds_p_001.mat')
from SlidingWindowDecoder.src.build_circuit import build_circuit
from SlidingWindowDecoder.src.codes_q import create_bivariate_bicycle_codes, create_circulant_matrix
import BPOTF
from BPOTF2 import UFCLN
import numpy as np
from stimbposd import detector_error_model_to_check_matrices
from timeit import default_timer as timer
from ldpc import bp_decoder

def checkDifferences(mat1, mat2):
    if mat1.shape[0] != mat2.shape[0] or mat1.shape[1] != mat2.shape[1]:
        print("different shapes!")
        return

    for i in range(mat1.shape[0]):
        if mat1[i] != mat2[i]:
            print(f"Difference in idx({i}): {mat1[i]} vs {mat2[i]}")


dem = conts['dem']
priors = conts['priors']
transf_M = conts['transfMatFull']
H_phen = conts['Hphen']
obs = conts['obs']
hz = conts['hz']
code, A_list, B_list = create_bivariate_bicycle_codes(12, 6, [3], [1,2], [1,2], [3])
d = 12
p = 0.003
NMC = 10**4

circuit = build_circuit(code, A_list, B_list, 
                        p=p, # physical error rate
                        num_repeat=d, # usually set to code distance
                        z_basis=True,   # whether in the z-basis or x-basis
                        use_both=False, # whether use measurement results in both basis to decode one basis
                        )
dem = circuit.detector_error_model()

cm = detector_error_model_to_check_matrices(dem, True)

transfer_mat = sio.loadmat('transfermatrices/transferMatrixcodel12m6.mat')['transfMat']
transfer_mat_T = transfer_mat.transpose()
print(transfer_mat.shape)

start_bpotf = timer()
bpotf = BPOTF.OBPOTF(dem, 0.01, BPOTF.NoiseType.CLN, transfer_mat.astype('uint8'))
stop_bpotf = timer()

start_ton = timer()
bpotf_ton = UFCLN(dem, d=12)
stop_ton = timer()

print(f"Elapsed time creating C++ object: {stop_bpotf-start_bpotf}")
print(f"Elapsed time creating UFCLN object: {stop_ton-start_ton}")

bpotf_H = bpotf.get_pcm().reshape(bpotf.get_rows(), bpotf.get_cols())
if not np.all(bpotf_ton.H == bpotf_H):
    print("PCMS are not equal...")
else:
    print("PCM SUCCESS")

bpotf_H_phen = bpotf.get_phen_pcm().reshape(bpotf.get_phen_pcm_rows(), bpotf.get_phen_pcm_cols())
if not np.array_equal(bpotf_ton.H_phen, bpotf_H_phen):
    print("PHEN PCMS are not equal...")
    print(bpotf_ton.H_phen.shape)
    print(bpotf_H_phen.shape)
else:
    print("PHEN PCM SUCCESS")

bpotf_obs = bpotf.get_obs().reshape(bpotf.get_obs_rows(), bpotf.get_obs_cols())
if not np.array_equal(bpotf_ton.obs, bpotf_obs):
    print("Observable matrices are not equal...")
else:
    print("OBS SUCCESS")

bpotf_transf = bpotf.get_transf().reshape(bpotf.get_transf_rows(), bpotf.get_transf_cols())
if not np.array_equal(transfer_mat.astype('uint8'), bpotf_transf):
    print("Transfer matrices are not equal...")
else:
    print("Transfer SUCCESS")

print(cm.priors)
print(bpotf_ton.priors)
print(bpotf.get_priors())