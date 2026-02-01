import scipy.io as sio
from scipy.sparse import csc_matrix
# from BPOTFog import UFCLN as UFCLN2
from SlidingWindowDecoder.src.build_circuit import build_circuit
from SlidingWindowDecoder.src.codes_q import create_bivariate_bicycle_codes, create_circulant_matrix
import numpy as np
# from ldpc import bposd_decoder
from beliefmatching import detector_error_model_to_check_matrices
# from ldpc import bp_decoder

from packaging.version import Version
from ldpc import __version__ as ldpc_version
ldpc_v2 = Version(ldpc_version) >= Version("2.0.0")
print("Using LDPC version v{}".format(ldpc_version))
if ldpc_v2 is True:
    from ldpc import BpDecoder as bp_decoder
    from ldpc import BpOsdDecoder as bposd_decoder
else:
    from ldpc import bp_decoder
    from ldpc import bposd_decoder
import BPOTF
from timeit import default_timer as timer

import os
import psutil

from BPOTF import __version__ as bpotf_version
print(f"BPOTF version is v{bpotf_version}")
import stim, sinter
import numpy as np
# from BPOTF import SinterBpOtfDecoder # TODO ->  SinterBpOtfDecoder not successfully loaded from BPOTF
from sinter_bpotf import SinterBpOtfDecoder
from BPOTF import NoiseType

#TODO -> Remove unnecessary imports.
#TODO -> This test should be done with surface codes.

def main():
    PRINTING = False



    BB_TYPE =  72

    ps = [1e-3, 1.5e-3, 2e-3, 3e-3][::-1]
    NMC = 10**4
    p = ps[0]

    if BB_TYPE == 72:
        # [72, 12, 6] último número es el numero de rondas
        code, A_list, B_list = create_bivariate_bicycle_codes(6, 6, [3], [1,2], [1,2], [3])
        d = 6
        transfer_mat = sio.loadmat('transfermatrices/transferMatrixcodel6m6.mat')['transfMat']
    elif BB_TYPE == 108:
        # [108, 8, 10]
        code, A_list, B_list = create_bivariate_bicycle_codes(9, 6, [3], [1,2], [1,2], [3])
        d = 10
        transfer_mat = sio.loadmat('transfermatrices/transferMatrixcodel9m6.mat')['transfMat']
    elif BB_TYPE == 144:
        # [144, 12, 12]
        code, A_list, B_list = create_bivariate_bicycle_codes(12, 6, [3], [1,2], [1,2], [3])
        d = 12
        transfer_mat = sio.loadmat('transfermatrices/transferMatrixcodel12m6.mat')['transfMat']
    elif BB_TYPE == 288:
        #  [288, 12, 18]
        code, A_list, B_list = create_bivariate_bicycle_codes(12, 12, [3], [2,7], [1,2], [3])
        d = 18
        # transfer_mat = sio.loadmat('transfermatrices/BB288TransfDemsObs.mat')['transfMat']
        transfer_mat = csc_matrix(sio.loadmat('transfermatrices/BB288CSC.mat')['transfMat'])
    else:
        raise Exception("No such option!")

    ps = 3e-3
    NMC = 10**4




    circuit = build_circuit(code, A_list, B_list, 
                            p=ps, # physical error rate
                            num_repeat=d, # usually set to code distance
                            z_basis=True,   # whether in the z-basis or x-basis
                            use_both=False, # whether use measurement results in both basis to decode one basis
                            )
    dem = circuit.detector_error_model()
    bm = detector_error_model_to_check_matrices(dem, True)
    sampler = circuit.compile_detector_sampler()
    # myDecoder = UFCLN(dem, d=d)

    dem = circuit.detector_error_model()
    bm = detector_error_model_to_check_matrices(dem, True)
    sampler = circuit.compile_detector_sampler()
    # myDecoder = UFCLN(dem, d=d)
    # bpbpotf_cpp = BPOTF.OBPOTF(dem, p, BPOTF.OBPOTF.NoiseType.E_CLN, transfer_mat.astype('uint8'))
    dem_data = BPOTF.DemData()
    dem_data.priors = bm.priors
    dem_data.obs_matrix = bm.observables_matrix.toarray('F').astype(np.uint8)
    # Make the phenomenological matrix
    nonzero_counts = bm.check_matrix.sum(axis=0).A1
    selected_cols = np.where(nonzero_counts <= 3)[0]
    temp = bm.observables_matrix[:, selected_cols]
    dem_data.phen_obs_matrix = temp.toarray('F').astype(np.uint8)
    # Make the observables phenomenological matrix
    temp2 = bm.check_matrix[:, selected_cols]
    dem_data.phen_check_matrix = temp2.toarray('F').astype(np.uint8)
    dem_data.transfer_matrix = transfer_mat.astype(np.uint8)

    interval_value = d  # You can change this value to your desired interval
    max_iters_dem = 100
    max_iters_sparse = 400
    max_iters_forest = 100

    bp_iterations = [max_iters_dem, max_iters_sparse, max_iters_forest]


    bpbp_otf_decoder = BPOTF.OBPOTF(
        bm.check_matrix, 
        p, 
        BPOTF.NoiseType.E_CLN,
        ps_ext_dem_data=dem_data,
        po_ext_bp_iters = bp_iterations
    )
        



    # TODO ADAPT sinter to next



    task = sinter.Task(
        circuit=circuit,
        collection_options=sinter.CollectionOptions(max_shots=10, max_errors=10),
    )

    # iterations: 35, default and 800
    sinter.collect(
        num_workers=1,
        tasks=[task],
        decoders=["bpotf"],
        custom_decoders={
            "bpotf": SinterBpOtfDecoder(
                # p=0.002, # Esto deberían ser las priors, no un float. Creo que sinter debería ser solo para cln. 
                # noise_type= NoiseType.E_CC,
                # po_ext_bp_iters=bp_iterations,
                bm.check_matrix, 
                p, 
                BPOTF.NoiseType.E_CLN,
                ps_ext_dem_data=dem_data,
                po_ext_bp_iters = bp_iterations
            )
            },
        print_progress=True,
        save_resume_filepath="sinter_results/results.csv"
    )

# This block seems mandatory on Windows when using siner/multiprocessing
if __name__ == "__main__":
    main()