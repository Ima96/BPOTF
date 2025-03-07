import scipy.io as sio
from scipy.sparse import csc_matrix
conts = sio.loadmat('testBBCLNmap144_12_12_12rounds_p_001.mat')
from BPOTF2 import UFCLN
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

PRINTING = False

dem = conts['dem']
priors = conts['priors']
transf_M = conts['transfMatFull']
H_phen = conts['Hphen']
obs = conts['obs']
hz = conts['hz']

BB_TYPE =  72
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

name_file_bpbpotf = f'results/bpbpotf{d}bbcode.txt'
name_file_bposd = f'results/bposd{d}bbcode.txt'
ps = [ 2e-3, 3e-3, 4e-3, 5e-3]
NMC = 10**4


for p in ps:

    number_of_iterations = 0
    circuit = build_circuit(code, A_list, B_list, 
                            p=p, # physical error rate
                            num_repeat=d, # usually set to code distance
                            z_basis=True,   # whether in the z-basis or x-basis
                            use_both=False, # whether use measurement results in both basis to decode one basis
                            )
    dem = circuit.detector_error_model()
    bm = detector_error_model_to_check_matrices(dem, True)
    sampler = circuit.compile_detector_sampler()
    # myDecoder = UFCLN(dem, d=d)
    bpbpotf_cpp = BPOTF.OBPOTF(dem, p, BPOTF.OBPOTF.NoiseType.E_CLN, transfer_mat.astype('uint8'))
    bposd = bposd_decoder(
        bm.check_matrix,
        channel_probs = bm.priors,
        max_iter = 1000,
        osd_method = "osd_0"
    )
    process = psutil.Process(os.getpid())
    print(f"Memory used: {process.memory_info().rss / 1024 ** 2:.2f} MB")  # Memory in MB

    Pl_bposd = 0
    Pl_py_otf = 0
    Pl_cpp_otf = 0


    sum_ton = 0
    sum_cpp = 0
    sum_bposd = 0
    worst_times = [0, 0, 0] # In this order, BPOSD, PyOTF, CppOTF
    # stages_list = []
    updates = 0

    while (Pl_bposd < 100) or (Pl_cpp_otf < 100):
        print(f'\nNumber of iterations {number_of_iterations}')
        number_of_iterations += NMC
        print(f'Plbpbpotf {Pl_cpp_otf}')
        print(f'Plbposd {Pl_bposd}\n')
        
        detection_events, observable_flips = sampler.sample(NMC, separate_observables=True)
        # detection_events = np.concatenate([detection_events, error_case_syndrome])
        # observable_flips = np.concatenate([observable_flips, error_case_observable])
        # mdic = {'demPcm' : bm.check_matrix.toarray(), "label": " demPcm"}
        # sio.savemat("dem_pcm.mat", mdic)
        
        
        for index, detection_event in enumerate(detection_events):
            
            
            
            
            observable_flip = observable_flips[index]

            finished = 100 * (index / NMC)

            bposd_failed = False
            py_otf_failed = False
            cpp_otf_failed = False

            # start_ton = timer()
            # recovered_error, stage = myDecoder.decode(detection_event)#[0]
            # stop_ton = timer()
            # py_otf_time = (stop_ton-start_ton)
            
            # start_cpp = timer()
            recovered_error_cpp = bpbpotf_cpp.decode(detection_event.astype(np.uint8))
            # stop_cpp = timer()
            # cpp_otf_time = (stop_cpp-start_cpp)

            # start_bposd = timer()
            recovered_error2 = (bm.observables_matrix @ bposd.decode(detection_event)) %2
            # end_bposd = timer()
            # bposd_time = (end_bposd-start_bposd)

            # sum_ton += py_otf_time
            # sum_cpp += cpp_otf_time
            # sum_bposd += bposd_time

            # if bposd_time > worst_times[0]:
            #     worst_times[0] = bposd_time
            # # if py_otf_time > worst_times[1]:
            # #     worst_times[1] = py_otf_time
            # if cpp_otf_time > worst_times[2]:
            #     worst_times[2] = cpp_otf_time

            # stages_list.append(stage)
            
            if not np.all(recovered_error2 == observable_flip):
                bposd_failed = True
                Pl_bposd += 1
                # print('BPOSD failed')
            # if not np.all(recovered_error == observable_flip):
            #     py_otf_failed = True
            #     Pl_py_otf += 1
                # print('PyOTF failed')
            if not np.all(recovered_error_cpp == observable_flip):
                cpp_otf_failed = True
                Pl_cpp_otf += 1
                # print('CppOTF failed')

            if ((bposd_failed is True) or (py_otf_failed is True) or (cpp_otf_failed is True)) and PRINTING:
                print(f'Iteration number {index}')
                print(f'Pl BPOSD = {Pl_bposd/index}')
                # print(f'Pl PyOTF = {Pl_py_otf/index}')
                print(f'Pl CppOTF = {Pl_cpp_otf/index}')
                print('\n')

            if divmod(finished, 10) == (updates, 0):
                updates += 1
                print('Finished processing {} % of all events'.format(int(finished)))
                process = psutil.Process(os.getpid())
                print(f"Memory used: {process.memory_info().rss / 1024 ** 2:.2f} MB")  # Memory in MB
            

    # if (1 in stages_list) and (2 in stages_list) and (3 in stages_list):
    #     print("All stages tested!")
    # else:
    #     print("Not all stages tested...")

    pe_bposd = Pl_bposd/number_of_iterations
    ped_bposd = pe_bposd/d
    # pe_py_otf = Pl_py_otf/NMC
    # ped_py_otf = pe_py_otf/d
    pe_cpp_otf = Pl_cpp_otf/number_of_iterations
    ped_cpp_otf = pe_cpp_otf/d

    print(f"Simulation data:")
    print(f" - p = {p}")
    print(f" - NMC = {NMC}")
    print(f" - d = {d}")
    print("")
    print(f"============ Probability results ============")
    print(f"BPOSD results:")
    print(f" * p(e): {pe_bposd}")
    print(f" * p(e)/d: {ped_bposd}")
    # print(f"PyOTF results:")
    # print(f" * p(e): {pe_py_otf}")
    # print(f" * p(e)/d: {ped_py_otf}")
    print(f"CppOTF results:")
    print(f" * p(e): {pe_cpp_otf}")
    print(f" * p(e)/d: {ped_cpp_otf}")
    print("")
    print(f"============ Timing results ============")
    # print(f"BPOSD times:")
    # print(f" * Average time: {sum_bposd/NMC} s")
    # print(f" * Worst time: {worst_times[0]} s")
    # # print(f"PyOTF times:")
    # # print(f" * Average time: {sum_ton/NMC} s")
    # # print(f" * Worst time: {worst_times[1]} s")
    # print(f"CppOTF times:")
    # print(f" * Average time: {sum_cpp/NMC} s")
    # print(f" * Worst time: {worst_times[2]} s")

    # Here I write the results:
    
    with open(name_file_bpbpotf, "a") as file:
        file.write(f"{p}\t\t {ped_cpp_otf} \t\t {pe_cpp_otf} \t\t {number_of_iterations}\n")

    with open(name_file_bposd, "a") as file:
        file.write(f"{p}\t\t {ped_bposd} \t\t {pe_bposd} \t\t {number_of_iterations}\n")
    # TODO Print error rate, logical error rate, logical error rate / d and number of iterations
