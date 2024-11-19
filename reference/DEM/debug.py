import scipy.io as sio
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

from BPOTF import __version__ as bpotf_version
print(f"BPOTF version is v{bpotf_version}")

dem = conts['dem']
priors = conts['priors']
transf_M = conts['transfMatFull']
H_phen = conts['Hphen']
obs = conts['obs']
hz = conts['hz']
code, A_list, B_list = create_bivariate_bicycle_codes(12, 6, [3], [1,2], [1,2], [3])
d = 12
p = 0.0003
NMC = 10**6

circuit = build_circuit(code, A_list, B_list, 
                        p=p, # physical error rate
                        num_repeat=d, # usually set to code distance
                        z_basis=True,   # whether in the z-basis or x-basis
                        use_both=False, # whether use measurement results in both basis to decode one basis
                        )
dem = circuit.detector_error_model()
bm = detector_error_model_to_check_matrices(dem, True)
sampler = circuit.compile_detector_sampler()
myDecoder = UFCLN(dem, d=12)
transfer_mat = sio.loadmat('transfermatrices/transferMatrixcodel12m6.mat')['transfMat']
transfer_mat_T = transfer_mat.transpose()
bpbpotf_cpp = BPOTF.OBPOTF(dem, p, BPOTF.OBPOTF.NoiseType.E_CLN, transfer_mat.astype('uint8'))
bposd = bposd_decoder(
    bm.check_matrix,
    channel_probs = bm.priors,
    max_iter = 100,
    osd_method = "osd_0"
)

detection_events, observable_flips = sampler.sample(NMC, separate_observables=True)
Pl_bposd = 0
Pl_py_otf = 0
Pl_cpp_otf = 0

# detection_events = np.concatenate([detection_events, error_case_syndrome])
# observable_flips = np.concatenate([observable_flips, error_case_observable])

sum_ton = 0
sum_cpp = 0
sum_bposd = 0
worst_times = [0, 0, 0] # In this order, BPOSD, PyOTF, CppOTF
stages_list = []
for index, detection_event in enumerate(detection_events):
    observable_flip = observable_flips[index]

    start_ton = timer()
    recovered_error, stage = myDecoder.decode(detection_event)#[0]
    stop_ton = timer()
    py_otf_time = (stop_ton-start_ton)
    
    start_cpp = timer()
    recovered_error_cpp = bpbpotf_cpp.decode(detection_event.astype(np.uint8))
    stop_cpp = timer()
    cpp_otf_time = (stop_cpp-start_cpp)

    start_bposd = timer()
    recovered_error2 = (bm.observables_matrix @ bposd.decode(detection_event)) %2
    end_bposd = timer()
    bposd_time = (end_bposd-start_bposd)

    sum_ton += py_otf_time
    sum_cpp += cpp_otf_time
    sum_bposd += bposd_time

    if bposd_time > worst_times[0]:
        worst_times[0] = bposd_time
    if py_otf_time > worst_times[1]:
        worst_times[1] = py_otf_time
    if cpp_otf_time > worst_times[2]:
        worst_times[2] = cpp_otf_time

    stages_list.append(stage)

    if not np.all(recovered_error2 == recovered_error) or not np.all(recovered_error2 == recovered_error_cpp):
        if not np.all(recovered_error2 == observable_flip):
            bposd_failed = True
            Pl_bposd += 1
            print('BPOSD failed')
        if not np.all(recovered_error == observable_flip):
            py_otf_failed = True
            Pl_py_otf += 1
            print('PyOTF failed')
        if not np.all(recovered_error_cpp == observable_flip):
            cpp_otf_failed = True
            Pl_cpp_otf += 1
            print('CppOTF failed')

        print(f'Iteration number {index}')
        print(f'Pl BPOSD = {Pl_bposd/index}')
        print(f'Pl PyOTF = {Pl_py_otf/index}')
        print(f'Pl CppOTF = {Pl_cpp_otf/index}')
        print('\n')
        

print("Finished!")
print(f"Final BPOSD: {Pl_bposd/NMC}")
print(f"Final PyOTF: {Pl_py_otf/NMC}")
print(f"Final CppOTF: {Pl_cpp_otf/NMC}")
if (1 in stages_list) and (2 in stages_list) and (3 in stages_list):
    print("All stages tested!")
else:
    print("Not all stages tested...")

print(f"============ Timing results ============")
print(f"BPOSD times:")
print(f" * Average time: {sum_bposd/NMC} s")
print(f" * Worst time: {worst_times[0]} s")
print(f"PyOTF times:")
print(f" * Average time: {sum_ton/NMC} s")
print(f" * Worst time: {worst_times[1]} s")
print(f"CppOTF times:")
print(f" * Average time: {sum_cpp/NMC} s")
print(f" * Worst time: {worst_times[2]} s")