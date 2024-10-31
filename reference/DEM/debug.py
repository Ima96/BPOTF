import scipy.io as sio
conts = sio.loadmat('testBBCLNmap144_12_12_12rounds_p_001.mat')
from BPOTF2 import UFCLN
# from BPOTFog import UFCLN as UFCLN2
from SlidingWindowDecoder.src.build_circuit import build_circuit
from SlidingWindowDecoder.src.codes_q import create_bivariate_bicycle_codes, create_circulant_matrix
import numpy as np
from ldpc import bposd_decoder
from beliefmatching import detector_error_model_to_check_matrices
from ldpc import bposd_decoder

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
bm = detector_error_model_to_check_matrices(dem, True)
sampler = circuit.compile_detector_sampler()
myDecoder = UFCLN(dem, d=12)
bposd = bposd_decoder(
    bm.check_matrix,
    channel_probs = bm.priors,
    max_iter = 100,
    osd_method = "osd_0"
)

detection_events, observable_flips = sampler.sample(NMC, separate_observables=True)
Pl = 0
Plog = 0

for index, detection_event in enumerate(detection_events):
    observable_flip = observable_flips[index]
    # print(index)
    recovered_error = myDecoder.decode(detection_event)[0]
    recovered_error2 = (bm.observables_matrix @ bposd.decode(detection_event)) %2
    
    if not np.all(recovered_error2 == recovered_error):
        if not np.all(recovered_error2 == observable_flip):
            Plog += 1
            print('BPOSD failed')
            print(f'Iteration number {index}')
            print(f'Pl BPOSD = {Plog/index}')
            print(f'Pl BPBP = {Pl/index}')
            print('\n')
    if not np.all(recovered_error == observable_flip):
        Pl += 1
        print('BPOTF failed')
        print(f'Iteration number {index}')
        print(f'Pl BPOSD = {Plog/index}')
        print(f'Pl BPBP = {Pl/index}')
        print('\n')
    else:
        if not np.all(recovered_error2 == observable_flip):
            print('Both failed')
            Plog += 1
            Pl += 1
            print(f'Iteration number {index}')
            print(f'Pl BPOSD = {Plog/index}')
            print(f'Pl new = {Pl/index}')
            print('\n')
print(Plog/NMC)
print(Pl/NMC)