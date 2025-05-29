from beliefmatching import BeliefMatching, detector_error_model_to_check_matrices
import stim
import numpy as np
import BPOTF
import pymatching
from parser import parser
import stim
import numpy as np
from beliefmatching import BeliefMatching
import parser

NMC = 1000
d = 11
ps = [1e-3, 2.5e-3, 5e-3, 7.5e-3, 1e-2]
bp_iters = 1000
bp_iters_bpbpotf = np.array([bp_iters, bp_iters, 100])

for p in ps:
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        rounds=d,
        distance=d,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p,
        after_clifford_depolarization=p
    )
    
    Pl_bpbpotf = 0
    Pl_pymatching = 0
    Pl_beliefmatching = 0
    
    dem = circuit.detector_error_model(decompose_errors=True)
    bm = detector_error_model_to_check_matrices(dem, allow_undecomposed_hyperedges=False)
    
    
    otf_matrix = parser(circuit, bm.edge_check_matrix)   # Esta es la matriz de paridad con checks virtuales para x y z checks.
    
    dem_data = BPOTF.DemData()
    dem_data.priors = bm.priors
    dem_data.obs_matrix = bm.observables_matrix.toarray('F').astype(np.uint8)


    dem_data.phen_obs_matrix = bm.edge_observables_matrix.toarray('F').astype(np.uint8)
        # Make the observables phenomenological matrix
    dem_data.phen_check_matrix = bm.edge_check_matrix.toarray('F').astype(np.uint8)
    dem_data.transfer_matrix = bm.hyperedge_to_edge_matrix.toarray().astype(np.uint8)


    # Decoders being considered
    bpotf = BPOTF.OBPOTF(
        bm.check_matrix.toarray(), 
        p, 
        BPOTF.NoiseType.E_CLN,
        ps_ext_dem_data=dem_data, 
        po_ext_bp_iters= bp_iters_bpbpotf
        )
    
    bm = BeliefMatching(circuit, max_bp_iters=bp_iters)
    pm = pymatching.Matching.from_detector_error_model(dem)
    
    # A partir de aquí empezaremos a samplear errores.
    sampler = circuit.compile_detector_sampler()
    number_of_iters = 0
    
    while min([Pl_bpbpotf, Pl_pymatching, Pl_beliefmatching]) < 100:
        shots, observables = sampler.sample(NMC, separate_observables=True)
        number_of_iters += NMC
        for index, shot in enumerate(shots):
            corrected = True
            match_recovery = pm.decode(shot)
            bpbp_otf_recovery = bpotf.decode(shot)
            bm_recovery = bm.decode(shot)
            
            if not np.all(match_recovery == observables[index]):
                Pl_pymatching += 1
                corrected = False
            
            if not np.all(bpbp_otf_recovery == observables[index]) or not bpotf.has_converged():
                Pl_bpbpotf += 1
                corrected = False
                
            if not np.all(bm_recovery == observables[index]):
                Pl_beliefmatching += 1
                corrected = False
            
            if not corrected:
                print(f'Number of iters {number_of_iters}')
                print(f'MWPM error numbers {Pl_pymatching}')
                print(f"MWPM Error rate = {Pl_pymatching/number_of_iters}")
                print(f'BPBPOTF error numbers {Pl_bpbpotf}')
                print(f"BPBPOTF Error rate = {Pl_bpbpotf/number_of_iters}")
                print(f'BM error numbers {Pl_beliefmatching}')
                print(f"BM Error rate = {Pl_beliefmatching/number_of_iters}")
                print('\n')
                
    # Save results to separate files for each decoding process
    with open(f"results/bpbpotf_d{d}.txt", "a") as bpbpotf_file:
        bpbpotf_file.write(f"{p}\t\t{Pl_bpbpotf/number_of_iters}\t\t{number_of_iters}\n")
    
    with open(f"results/mwpm_d{d}.txt", "a") as mwpm_file:
        mwpm_file.write(f"{p}\t\t{Pl_pymatching/number_of_iters}\t\t{number_of_iters}\n")
    
    with open(f"results/bm_d{d}.txt", "a") as bm_file:
        bm_file.write(f"{p}\t\t{Pl_beliefmatching/number_of_iters}\t\t{number_of_iters}\n")