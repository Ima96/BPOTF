import stim
from ldpc import bposd_decoder
from beliefmatching import detector_error_model_to_check_matrices
import numpy as np
from tqdm import tqdm

def preprocessing(demcln : stim.DetectorErrorModel, demphen : stim.DetectorErrorModel):
    """
    This function should input two dems. The first one is with circuit-level noise and the second one with phenomenological noise.
    
    1. We compute both parity check matrices with "detector_error_model_to_check_matrices"
    
    2. We iterate over the dem from the cln and decode each column with the phenomenological pcm.
    
    3. Each result will be inputted as a column in the new matrix A, which will serve as a mapping from Hcln to Hphenom.
    
    Return A.
    
    """
    DEMCLN = detector_error_model_to_check_matrices(demcln)
    DEMphen = detector_error_model_to_check_matrices(demphen)
    
    pcm_cln = DEMCLN.check_matrix.toarray()
    
    bp_osd = bposd_decoder(
        DEMphen.check_matrix,
        error_rate = .1,
        osd_method = "osd0"
    )
    
    A = np.zeros((DEMphen.check_matrix.shape[1], DEMCLN.check_matrix.shape[1]), dtype = bool)
    
    for column in tqdm(range(pcm_cln.shape[1])):
        recovered_error = bp_osd.decode(pcm_cln[:,column])
        A[np.where(recovered_error == 1)[0],column] = True
        
    return A


if __name__ == "__main__":
    from quasicyclic import quasi_cyclic_from_params, QuasiCyclic
    from quasicyclic.quasi_cyclic_phenom import quasi_cyclic_from_params_phenom, QuasiCyclic_phenom

    l = 6
    
    qc: QuasiCyclic = quasi_cyclic_from_params(l=l, m=6, A_poly="x^3 + y + y^2", B_poly="y^3 + x + x^2")
    circuit: stim.Circuit = qc.generate_circuit(measure_basis="Z", num_rounds=10, p=0.001)
    
    qc_2: QuasiCyclic_phenom = quasi_cyclic_from_params_phenom(l=l, m=6, A_poly="x^3 + y + y^2", B_poly="y^3 + x + x^2")
    circuit_2: stim.Circuit = qc_2.generate_circuit(measure_basis="Z", num_rounds=10, p=0.001)
    
    demcln = circuit.detector_error_model(decompose_errors=False)
    demphen = circuit_2.detector_error_model(decompose_errors=False)
    
    A = preprocessing(demcln , demphen )
    
    pass