from packaging.version import Version
from ldpc import __version__ as ldpc_version
ldpc_v2 = Version(ldpc_version) >= Version("2.0.0")
import numpy as np
print("Using LDPC version v{}".format(ldpc_version))
if ldpc_v2 is True:
    from ldpc import BpDecoder as bp_decoder
    from ldpc import BpOsdDecoder as bposd_decoder
else:
    from ldpc import bp_decoder
    from ldpc import bposd_decoder
from scipy.sparse import csr_matrix, csc_matrix



class BPBPOSD:
    def __init__(self, transf_mat, H, H_sparse, priors, obs, sobs, max_iter=100, bp_type = "ms", ms_scaling_factor=1):
        self.bp = bp_decoder(
            H,
            channel_probs = priors,
            bp_method = bp_type,
            max_iter = max_iter,
            ms_scaling_factor = ms_scaling_factor
        )
        self.transf_mat = transf_mat
        self.transf_M_red = csr_matrix(self.transf_mat)
        self.H = H
        self.H_sparse = H_sparse
        self.bposd = bposd_decoder(
            H_sparse,
            max_iter = 0,
            channel_probs = self.propagation(priors),
            osd_method = "osd_0"
        )
        self.obs = obs
        self.sobs = sobs
    
    def propagation(self, p):
        p_ph = np.zeros(self.transf_M_red.shape[0])
        for row in range(len(p_ph)):
            columns = self.transf_M_red[row,:].indices
            pphen_prod = 1
            for _, col in enumerate(columns):
                # Exclude the current element from the product calculation
                pphen_prod *= (1-(2*p[col]))
            p_ph[row] = .5*(1-pphen_prod)
        return p_ph
    
    def decode(self, syndrome):
        
        # BP
        recovery = self.bp.decode(syndrome)
        
        if self.bp.converge:
            return (self.obs @ recovery) % 2
        
        # BPOSD
        llrs = self.bp.log_prob_ratios
        llrs[llrs < 1e14]   = 1e14
        llrs[llrs > -1e14] = -1e14
        
        ps = 1/(1+np.exp(llrs))
        ps = np.clip(ps, 1e-14, 1-1e-14)
        
        p_ph = self.propagation(ps)
        
        self.bposd.update_channel_probs(p_ph)
        recovery = self.bposd.decode(syndrome)
        
        return (self.sobs @ recovery ) % 2
    
    
    
if __name__ == "__main__":
    from SlidingWindowDecoder.src.codes_q import create_bivariate_bicycle_codes
    from SlidingWindowDecoder.src.build_circuit import build_circuit
    from beliefmatching import detector_error_model_to_check_matrices
    import scipy.io as sio
    
    
    BB_TYPE = 108


    ps = [1e-3, 1.5e-3, 2e-3, 3e-3][::-1]
    NMC = 10**4
    max_iter = 100
    
    for p in ps:
    
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
            transfer_mat = csc_matrix(sio.loadmat('transfermatrices/BB288CSC.mat')['transfMat']).toarray('F')
        else:
            raise Exception("No such option!")
        
        print(f"Running for p = {p}")
        circuit = build_circuit(code, A_list, B_list, 
                                p=p, # physical error rate
                                num_repeat=d, # usually set to code distance
                                z_basis=True,   # whether in the z-basis or x-basis
                                use_both=False, # whether use measurement results in both basis to decode one basis
                                )
        dem = circuit.detector_error_model()
        bm = detector_error_model_to_check_matrices(dem, True)
        sampler = circuit.compile_detector_sampler()
        
        H = bm.check_matrix.toarray()
        priors = bm.priors
        obs = bm.observables_matrix.toarray('F').astype(np.uint8)
        nonzero_counts = bm.check_matrix.sum(axis=0).A1
        selected_cols = np.where(nonzero_counts <= 3)[0] # Ima aquí pilla las columnas de mínimo peso que son 3
        sobs = bm.observables_matrix[:, selected_cols].toarray('F').astype(np.uint8)
        H_sparse = H[:, selected_cols]
        

        
        BPBPOSD = BPBPOSD(transfer_mat, H, H_sparse, priors, obs, sobs, max_iter=100, bp_type = "ms", ms_scaling_factor=1)
        
        BPOSD =  bposd_decoder(
            H,
            max_iter = 100,
            channel_probs = priors,
            osd_method = "osd_0"
        )
        
        Pl_bposd = 0
        Pl_bpbposd = 0
        number_of_iters = 0
        
        while min([Pl_bposd, Pl_bpbposd]) < 100:
            detection_events, observable_flips = sampler.sample(NMC, separate_observables=True)
            number_of_iters += NMC
            for index, detection_event in enumerate(detection_events):
                observable_flip = observable_flips[index]

                finished = 100 * (index / NMC)

                bposd_failed = False
                bpbposd_failed = False
                
                recovered_error_bposd = (bm.observables_matrix @ BPOSD.decode(detection_event)) %2
        
                if not np.all(recovered_error_bposd == observable_flip):
                    Pl_bposd += 1
        
                recovered_error_bpbposd = BPBPOSD.decode(detection_event)
                if not np.all(recovered_error_bpbposd == observable_flip):
                    Pl_bpbposd += 1
                    
            print(f"Pl_bposd = {Pl_bposd/(number_of_iters*d)}")
            print(f'Pl_bpbposd = {Pl_bpbposd/(number_of_iters*d)}')
            
        # TODO Josu, aquí puedes escribir los resultados en un fichero o algo así
        # print(f"Pl_bposd = {Pl_bposd/(number_of_iters*d)}")
        # print(f'Pl_bpbposd = {Pl_bpbposd/(number_of_iters*d)}')
        