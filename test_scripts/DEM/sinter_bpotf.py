import pathlib
import numpy as np
import sinter
import stim
from BPOTF import OBPOTF, NoiseType, DemData
from beliefmatching import detector_error_model_to_check_matrices

from ldpc.ckt_noise.dem_matrices import detector_error_model_to_check_matrices

class SinterBpOtfDecoder(sinter.Decoder):
    """
    TODO: Insert here a description of what is this for.
    TODO: It inputs too many arguments. Maybe we should place most of them in DemData.
    """

    def __init__(
            self,
            p, # TODO, no veo de que sirve tener la probabilidad aquí.
            dem_data,
            noise_type = NoiseType.E_CLN,
            po_otf_csc_mat = None,
            po_ext_bp_iters = None,
            transfer_matrix = None,
            otf_matrix = None,
            phen_check_matrix = None,
            phen_obs_matrix = None,
            allow_undecomposed_hyperedges = True,
            pcm = None
            ):
        
        self.m_p = p
        self.m_noise_type = noise_type
        self.m_po_otf_csc_mat = po_otf_csc_mat
        self.m_po_ext_bp_iters = po_ext_bp_iters
        self.m_dem_data = dem_data
        # self.dem_data_phen_check_matrix = phen_check_matrix => In DemData
        # self.transfer_matrix = transfer_matrix => In DemData
        self.m_otf_matrix = otf_matrix
        self.m_phen_obs_matrix = phen_obs_matrix
        self.m_allow_undecomposed_hyperedges = allow_undecomposed_hyperedges
        self.bpotf = OBPOTF(
            pcm, 
            self.m_p, 
            self.m_noise_type,
            ps_ext_dem_data=self.m_dem_data, 
            po_ext_bp_iters= self.m_po_ext_bp_iters
        )

    def decode_via_files(
            self,
            *,
            num_shots: int,
            num_dets: int,
            num_obs: int,
            dem_path: pathlib.Path,
            dets_b8_in_path: pathlib.Path,
            obs_predictions_b8_out_path: pathlib.Path,
            tmp_dir: pathlib.Path
            ) -> None:

        self.m_dem = stim.DetectorErrorModel.from_file(dem_path)
        self.m_matrices = detector_error_model_to_check_matrices(
            self.m_dem, 
            allow_undecomposed_hyperedges=self.allow_undecomposed_hyperedges
        )
        self.m_pcm = self.m_matrices.check_matrix
        bm = detector_error_model_to_check_matrices(
            self.m_dem, 
            allow_undecomposed_hyperedges = self.m_allow_undecomposed_hyperedges
        ) # Esto lo tendremos que cambiar


        # dem_data = DemData()
        # dem_data.priors = bm.priors
        # dem_data.obs_matrix = bm.observables_matrix.toarray('F').astype(np.uint8)


        # dem_data.phen_obs_matrix = bm.edge_observables_matrix.toarray('F').astype(np.uint8)
        #     # Make the observables phenomenological matrix
        # dem_data.phen_check_matrix = self.dem_data_phen_check_matrix
        # dem_data.transfer_matrix = self.transfer_matrix

        self.bpotf = OBPOTF(
            bm.check_matrix, 
            self.m_p, 
            self.m_noise_type,
            ps_ext_dem_data=self.m_dem_data, 
            po_ext_bp_iters= self.m_po_ext_bp_iters
        )

        shots = stim.read_shot_data_file(path=dets_b8_in_path, format="b8", num_detectors=num_dets)
        predictions = np.zeros((num_shots, num_obs), dtype=bool)

        for i in range(num_shots):
            predictions[i, :] = self.decode(shots[i, :])

        stim.write_shot_data_file(data=predictions,
                                  path=obs_predictions_b8_out_path,
                                  format="b8",
                                  num_observables=num_obs
                                  )
    
    def decode(self, syndrome: np.ndarray) -> np.ndarray:
        return self.bpotf.decode(syndrome)