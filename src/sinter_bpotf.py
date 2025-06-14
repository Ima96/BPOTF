import pathlib
import numpy as np
import sinter
import stim
from BPOTF import OBPOTF, NoiseType

from ldpc.ckt_noise.dem_matrices import detector_error_model_to_check_matrices

class SinterBpOtfDecoder(sinter.Decoder):
    """
    TODO: Insert here a description of what is this for.
    """

    def __init__(
            self,
            p,
            noise_type = NoiseType.E_CLN,
            po_otf_csc_mat = None,
            po_ext_bp_iters = None,
            ps_ext_dem_data = None
            ):
        
        self.m_p = p
        self.m_noise_type = noise_type
        self.m_po_otf_csc_mat = po_otf_csc_mat
        self.m_po_ext_bp_iters = po_ext_bp_iters
        self.m_ps_ext_dem_data = ps_ext_dem_data

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
        self.m_matrices = detector_error_model_to_check_matrices(self.m_dem, allow_undecomposed_hyperedges=True)
        self.m_pcm = self.m_matrices.check_matrix

        self.bpotf = OBPOTF(
            self.m_pcm,
            self.m_p,
            self.m_noise_type,
            self.m_po_otf_csc_mat,
            self.m_po_ext_bp_iters,
            self.m_ps_ext_dem_data
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