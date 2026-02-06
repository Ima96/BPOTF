import pathlib
import numpy as np
import sinter
import stim
from BPOTF import OBPOTF, NoiseType, DemData

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
			transfer_matrix = None,
			decimation = 1e-9
			):
		
		if noise_type != NoiseType.E_CLN:
			raise ValueError(
				f"Configuration Error: sinter only should be used for NoiseType.E_CLN..."
			)

		self.m_p = p
		self.m_noise_type = noise_type
		self.m_po_otf_csc_mat = po_otf_csc_mat
		self.m_po_ext_bp_iters = po_ext_bp_iters
		self.m_transfer_matrix = transfer_matrix
		self.m_decimation = decimation

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
		self.process_dem()

		self.bpotf = OBPOTF(self.m_pcm,
									self.m_p,
									self.m_noise_type,
									self.m_po_otf_csc_mat,
									self.m_po_ext_bp_iters,
									self.m_dem_data,
									self.m_decimation)

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

	# TODO: I think that we should offer to register a function by the user so that the method to generate the
	# phenomenological matrices can be altered, and also tested. Further investigation for this, as it can be limited by
	# the sinter parallelization etc.
	def process_dem(self):
		self.m_dem_data = DemData()
		bm = detector_error_model_to_check_matrices(self.m_dem, allow_undecomposed_hyperedges=True)
		
		self.m_pcm = bm.check_matrix
		self.m_dem_data.transfer_matrix = self.m_transfer_matrix
		self.m_dem_data.priors = bm.priors
		self.m_dem_data.obs_matrix = bm.observables_matrix.toarray('F').astype(np.uint8)

		# Previosly, these two were generated like this...
		# self.m_dem_data.phen_obs_matrix = bm.edge_observables_matrix.toarray('F').astype(np.uint8)
		# self.m_dem_data.phen_check_matrix = self.dem_data_phen_check_matrix

		# This is how I have observed that these two were generated in the testing scripts...
		nonzero_counts = bm.check_matrix.sum(axis=0).A1
		selected_cols = np.where(nonzero_counts <= 3)[0]
		temp = bm.observables_matrix[:, selected_cols]
		self.m_dem_data.phen_obs_matrix = temp.toarray('F').astype(np.uint8)
		#ssake the observables phenomenological matrix
		temp2 = bm.check_matrix[:, selected_cols]
		self.m_dem_data.phen_check_matrix = temp2.toarray('F').astype(np.uint8)



