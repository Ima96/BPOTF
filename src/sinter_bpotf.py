import pathlib
import numpy as np
import sinter
import stim
from BPOTF import OBPOTF, NoiseType, DemData

from ldpc.ckt_noise.dem_matrices import detector_error_model_to_check_matrices

class SinterBpOtfDecoder(sinter.Decoder):
	"""Sinter decoder adapter for the BPOTF (BP+BP+OTF) decoder.

	Wraps OBPOTF so it can be used with sinter.collect() for benchmarking
	under circuit-level noise.

	The DEM is the only mandatory input. Everything else is derived from it
	automatically, but can be overridden by the user when pre-computed or
	custom matrices are available.

	The three-stage BP+BP+OTF pipeline is:
	  Stage 1: BP on the full DEM check matrix (pcm).
	  Stage 2: BP on the sparsified (edge) check matrix.
	  Stage 3: OTF (Kruskal) + BP on the OTF matrix.

	Note: OBPOTF is built lazily on the first call to decode_via_files,
	because sinter workers need to pickle the decoder object and the C++
	OBPOTF/DemData objects are not picklable.
	"""

	def __init__(
			self,
			dem,
			p,
			pcm = None,
			obs_matrix = None,
			transfer_matrix = None,
			phen_check_matrix = None,
			phen_obs_matrix = None,
			otf_matrix = None,
			bp_iters = None,
			decimation = 1e-9
			):
		"""
		Parameters
		----------
		dem : stim.DetectorErrorModel
			The detector error model. Used to derive the check matrix,
			observable matrix, priors, and sparsified matrices when they
			are not explicitly provided.
		p : float
			Physical error probability used for the OTF channel.
		pcm : scipy.sparse matrix or None
			Parity check matrix (detectors x error mechanisms). If None,
			derived from the DEM.
		obs_matrix : np.ndarray or None
			Observable matrix. If None, derived from the DEM.
		transfer_matrix : np.ndarray or None
			Transfer matrix mapping DEM to the sparsified detector graph.
			If provided, enables the three-stage BP+BP+OTF decode path.
		phen_check_matrix : np.ndarray or None
			Sparsified (phenomenological) parity check matrix. If None,
			derived from the DEM via the edge_check_matrix.
		phen_obs_matrix : np.ndarray or None
			Sparsified (phenomenological) observable matrix. If None,
			derived from the DEM via the edge_observables_matrix.
		otf_matrix : np.ndarray or None
			Matrix on which the OTF (Kruskal) stage is applied. If None,
			OBPOTF uses the sparsified check matrix as default.
		bp_iters : list of 3 ints or None
			Max BP iterations as [stage1_dem, stage2_phen, stage3_otf].
			If None, OBPOTF uses its internal defaults (30, 100, 100).
		decimation : float
			Decimation parameter for the OTF stage.
		"""
		self.m_p = p
		self.m_decimation = decimation
		self.m_otf_matrix = otf_matrix

		# Process DEM to get default matrices
		bm = detector_error_model_to_check_matrices(dem, allow_undecomposed_hyperedges=True)

		# PCM: use provided or derive from DEM (store as scipy sparse)
		self.m_pcm = pcm if pcm is not None else bm.check_matrix

		# Priors always come from the DEM
		self.m_priors = bm.priors

		# Observable matrix: use provided or derive from DEM
		if obs_matrix is not None:
			self.m_obs_matrix = obs_matrix
		else:
			self.m_obs_matrix = bm.observables_matrix.toarray('F').astype(np.uint8)

		# Transfer matrix
		self.m_transfer_matrix = transfer_matrix

		# Sparsified check matrix: use provided or edge_check_matrix from DEM
		if phen_check_matrix is not None:
			self.m_phen_check_matrix = phen_check_matrix
		else:
			self.m_phen_check_matrix = bm.edge_check_matrix.toarray('F').astype(np.uint8)

		# Sparsified observable matrix: use provided or edge_observables_matrix from DEM
		if phen_obs_matrix is not None:
			self.m_phen_obs_matrix = phen_obs_matrix
		else:
			self.m_phen_obs_matrix = bm.edge_observables_matrix.toarray('F').astype(np.uint8)

		# BP iterations: convert to numpy int32 array or None
		if bp_iters is not None:
			self.m_bp_iters = np.array(bp_iters, dtype=np.int32)
		else:
			self.m_bp_iters = None

		# OBPOTF is built lazily in _build_decoder(), not here,
		# because C++ objects (DemData, OBPOTF) are not picklable
		# and sinter needs to pickle this object for worker processes.
		self.bpotf = None

	def _build_decoder(self):
		"""Build the OBPOTF decoder from the stored configuration.

		Called lazily on the first decode_via_files call, inside the
		sinter worker process (after pickling).
		"""
		dem_data = DemData()
		dem_data.priors = self.m_priors
		dem_data.obs_matrix = self.m_obs_matrix
		dem_data.transfer_matrix = self.m_transfer_matrix
		dem_data.phen_check_matrix = self.m_phen_check_matrix
		dem_data.phen_obs_matrix = self.m_phen_obs_matrix

		self.bpotf = OBPOTF(
			self.m_pcm,
			self.m_p,
			NoiseType.E_CLN,
			self.m_otf_matrix,
			self.m_bp_iters,
			dem_data,
			self.m_decimation
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

		if self.bpotf is None:
			self._build_decoder()

		shots = stim.read_shot_data_file(path=dets_b8_in_path, format="b8", num_detectors=num_dets)
		predictions = np.zeros((num_shots, num_obs), dtype=bool)

		for i in range(num_shots):
			predictions[i, :] = self.bpotf.decode(shots[i, :])

		stim.write_shot_data_file(data=predictions,
								  path=obs_predictions_b8_out_path,
								  format="b8",
								  num_observables=num_obs
								  )

	def decode(self, syndrome: np.ndarray) -> np.ndarray:
		if self.bpotf is None:
			self._build_decoder()
		return self.bpotf.decode(syndrome)
