import scipy.io as sio
from scipy.sparse import csc_matrix
from SlidingWindowDecoder.src.build_circuit import build_circuit
from SlidingWindowDecoder.src.codes_q import create_bivariate_bicycle_codes#, create_circulant_matrix
from beliefmatching import detector_error_model_to_check_matrices
import numpy as np

import sinter
from BPOTF import __version__ as bpotf_version
print(f"BPOTF version is v{bpotf_version}")
from BPOTF import NoiseType, DemData
from sinter_bpotf import SinterBpOtfDecoder

#TODO -> Remove unnecessary imports.
#TODO -> This test should be done with surface codes.

def main():
	BB_TYPE =  108

	ps = [1e-3, 1.5e-3, 2e-3, 3e-3][::-1]
	p = ps[0]

	if BB_TYPE == 72:
		# [72, 12, 6] último número es el numero de rondas
		code, A_list, B_list = create_bivariate_bicycle_codes(6, 6, [3], [1,2], [1,2], [3])
		d = 6
		# transfer_mat = sio.loadmat('transfermatrices/transferMatrixcodel6m6.mat')['transfMat']
		transfer_mat = sio.loadmat('transfermatrices/PhenoTransf72Test.mat')['transfMatDEMtoPheno']
		phenoDEM = sio.loadmat('transfermatrices/PhenoTransf72Test.mat')['dem_pheno']
		phenoOBS = sio.loadmat('transfermatrices/PhenoTransf72Test.mat')['obsphen']
	elif BB_TYPE == 90:
		# [[90,8,10]]
		code, A_list, B_list = create_bivariate_bicycle_codes(15, 3, [9], [1,2], [2,7], [0])
		d = 10
		transfer_mat = sio.loadmat('transfermatrices/PhenoTransf90Test.mat')['transfMatDEMtoPheno']
		phenoDEM = sio.loadmat('transfermatrices/PhenoTransf90Test.mat')['dem_pheno']
		phenoOBS = sio.loadmat('transfermatrices/PhenoTransf90Test.mat')['obsphen']
	elif BB_TYPE == 108:
		# [108, 8, 10]
		code, A_list, B_list = create_bivariate_bicycle_codes(9, 6, [3], [1,2], [1,2], [3])
		d = 10
		# transfer_mat = sio.loadmat('transfermatrices/transferMatrixcodel9m6.mat')['transfMat']
		transfer_mat = sio.loadmat('transfermatrices/PhenoTransf108Test.mat')['transfMatDEMtoPheno']
		phenoDEM = sio.loadmat('transfermatrices/PhenoTransf108Test.mat')['dem_pheno']
		phenoOBS = sio.loadmat('transfermatrices/PhenoTransf108Test.mat')['obsphen']
	elif BB_TYPE == 144:
		# [144, 12, 12]
		code, A_list, B_list = create_bivariate_bicycle_codes(12, 6, [3], [1,2], [1,2], [3])
		d = 12
		# transfer_mat = sio.loadmat('transfermatrices/transferMatrixcodel12m6.mat')['transfMat']
		transfer_mat = sio.loadmat('transfermatrices/PhenoTransf144Test.mat')['transfMatDEMtoPheno']
		phenoDEM = sio.loadmat('transfermatrices/PhenoTransf144Test.mat')['dem_pheno']
		phenoOBS = sio.loadmat('transfermatrices/PhenoTransf144Test.mat')['obsphen']
	elif BB_TYPE == 288:
		#  [288, 12, 18]
		code, A_list, B_list = create_bivariate_bicycle_codes(12, 12, [3], [2,7], [1,2], [3])
		d = 18
		# transfer_mat = sio.loadmat('transfermatrices/BB288TransfDemsObs.mat')['transfMat']
		# transfer_mat = csc_matrix(sio.loadmat('transfermatrices/BB288CSC.mat')['transfMat']).toarray('F')
		transfer_mat = sio.loadmat('transfermatrices/PhenoDEMSparseNEW.mat')['transfMatDEMtoPheno'].toarray('F')
		phenoDEM = sio.loadmat('transfermatrices/PhenoDEMSparseNEW.mat')['dem_pheno']
		phenoOBS = sio.loadmat('transfermatrices/PhenoDEMSparseNEW.mat')['obsphen']
	else:
		raise Exception("No such option!")



	ps = 3e-3
	circuit = build_circuit(code, A_list, B_list, 
									p=ps, # physical error rate
									num_repeat=d, # usually set to code distance
									z_basis=True,   # whether in the z-basis or x-basis
									use_both=False, # whether use measurement results in both basis to decode one basis
	)
	dem = circuit.detector_error_model()
	bm = detector_error_model_to_check_matrices(dem, allow_undecomposed_hyperedges = True)
	
	dem_data = DemData()
	dem_data.priors = bm.priors
	dem_data.obs_matrix = bm.observables_matrix.toarray('F').astype(np.uint8)
	dem_data.transfer_matrix = transfer_mat.astype(np.uint8)
	dem_data.phen_obs_matrix = phenoOBS
	dem_data.phen_check_matrix = phenoDEM

	max_iters_dem = 100
	max_iters_sparse = 400
	max_iters_forest = 100

	bp_iterations = [max_iters_dem, max_iters_sparse, max_iters_forest]

########################################################################################################################

	# LAS SIGUIENTES LINEAS LAS USO PARA DEBUGGEAR, BORRAR CUANDO SINTER FUNCIONE
 
	sinter_BPOTF = SinterBpOtfDecoder(
							ps, 
							dem_data,
							po_ext_bp_iters = bp_iterations,
							pcm = bm.check_matrix
							)
########################################################################################################################


	# TODO ADAPT sinter to next
	task = sinter.Task(
		circuit=circuit,
		collection_options=sinter.CollectionOptions(max_shots=1_000_000, max_errors=1_000),
	)

	# iterations: 35, default and 800
	samples = sinter.collect(
					num_workers=1,
					tasks=[task],
					decoders=["bpotf"],
					custom_decoders={
						"bpotf": SinterBpOtfDecoder(
							ps, 
							dem_data,
							po_ext_bp_iters = bp_iterations,
							)
						},
					print_progress=True,
					save_resume_filepath="sinter_results/results.csv"
				)

	# Print samples as CSV data.
	print(sinter.CSV_HEADER)
	for sample in samples:
		print(sample.to_csv_line())

# This block seems mandatory on Windows when using siner/multiprocessing
if __name__ == "__main__":
	main()




