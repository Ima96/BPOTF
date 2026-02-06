import scipy.io as sio
from scipy.sparse import csc_matrix
from SlidingWindowDecoder.src.build_circuit import build_circuit
from SlidingWindowDecoder.src.codes_q import create_bivariate_bicycle_codes#, create_circulant_matrix

# from packaging.version import Version
# from ldpc import __version__ as ldpc_version
# ldpc_v2 = Version(ldpc_version) >= Version("2.0.0")
# print("Using LDPC version v{}".format(ldpc_version))
# if ldpc_v2 is True:
# 	from ldpc import BpDecoder as bp_decoder
# 	from ldpc import BpOsdDecoder as bposd_decoder
# else:
# 	from ldpc import bp_decoder
# 	from ldpc import bposd_decoder

import sinter
from BPOTF import __version__ as bpotf_version
print(f"BPOTF version is v{bpotf_version}")
from BPOTF import SinterBpOtfDecoder, NoiseType

#TODO -> Remove unnecessary imports.
#TODO -> This test should be done with surface codes.

def main():
	BB_TYPE =  72

	ps = [1e-3, 1.5e-3, 2e-3, 3e-3][::-1]
	p = ps[0]

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
		transfer_mat = csc_matrix(sio.loadmat('transfermatrices/BB288CSC.mat')['transfMat'])
	else:
		raise Exception("No such option!")

	ps = 3e-3
	circuit = build_circuit(code, A_list, B_list, 
									p=ps, # physical error rate
									num_repeat=d, # usually set to code distance
									z_basis=True,   # whether in the z-basis or x-basis
									use_both=False, # whether use measurement results in both basis to decode one basis
									)

	max_iters_dem = 100
	max_iters_sparse = 400
	max_iters_forest = 100

	bp_iterations = [max_iters_dem, max_iters_sparse, max_iters_forest]

	# TODO ADAPT sinter to next
	task = sinter.Task(
		circuit=circuit,
		collection_options=sinter.CollectionOptions(max_shots=10, max_errors=10),
	)

	# iterations: 35, default and 800
	samples = sinter.collect(
					num_workers=1,
					tasks=[task],
					decoders=["bpotf"],
					custom_decoders={
						"bpotf": SinterBpOtfDecoder(
							p, 
							NoiseType.E_CLN,
							po_ext_bp_iters = bp_iterations,
							transfer_matrix = transfer_mat
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




