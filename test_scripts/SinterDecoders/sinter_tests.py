import stim, sinter
import numpy as np
# from BPOTF import SinterBpOtfDecoder # TODO ->  SinterBpOtfDecoder not successfully loaded from BPOTF
from sinter_bpotf import SinterBpOtfDecoder


def main():
    d = 3
    p = 0.001
    circuit = stim.Circuit.generated(
        "color_code:memory_xyz",
        rounds=d,
        distance=d,
        after_clifford_depolarization=p,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p,
        # code_task="color_code:rotated_memory_z",
    )

    task = sinter.Task(
        circuit=circuit,
        collection_options=sinter.CollectionOptions(max_shots=10, max_errors=10),
    )

    # iterations: 35, default and 800
    bp_iterations = np.array([35, -1, 800], dtype=np.int32)
    sinter.collect(
        num_workers=1,
        tasks=[task],
        decoders=["bpotf"],
        custom_decoders={
            "bpotf": SinterBpOtfDecoder(
                p=0.002,
                po_ext_bp_iters=bp_iterations,
                
            )
            },
        print_progress=True,
        save_resume_filepath="results.csv"
    )

# This block seems mandatory on Windows when using siner/multiprocessing
if __name__ == "__main__":
    main()