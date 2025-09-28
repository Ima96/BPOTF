import numpy as np

def otf_matrix_computer(circuit, parity_check_matrix, extraction_rounds):
    """ This function should input the detector error model of a surface code and return a matrix desired for OTF computation. The arguments which are inputted are a circuit which produces 
    the syndrome extraction for the surface code and the parity check matrix of the code.
    circuit: is the stim circuit for the surface code.
    parity_check_matrix: is the parity check matrix that will be considered for the OTF computation."""
    
    # Instead of indices, we are going to use a boolean array to indicate the checks for each row.
    
    significant_checks = np.zeros(parity_check_matrix.shape[0], dtype=bool)
    
    lines = str(circuit).splitlines()
    zchecks = []
    
    for line in lines:
        if line.startswith("H "):
            number_strings = line[2:].split()  # line[2:] removes the "H " part
            zchecks = [int(num) for num in number_strings] 
            break
    assert len(zchecks) > 0, " Z checks not detected."
    indices = []
    
    true_indices = False # State can only have three states 0: first round, 1: bulk d-1 rounds, 2: last round.
    
    for line in lines:
        if line.startswith("MR "):
            number_strings = line[2:].split()  # line[2:] removes the "MR " part
            numbers_2 =  [int(num) for num in number_strings]
            for index,number in enumerate(numbers_2):
                if number in zchecks:
                    indices.append(index)
            break
    
    for line in lines:
        if line.startswith("DETECTOR"):
            rec_index = line.split("rec[")[1].split("]")[0]
            negative_value = int(rec_index)
            if numbers_2[negative_value] in zchecks:
                true_indices = True
            break
    
    
    if true_indices:
        significant_checks[:len(zchecks)] = True
        significant_checks[-len(zchecks):] = True
        for round in range(extraction_rounds - 1):
            for index,number in enumerate(numbers_2):
                if number not in zchecks:
                    significant_checks[len(zchecks) + round * (len(numbers_2)) + index] = True
    else:
        for round in range(extraction_rounds - 1):
            for index,number in enumerate(numbers_2):
                if number in zchecks:
                    significant_checks[len(zchecks) + round * (len(numbers_2)) + index] = True
    

    zero_rows = np.zeros((2, parity_check_matrix.shape[1]))

    otf_matrix = np.vstack((parity_check_matrix, zero_rows))

    for column in range(parity_check_matrix.shape[1]):
        pair = np.where(parity_check_matrix[:, column] == 1)[0]
        if len(pair) == 1:
            #First extraction round, there are only x checks
            if significant_checks[pair[0]]:
                otf_matrix[-1, column] = 1
            else:
                otf_matrix[-2, column] = 1

    return otf_matrix




# zero_rows = np.zeros((2, self.H_phen.shape[1]))

# self.otf_matrix = np.vstack((otf_matrix, zero_rows))

# for column in range(self.H_phen.shape[1]):
#     pair = np.where(self.H_phen[:,column]==1)[0]
#     if len(pair)==1:
#         # First extraction round, only x checks
#         if pair[0] < int(((d**2)-1)//2):
#             self.otf_matrix[self.H_phen.shape[0], column] =  1
#         # Last extraction round, only x checks
#         elif pair[0] > self.H_phen.shape[0]-int(((d**2)-1)//2):
#             self.otf_matrix[self.H_phen.shape[0], column] =  1
#         # Bulk extraction rounds, both xchecks and zchecks.
#         else:
#             # round = (pair[0]- int(((d**2)-1)//2)) // ((d**2)-1)
#             check = (pair[0]- int(((d**2)-1)//2)) % ((d**2)-1)
#             if check in indices:
#                 self.otf_matrix[self.H_phen.shape[0]+1, column] =  1
#             else:
#                 self.otf_matrix[self.H_phen.shape[0], column] =  1