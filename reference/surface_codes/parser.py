import numpy as np

def parser(circuit, parity_check_matrix):
    """ This function should input the detector error model of a surface code and return a matrix desired for OTF computation. The arguments which are inputted are a circuit which produces 
    the syndrome extraction for the surface code and the parity check matrix of the code.
    circuit: is the stim circuit for the surface code.
    parity_check_matrix: is the parity check matrix that will be considered for the OTF computation."""
    lines = str(circuit).splitlines()
    zchecks = []
    for line in lines:
        if line.startswith("H "):
            number_strings = line[2:].split()  # line[2:] removes the "H " part
            zchecks = [int(num) for num in number_strings] 
            break
    assert len(zchecks) > 0, " Z checks not detected."
    indices = []
    
    for line in lines:
        if line.startswith("MR "):
            number_strings = line[2:].split()  # line[2:] removes the "MR " part
            numbers_2 =  [int(num) for num in number_strings] 
            if zchecks[0] not in numbers_2:
                continue
            for zcheck in zchecks:
                indices.append(numbers_2.index(zcheck))
            break


    zero_rows = np.zeros((2, parity_check_matrix.shape[1]))

    otf_matrix = np.vstack((parity_check_matrix, zero_rows))

    for column in range(parity_check_matrix.shape[1]):
        pair = np.where(parity_check_matrix[:, column] == 1)[0]
        if len(pair) == 1:
            #First extraction round, there are only x checks
            if pair[0] < int(((parity_check_matrix.shape[0]-2) - 1) // 2):
                otf_matrix[parity_check_matrix.shape[0], column] = 1
            # On the last extraction round, there are only x checks
            elif pair[0] > parity_check_matrix.shape[0] - int(((parity_check_matrix.shape[0]-2) - 1) // 2):
                otf_matrix[parity_check_matrix.shape[0], column] = 1
            # Bulk extraction rounds, both xchecks and zchecks.
            else:
                # round = (pair[0]- int(((d**2)-1)//2)) // ((d**2)-1)
                check = (pair[0] - int(((parity_check_matrix.shape[0]-2) - 1) // 2)) % (parity_check_matrix.shape[0]-2)
                if check in indices:
                    otf_matrix[parity_check_matrix.shape[0] + 1, column] = 1
                else:
                    otf_matrix[parity_check_matrix.shape[0], column] = 1

    return otf_matrix


def sc_otf_matrix_computation(H):
    """
    Inputs a sparsified detector error model for a surface code and returns an OTFF matrix where all columns have two ones.

    Args:
        H (_type_): _description_
    """
    pass



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