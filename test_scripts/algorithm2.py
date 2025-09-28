from qecsim.models.rotatedplanar import RotatedPlanarCode
from packaging.version import Version
from ldpc import __version__ as ldpc_version
ldpc_v2 = Version(ldpc_version) >= Version("2.0.0")
print("Using LDPC version v{}".format(ldpc_version))
if ldpc_v2 is True:
    from ldpc import BpDecoder
    from ldpc import BpOsdDecoder
else:
    from ldpc import bp_decoder
    from ldpc import bposd_decoder
import numpy as np
import random
from qecsim import paulitools as pt
# from networkx import find_cycle
import networkx as nt
from copy import deepcopy
#from UF import UF
import time

import csv

#from module import BPOTF #bpbp
import BPOTF

from scipy import sparse

def error_generation(p, n):
    """Depolarizing error generation

    Args:
        p (float): Probability of error
        n (int): number of qubits
    """
    error = np.zeros(2*n)
    probabilities = [(1 - p), p/3, p/3, p/3]
    results = ['I', 'X', 'Y', 'Z']

    # Generate n random realizations based on the probabilities
    realizations = random.choices(results, probabilities, k=n)
    for index, realization in enumerate(realizations):
        if realization == 'X':
            error[index] = 1
        elif realization == 'Z':
            error[index+n] = 1
        elif realization == 'Y':
            error[index] = 1
            error[index+n] = 1
    return error.astype(int)

def order_matrix_by_vector(vector, matrix):
    # Sort the indices of the vector based on its values
    sorted_indices = np.argsort(vector)
    # Reorder the columns of the matrix based on the sorted indices
    ordered_matrix = matrix[:, sorted_indices]
    return ordered_matrix, sorted_indices


def kruskal_on_hypergraph(Hog):
    """
    We now need to produce the Tanner graph via nt.graph

    In contination follow the following route:
    1 add first hyperedge to the graph as different hyperedges
    2 initiate following loop:
    initial_graph = nt.graph()
    for column in H:
        graph_to_try = initial_graph.copy()
        for i in np.where(H[:,column]==0):
            graph_to_try.add_edge((i+H.shape[1],column))
        if graph_to_try.has_loops():
            continue
        else:
            initial_graph = graph_to_try.copy()
            # Checkear que no hayan ya más de n-k columnas
            
     """
     
    initial_graph = nt.Graph()
    rows, columns = Hog.shape
    column_to_square = np.zeros(rows)
    
    # Just before initializing the process we introduce two additional rows on H, introducing virtual checks.
    zeros_rows = np.zeros((2, columns))
    H =  np.vstack((Hog, zeros_rows))
    
    
    for i in range(columns):
        if len(np.where(H[:,i] == 1)[0]) == 1:
            if np.where(H[:,i] == 1)[0][0] < rows:
                H[-2,i] = 1
            else:
                H[-1,i] = 1
                
    
    column_number = 1
    # Rows (Checks) are numbers 0 to rows-1
    # Columns (hypergraphs) are numbers rows-1 to rows+columns-1
    # column_to_square is a vector with the columns that will be considered for the square matrix:
    # final_matrix = H[:, column_to_square]
    
    # We begin by adding the first column:
    for edge in np.where(H[:,0] == 1)[0]:
        initial_graph.add_edge(edge, rows+2)
    
    
    for i in range(1,columns):
        # column_to_consider = H[:,i]
        
        Graph_to_check = deepcopy(initial_graph)
        # edges_to_add = []
        for edge in np.where(H[:,i] == 1)[0]:
            # edges_to_add.append((edge, i+rows+2))
            Graph_to_check.add_edge(edge, i+rows+2)
        # initial_graph.add_edges_from(edges_to_add)
        # if len(list(nt.simple_cycles(initial_graph))) > 0:
        # if len(list(nt.cycle_basis(Graph_to_check))) > 0:
        #     # initial_graph.remove_edges_from(edges_to_add)
        #     continue
        # initial_graph = Graph_to_check
        # column_to_square[column_number] = i
        # column_number += 1
        # if column_number == rows:
        #     break
        try: 
            nt.find_cycle(Graph_to_check)
            continue
        except Exception:
            initial_graph = Graph_to_check
            column_to_square[column_number] = i
            column_number += 1
            if column_number == rows:
                break
    

    assert column_to_square[-1] != 0, " Ha habido un error, no ha encontrado n-k columnas independientes"
    
    H_square = H[:,column_to_square.astype(int)]
    
    empty_column = np.zeros((H_square.shape[0], 1))  # Create an empty column filled with zeros
    matrix_with_empty_column = np.hstack((H_square, empty_column))

    return matrix_with_empty_column[:-2,:], column_to_square
    
        
    
if __name__ == "__main__":
    # distances = [3]
    #distances = [11]
    distances = [3, 5, 7, 9, 11]
    NMCs = [10**4, 10**4, 10**4, 10**4, 10**4, 10**4, 10**4, 10**4, 10**4, 10**3, 10**3, 10**3, 10**3]
    ps = np.linspace(0.01, 0.13, num=13)
    PlsBP = {}
    PlsBPOSD = {}
    PlsBPBP = {}

    times_BPOSD = {}
    times_BPBP = {}

    # f = open("times.csv", "w")

    for distance in distances:
        myCode = RotatedPlanarCode(distance, distance)
        pcm = np.zeros((myCode.stabilizers.shape[0], myCode.stabilizers.shape[1]), dtype = int)

        pcm[:,:pcm.shape[1]//2] = myCode.stabilizers[:,pcm.shape[1]//2:]
        pcm[:,pcm.shape[1]//2:] = myCode.stabilizers[:,:pcm.shape[1]//2]
        
        PlsBP[distance] = []
        PlsBPBP[distance] = []
        PlsBPOSD[distance] = []

        times_BPOSD[distance] = []
        times_BPBP[distance] = []

        print(f'Distance: {distance}')
        print('-------------------------------------------------')
        for index, p in enumerate(ps):
        
            if ldpc_v2 is True:
                _bp = BpDecoder(
                        pcm,
                        error_rate = float(p)
                )

                _bposd = BpOsdDecoder(
                    pcm,
                    error_rate= float(p)
                )
            else:
                _bp = bp_decoder(
                    pcm,
                    max_iter=30,
                    error_rate = p
                )
                
                _bposd = bposd_decoder(
                    pcm,
                    max_iter=30,
                    error_rate = p,
                    osd_method = "osd_0"
                )

            
            # _uf = UF(
            #     pcm,
            #     p
            # )

            # pcm = pcm.astype(np.uint8)
            pcm = sparse.csc_matrix(pcm, dtype=np.uint8)
            _bpotf = BPOTF.OBPOTF(
                pcm,
                p
            )
            
            PlBP = 0
            PlBPBP = 0
            PlBPOSD = 0
            time_av_BPOSD = 0
            time_max_BPOSD = 0
            time_av_BPBP = 0
            time_max_BPOTF = 0
            kruskal_time_list = []
            
            #print(pcm)
            #print(p)
            for iteration in range(NMCs[index]):
                error = error_generation(p, myCode.n_k_d[0])
                # error = np.zeros(myCode.n_k_d[0]*2, dtype = int)
                # error[3] = 1
                syndrome = pt.bsp(error, myCode.stabilizers.T)
                #print(syndrome)
                
                # syndrome = pt.bsp(error, pcm.T)
                # BPOSD decoder
                #----------------------------
                a = time.time()
                recovered_error_BPOSD = _bposd.decode(syndrome)
                b = time.time()
                time_av_BPOSD += (b-a)/NMCs[index]
                times_BPOSD[distance].append(b-a)
                time_max_BPOSD = max(time_max_BPOSD, (b-a))

                recovered_error = _bp.decode(syndrome)

                a = time.time()
                #recovered_error_BPBP, kruskal_time = _uf.decode(syndrome)
                # kruskal_time = 0.0
                # recovered_error_BPBP = _bpotf.decode(syndrome.astype(np.int32))
                recovered_error_BPBP = _bpotf.decode(syndrome.astype(np.uint8))
                b = time.time()
                time_av_BPBP += (b-a)/NMCs[index]
                times_BPBP[distance].append(b-a)
                time_max_BPOTF = max(time_max_BPOTF, (b-a))
                
                if np.any(pt.bsp(recovered_error_BPOSD ^ error, myCode.logicals.T) == 1):
                    PlBPOSD += 1/NMCs[index]
                    
                if np.any(pt.bsp(recovered_error_BPBP ^ error, myCode.logicals.T) == 1):
                    PlBPBP += 1/NMCs[index]
                #----------------------------
                
                if _bp.converge:
                    if np.any(pt.bsp(recovered_error ^ error, myCode.logicals.T) == 1):
                        PlBP += 1/NMCs[index]
                    continue
                else:
                    PlBP += 1/NMCs[index]

                # if kruskal_time > 0:
                #     kruskal_time_list.append(kruskal_time)
                
            
            PlsBP[distance].append(PlBP)
            PlsBPOSD[distance].append(PlBPOSD)
            PlsBPBP[distance].append(PlBPBP)
                
            print(f'Physical error: {p}')
            print(f'Error BP: {PlBP}')
            print(f'Error BPOSD: {PlBPOSD} with average time {time_av_BPOSD} and max time {time_max_BPOSD}')
            print(f'Error BPBP: {PlBPBP} with average time {time_av_BPBP} and max time {time_max_BPOTF}')
            
            # f.write(f'Distance: {distance} and physical error: {p}\n')
            # f.write('-------------------------------------------------\n')
            # f.write('Max BPOSD time: {}\n'.format(max(times_BPOSD[distance])))
            # f.write('Max BPBP time: {}\n\n'.format(max(times_BPBP[distance])))
            
            # if len(kruskal_time_list) > 0:
            #     print(f'Average Kruskal time {sum(kruskal_time_list)/len(kruskal_time_list)}')
            # print('\n')
        
    # f.close()