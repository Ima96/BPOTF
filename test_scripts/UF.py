import numpy as np
from ldpc import bp_decoder
import copy
import time
import galois

from module import BPOTF

class UF:
    
    
    def __init__(self,
                 H: np.array,
                 p : float,
                 ):
        self.H = H
        self.p = p
        # Para hacer el primer BP
        self._bpd = bp_decoder(
            H,
            error_rate = p
        )
        # Para hacer el segundo BP
        self._bpd2 = bp_decoder(
            H,
            error_rate = p,
            #max_iter = 5
        )
        self.rows = H.shape[0]
        self.columns = H.shape[1]
        
        self.rank = np.linalg.matrix_rank(H).astype(np.uint64)
        
        self.Hog = copy.deepcopy(self.H)
        
        zeros_rows = np.zeros((2, self.columns))
        self.Hog =  np.vstack((self.Hog, zeros_rows))
        # self.Hog = self.Hog.astype(int)
        
        # We add virtual checks so all columns have at least two non-trivial elements.
        for i in range(self.columns):
            ones_in_col = np.where(H[:,i] == 1)[0]
            if len(ones_in_col) == 1:
                if ones_in_col[0] < self.rows//2:
                    self.Hog[-2,i] = 1
                else:
                    self.Hog[-1,i] = 1

        self.index_array = np.array([i for i in range(0, self.columns)], dtype=np.uint64)

        # Definimos el número máximo de índices por columna
        max_nontrivial_per_col = np.max(np.sum(self.Hog == 1, axis=0))
        # Definimos la matriz de índices. Cuando un valor es -1, es que ya no está incluido.
        self.index_matrix = np.full((max_nontrivial_per_col, self.columns), -1, dtype=np.int64)
        
        for column in range(self.columns):
            # Get the row indices where the value is 1 in the current column
            row_indices = np.where(self.Hog[:, column] == 1)[0]
            # Place these indices in the corresponding column of index_matrix
            self.index_matrix[:len(row_indices), column] = row_indices
        
        #self.index_matrix = self.index_matrix.astype(np.int8)
    
    def SetCluster(self):
        """Set the clusters which will be grown via Union Find.

        Returns:
            cluster: array of duples. First element of the duple is the pointer to the root of the tree.
                second element is the weight of the tree.
            
        """
        cluster =  np.zeros((self.rows+2, 2), dtype=int)
        cluster[:, 0] = np.arange(self.rows+2)
        return cluster
    
    
    def sort_matrix(self, llrs):
        sorted_indices = np.argsort(llrs)
        # H_sorted = self.Hog[:,sorted_indices]
        # return H_sorted, sorted_indices
        return  sorted_indices
    
    def Find(self, x :int, cluster_array = np.array):
        while x != cluster_array[x,0]:
            x = cluster_array[x,0]
        return x, cluster_array[x,1]
    
    # def Union(self, elements:list):
    #     pass
    
    def Kruskal_hypergraph(self, sorted_indices: np.ndarray):
        # Empezamos el cluster
        cluster_array = self.SetCluster()
        # La siguiente columna  indica qué columnas nos quedaremos como árbol.
        columns_chosen = np.zeros(self.columns, dtype = bool)
        # Iteramos sobre todas las columnas de la matriz self.Hog
        for column in range(self.columns):
            # Checker nos sirve para ever que checks coge cada evento.
            checker = np.zeros(self.rows+2, dtype = bool)
            # Depths: (root, depth, boolean about increasing tree size)
            depths = [0, -1, False]
            boolean_condition = True
            # Miramos para la columna "sorted_indices[column]]" de la matriz "self.Hog", que filas son no triviales.
            non_trivial_elements_in_column = self.index_matrix[:,sorted_indices[column]]
            # Las siguientes lineas son el método Union Find que hablamos el otro día:
            for non_trivial_element in non_trivial_elements_in_column:
                # Miramos la raiz y la profundidad del árbol al que corresponden.
                root, depth = self.Find(non_trivial_element, cluster_array)
                # Si otro valor no trivial te lleva a ese árbol omitimos la columna.
                if checker[root]:
                    boolean_condition = False
                    break
                # Si no lo apuntamos en el checker.
                checker[root] = True
                # Mantenemos el root de profundida más grande.
                if depth > depths[1]:
                    depths = [root, depth, False]
                # Si hay más de un árbol de profundidad máxima, la profundidad del árbol aumenta.
                if depth == depths[1]:
                    depths[2] = True
            # Si la hyperarista no ha incidido más de una vez a ningún árbol la añadimos.
            if boolean_condition:
                non_trivial_checker = np.where(checker == 1)[0]
                for element in non_trivial_checker:
                    cluster_array[element,0] = depths[0]
                columns_chosen[sorted_indices[column]] = True
                if  depths[2]:
                    cluster_array[depths[0]][1] += 1
        # La siguiente lista nos indica qué columnas han sido elegidas.
        indices_columns_chosen = np.where(columns_chosen==True)[0]
        # returned_H = H_sorted[:-2,indices_columns_chosen]
        # new_rows, new_cols = returned_H.shape
        # if new_rows >= new_cols:
        #     returned_H = np.hstack((returned_H,np.zeros((new_cols, new_rows-new_cols+1))))
        return indices_columns_chosen
    
    def Kruskal_hypergraph_v2(self, sorted_indices: np.ndarray):
        return BPOTF.kruskal_on_hypergraph_v2_uf(self.index_matrix, sorted_indices, self.rows+2, self.index_array, self.rank)
        #return BPOTF.kruskal_on_hypergraph_v2_classical_uf(self.index_matrix, sorted_indices, self.rows+2, self.index_array)
                
    def decode(self, syndrome : np.array):
        
        # Usamos el primer BP para encontrar el error.
        recovered_error = self._bpd.decode(syndrome)
        # Si converge devolvemos el error.
        if self._bpd.converge:
            return recovered_error, 0
        # Si no converge, nos quedamos con las llrs.
        llrs = self._bpd.log_prob_ratios
        # Ordenamos las llrs y nos quedamos con los indices por orden
        #sorted_indices = self.sort_matrix(llrs)
        #a = time.perf_counter()
        # Nos quedamos con las columnas linearmente independientes.
        #columns_chosen = self.Kruskal_hypergraph(sorted_indices)
        columns_chosen = self.Kruskal_hypergraph_v2(llrs)
        #b = time.perf_counter()
        #average_time = b-a
        # Para no tener que iniciar de nuevo una clase bp_decoder, usamos self._bp2.
        # Para que tenga en cuenta solo las columnas que nos interesan, lo iniciamos, le damos probabilidad 0 a las columnas que no usamos 
        # y probabilidad self.p a las que sí.
        updated_probs = np.zeros(self.columns)
        updated_probs[columns_chosen] = self.p
        self._bpd2.update_channel_probs(updated_probs)
        # Luego le damos a decode.
        second_recovered_error = self._bpd2.decode(syndrome)
        # if not self._bpd2.converge:
        #     print('Main error')

        average_time = 0
        return second_recovered_error, average_time
    
    
    
if __name__ == "__main__":
    pass