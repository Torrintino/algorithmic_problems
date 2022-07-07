import networkit as nk
import numpy as np
import scipy.sparse as sparse
from scipy.sparse.linalg import eigs
from scipy.sparse import csr_matrix


def csr_to_networkit(A):
        """ Converts an adjacency matrix to a NetworKit graph. """

        g=nk.graph.Graph(A.shape[0])
        row_idx, col_idx = np.nonzero(sparse.triu(A))
        for u, v in zip(row_idx, col_idx):
            g.addEdge(u, v)
        
        return g

class SpectralClustering:
    def compute_modularity_matrix(self, A):
        """ Computes the modularity matrix of the input adjacency matrix A. """

        A = csr_to_networkit(A)
        row = []
        col = []
        data = []

        for i in range(A.numberOfNodes()):
            for j in range(A.numberOfNodes()):
                aij = A.weight(i, j)
                row.append(i)
                col.append(j)
                bij = ( aij - (A.degree(i) * A.degree(j)) / (2 * A.numberOfEdges()) )
                data.append(bij)


        B = csr_matrix((data, (row, col)))
        return B
                 

        

    def compute_modularity(self, A, B, C):
        """ Computes the modularity of the partitioning of the vertices in the communities
            C and V \ C.
        """
        # TODO write your implementation here

        G = csr_to_networkit(A)
        m = G.numberOfEdges()
        x = np.array(C)
        xt = x.transpose()
 
        return xt.dot(B.dot(x))/(4*m)

    def run_spectral_clustering(self, A):
        """ Runs the spectral partitioning algorithm on the adjacency matrix A. """
        # TODO write your implementation here
       
        B = self.compute_modularity_matrix(A)

        v, w = eigs(B, k=B.shape[0]-2)

        best_ev = -999999999999999
        best_index = -1
        for i, ev in enumerate(v):

            if ev > best_ev:
                best_ev = ev;
                best_index = i
                z1 = w[:,best_index]

        C = [0 for _ in range(len(z1))]

        for i, z in enumerate(z1):
                    if z <= 0:
                        C[i] = 1
                    else:
                        C[i] = -1


        return C











