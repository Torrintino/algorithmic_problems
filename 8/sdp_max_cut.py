import cvxpy as cp
import numpy as np
import scipy


class SDP:

    def decomposition(self, A):
        """
        Computes the Cholesky decomposition of the Laplacian matrix of a
        regular graph.
        Input: adjacency matrix of the graph (numpy n x n matrix).
        Output: numpy n x n matrix.
        """
        n = A.shape[0]
        Y = cp.Variable(A.shape, symmetric=True)
        constraints = [Y >> 0]
        constraints += [Y[i,i] == 1 for i in range(n)]
        prob = cp.Problem(cp.Minimize(cp.trace(A @ Y)), constraints)
        prob.solve()
        Y = Y.value + np.eye(n) * 1e-9
        return scipy.linalg.cholesky(Y)

    def sdp(self, A, seed=2):
        """
        Runs the SDP for max-cut on the input graph.
        Input: adjacency matrix of the graph (numpy n x n matrix).
        Output: one of the two partitions of the cut.
        """
        n = A.shape[0]
        B = self.decomposition(A)
        r = np.random.uniform(size=n)
        r /= np.linalg.norm(r)
        partition = []
        for i, ui in enumerate(B.T):
            if (ui.T.dot(r)) >= 0:
                partition.append(i)
        return partition
