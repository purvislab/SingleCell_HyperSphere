import numpy as np 

def pca_new(X, L=None):
    """
    Principal Component Analysis for input data set.

    Parameters:
    - X: a NxD matrix which contains one D dimension data in each row.
    - L: the dimension to which the dataset will be reduced (default is min(N,D)).

    Returns:
    - U: a DxL matrix which stores principle components in columns.
    - lambda: a vector containing min(N,D) eigenvalues.
    - xc: the centroid of input data set.
    """
    N, D = X.shape  # N data points with dimension D.

    if L is None:
        L = min(N, D)

    xc = np.mean(X, axis=0)  # Obtain the centroid of the original dataset.
    X = X - np.tile(xc, (N, 1))  # Zero-mean data set.
    

    if D <= N:
        S = np.matmul(X.T, X)  # DxD Covariance matrix.
    else:
        S = np.matmul(X, X.T)  # NxN Gram matrix.


    # Assuming S is a NumPy array
    eigvalues, eigvectors = np.linalg.eig(S)

    # Extract the eigenvalues
    eigvalues, eigvectors = np.linalg.eig(S)

    # Sort the eigenvalues in descending order and get the corresponding indices
    sorted_indices = np.argsort(eigvalues)[::-1]
    lambda_ = eigvalues[sorted_indices]

    # Extract the corresponding eigenvectors
    U = eigvectors[:, sorted_indices[:L]]


# Note: lambda is a reserved keyword in Python, so I've used lambda_values instead

    if D > N:
        lv = lambda_[:L]
        z = np.where(lv > 0)
        lv[z] = 1.0 / np.sqrt(lv[z])
        U = np.matmul(X.T, U)
        U = np.matmul(U, np.diag(lv))

    lambda_ = lambda_ / N

    xc = np.array(xc)
    xc = xc.reshape(-1, 1)
    xc = xc.T

    return U, lambda_, xc

# Example usage:
# U, lambda_, xc = pca_new(X, L)

def Spherelets(X, d):
    # Find the best d-dimensional sphere to fit data X

    # Input: X = data matrix
    # Output: c = center of the spherelet
    #         r = radius of the spherelet
    #         V = the subspace where the sphere lies in the affine space c+V

    n, m = X.shape  # n = sample size


    if n > d + 1:  # If there are enough samples, fit the data by a sphere or a hyperplane

        # Do d+1 dimensional PCA first
        V, lambd, mu = pca_new(X, d+1)

        Y = np.ones((n, 1)) @ mu + (X - np.ones((n, 1)) @ mu) @ V @ V.T
        Y = np.array(Y)

        l = np.zeros((n,1))

        for i in range(n):
            l[i] = np.linalg.norm(Y[i, :]) ** 2

        lbar = np.mean(l)
        H = np.zeros((m, m))
        f = np.zeros((m,1))

        for i in range(n):
            H += np.matmul((mu - Y[i, :]).T, (mu - Y[i, :]))
            f += (l[i] - lbar) * (mu - Y[i, :]).T 

        H_inv = np.linalg.pinv(H)

        # Calculate the center c of the sphere
        c = mu.T + np.matmul(np.matmul(V, V.T), (-0.5 * np.matmul(H_inv, f) - mu.T))
        
        Riemd = np.zeros((n, 1))  # Initialize an array to store distances

        for i in range(n):

            distance =  np.sqrt(np.matmul((c.T - Y[i, :]),(c.T - Y[i, :]).T))
            Riemd[i] = distance  # Store the distance in the array

        # Calculate the mean of the distances
        r = np.mean(Riemd)
    else:
        c = None
        r = None
        V = None
    
    return c, r, V


# Example usage:
# c, r, V = Spherelets(X, d)


def Proj(X, r, c, V):
    ## X is n x p
    ## r is 1x1 
    ## c is p x 1
    ## V is p x (d+1)

    # Proj(X) = r*(X-C)*V./rownorm((X-C)*V).
    # Should use row-norm in the denominator. 
    # So rownorm((X-C)*V) gives you a n by 1 column vector, 
    # and ./ means divide by row. This Proj(X) will be n by d+1.

    num = ((X - c.T) @ V)
    denom = np.linalg.norm(num, axis = 1) #n x 1
    denom = denom[:, np.newaxis]
    # Divide by rown in denom
    res = r * num / denom

    #return dist
    return res

def ProjHyp(X, r, c, V):
    ## X is n x p
    ## r is 1x1 
    ## c is p x 1
    ## V is p x (d+1)

    ## returns projection of X into new space, but now the projection is n x p dimensions
    unos = np.ones((len(X), 1))
    center_ones = unos @ c.T

    num = ((X - (unos @ c.T )) @ (V @ V.T)) 
    denom = np.linalg.norm(num, axis = 1)
    denom = denom[:, np.newaxis]
    projd = center_ones + (r * (num / denom))

    #return dist
    return np.array(projd)