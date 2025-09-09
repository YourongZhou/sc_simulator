import numpy as np
from numba import njit, prange
from multiprocessing import Pool
from multiprocessing import Process, Manager,cpu_count

@njit(parallel=True, fastmath=True, nopython=True)
def CME_numba(X):

    # _, gene_num = X.shape
    gene_num, _ = X.shape
    min_result = np.zeros((gene_num, gene_num), dtype=np.int64)
    (i_ind, j_ind) = np.triu_indices(gene_num)

    for k in prange(len(i_ind)):
        i = i_ind[k]
        j = j_ind[k]
        min_ary = np.minimum(X[i,:], X[j,:])
        min_result[j,i] = min_result[i,j] = sum(min_ary)      

    return min_result  

def CME(X):

    min_res = CME_numba(X)
    
    sum_x = np.sum(X, axis=1)
    ratio_x = min_res / sum_x[:, None]
    ratio_y = min_res / sum_x[None, :]
    result = 1 - np.maximum(ratio_x, ratio_y)

    return result.T


@njit(parallel=True, fastmath=True, nopython=True)
def CME_numba_subset(X, a_indices):

    gene_num, _ = X.shape
    a_size = len(a_indices)

    min_result = np.zeros((gene_num, a_size), dtype=np.int64)
    
    for i in prange(gene_num):
        for k in range(a_size):
            j = a_indices[k]
            min_ary = np.minimum(X[i,:], X[j,:])
            min_sum = sum(min_ary)
            min_result[i, k] = min_sum
    
    return min_result  

def CME_subset(X, a):
    """
    计算所有基因与指定基因列表a之间的CME矩阵
    X: 输入矩阵
    a: 基因索引列表
    """

    a_indices = np.asarray(a, dtype=np.int64)

    min_res = CME_numba_subset(X, a_indices)
    sum_x = np.sum(X, axis=1)

    ratio_x = min_res / sum_x[:, None]
    ratio_y = min_res / sum_x[a_indices][None, :]
    result = 1 - np.maximum(ratio_x, ratio_y)
    
    return result