import numpy as np
from numba import njit, prange
from multiprocessing import Pool
from multiprocessing import Process, Manager,cpu_count

from sklearn.utils import compute_class_weight


# Compute CME score.
@njit
def compute_CME(X, i, j, sum_ary):
    min_ary = np.minimum(X[i,:], X[j,:])
    min_sum = min_ary.sum()

    ratio_x = min_sum / sum_ary[i]
    ratio_y = min_sum / sum_ary[j]

    cme = 1 - max(ratio_x, ratio_y)

    return cme


# Input X: gene * cell. A row is a gene
#@njit(parallel=False, fastmath=False, nopython=False)
@njit(parallel=True, fastmath=True, nopython=True)
def CME_sym_numba(X: np.ndarray[np.float32], feature_indices=None) -> np.ndarray[np.float32]:
    # Initialize the CME matrix.
    if feature_indices is None:
        feature_indices = np.arange(X.shape[0], dtype=np.int64)

    feature_size = len(feature_indices)

    cme_mtx = np.zeros((feature_size, feature_size), dtype=np.float32)

    # Prepare the index mapping
    (i_ind, j_ind) = np.triu_indices(feature_size)

    # Compute sum by gene    
    sum_ary = X.sum(axis=1)

    for k in prange(len(i_ind)):
        i = i_ind[k]
        j = j_ind[k]

        x_idx = feature_indices[i]
        y_idx = feature_indices[j]

        cme_mtx[i, j] = cme_mtx[j, i]  = compute_CME(X, x_idx, y_idx, sum_ary)
    
    return cme_mtx


# Input X: gene * cell. A row is a gene
@njit(parallel=True, fastmath=True, nopython=True)
def CME_asym_numba(X: np.ndarray[np.float32], data_indices=None, feature_indices=None) -> np.ndarray[np.float32]:

    # Initialize the CME matrix.
    if data_indices is None:
        data_indices = np.arange(X.shape[0], dtype=np.int32)

    if feature_indices is None:
        feature_indices = np.arange(X.shape[0], dtype=np.int32)

    data_size = len(data_indices)
    feature_size = len(data_indices)

    cme_mtx = np.zeros((data_size, feature_size), dtype=np.float32)
    
    # Compute sum by gene    
    sum_ary = X.sum(axis=1)

    # compute the minimum sum matrix
    #for i in prange(data_size):
    for i in range(data_size):
        for j in range(feature_size):
            data_idx = data_indices[i]
            feature_idx = feature_indices[j]

            cme_mtx[i, j] = compute_CME(X, data_idx, feature_idx, sum_ary)
            
    return cme_mtx

# input X has to be a numpy array of float32.
# Rows of X has to be genes. Column cells.
def CME(X: np.ndarray[np.float32], normalize=True, feature_indices=None) -> np.ndarray[np.float32]:

    #X = X.astype(X.shape, dtype=np.float32)

    if normalize:
        median_ary = np.apply_along_axis(lambda v: np.median(v[np.nonzero(v)]), 1, X)
        X_normed = X / median_ary[:, None]
    else:
        X_normed = X

    if feature_indices is None:
        feature_indices = np.arange(X.shape[0], dtype=np.int64)

    data_indices = np.arange(X.shape[0], dtype=np.int64)
    data_indices = data_indices[np.isin(data_indices, feature_indices, invert=True)]

    # Compute CME for symmetric part
    CME = CME_sym_numba(X_normed, feature_indices)

    # Compute CME for asymmetric part, if there is any.
    if len(data_indices) > 0:
        CME_asym = CME_asym_numba(X_normed, data_indices, feature_indices)

        CME = np.vstack((CME, CME_asym))

    return CME


def CME_correction(X: np.ndarray[np.float32], cme: np.ndarray[np.float32], iterations=10, normalize=True, feature_indices=None, tp_cutoff=0.95, fp_cutoff=0.05) -> np.ndarray[np.float32]:
    # Populate the null CME distribution by shuffling
    null_CME = []

    X_shuffled = X.copy()

    for row in X_shuffled:
        np.random.shuffle(row)

    X_shuffled = X_shuffled / X_shuffled.sum(axis=1)[None,:].T * 1e4
    
    cme_shuffled = CME(X_shuffled, normalize=normalize, feature_indices=feature_indices)

    null_CME.append(cme_shuffled)

    null_CME_ary = np.array(null_CME)

    cme_rank = np.zeros(cme.shape)
    for i in range(cme.shape[0]):
        for j in range(cme.shape[1]):
            cme_rank[i, j] = np.sum(null_CME_ary[:, i, j] < cme[i, j])

    tp_rank_thresh = iterations * tp_cutoff
    fp_rank_thresh = iterations * fp_cutoff

    cme_corrected = cme.copy()
    cme_corrected[cme_rank < fp_rank_thresh] = 0
    cme_corrected[(cme_rank > fp_rank_thresh) &
                  (cme_rank > tp_rank_thresh)] = 0.5
    
    return cme_corrected
