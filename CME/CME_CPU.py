import numpy as np
from numba import njit, prange


def CME_pva_cpu(cme: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                    null_cme_ary: np.ndarray[tuple[int, int, int], np.dtype[np.float32]]
                    ) -> np.ndarray[int, np.dtype[np.float32]]:
    # Populate the null CME distribution by shuffling
    count_matrix = np.sum(null_cme_ary < cme, axis=0)
    pval_mtx = count_matrix / float(null_cme_ary.shape[0])

    return pval_mtx


def shuffle_and_normalize_cpu(X: np.ndarray[tuple[int, int], np.dtype[np.float32]]
                               ) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
    
    # Shuffle and normalize
    rand_mat = np.random.rand(*X.shape)
    shuffled_idx = rand_mat.argsort(axis=1)
    shuffled_X = np.take_along_axis(X, shuffled_idx, axis=1)

    col_sum_ary = shuffled_X.sum(axis=0)
    shuffled_X = shuffled_X / col_sum_ary[None, :] * 1e4

    return shuffled_X


# input X has to be a numpy array of float32.
# Rows of X has to be genes. Column cells.
def CME_cpu(X: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                normalize: bool = False,
                feature_indices: np.ndarray[int, np.dtype[np.int32]] = None
                ) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
    # Normalize the input matrix if needed
    if normalize:
        X_nan = np.where(X == 0.0, np.nan, X)
        medians = np.nanmedian(X_nan, axis=1)
        X_normed = X / medians[:, np.newaxis]
    else:
        X_normed = X

    # Compute sum by gene 
    sum_ary = X_normed.sum(axis=1)

    # Populate feature and data indices
    if feature_indices is None:
        feature_indices = np.arange(X.shape[0], dtype=np.int32)

    data_indices = np.arange(X.shape[0], dtype=np.int32)
    data_indices = data_indices[np.isin(data_indices, feature_indices, invert=True)]

    # Compute the symmetric part [upper triangle of a squire matrix]
    cme = CME_sym_numba(X_normed, sum_ary, feature_indices)

    # # Compute the symmetric part [rectangular matrix]
    if len(data_indices) > 0:
        cme_asym = CME_asym_numba(X_normed, sum_ary, data_indices, feature_indices)
        cme = np.vstack((cme, cme_asym))

    return cme


# Compute CME score.
@njit
def compute_CME(X, i, j, sum_ary):
    #Get the sum of the min between genes.
    min_ary = np.minimum(X[i,:], X[j,:])
    min_sum = min_ary.sum()

    ratio_x = min_sum / sum_ary[i]
    ratio_y = min_sum / sum_ary[j]

    cme = 1 - max(ratio_x, ratio_y)

    return cme


# Input X: gene * cell. A row is a gene
#@njit(parallel=False, fastmath=False, nopython=False)
@njit(parallel=True, fastmath=True, nopython=True)
def CME_sym_numba(X: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                    sum_ary: np.ndarray[int, np.dtype[np.float32]],
                    feature_indices: np.ndarray[int, np.dtype[np.int32]]
                    ) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
    # Get the dimensions
    feature_size = len(feature_indices)
    cme_mtx = np.zeros((feature_size, feature_size), dtype=np.float32)

    # Prepare the index mapping for the upper-right triangle
    (i_ind, j_ind) = np.triu_indices(feature_size)

    # Populate the CME matrix
    for k in prange(len(i_ind)):
        # Get index mapping
        feature_i = i_ind[k]
        feature_j = j_ind[k]

        gene_i = feature_indices[feature_i]
        gene_j = feature_indices[feature_j]

        cme_mtx[feature_i, feature_j] = cme_mtx[feature_j, feature_i]  = compute_CME(X, gene_i, gene_j, sum_ary)
    
    return cme_mtx


# Input X: gene * cell. A row is a gene
@njit(parallel=True, fastmath=True, nopython=True)
def CME_asym_numba(X: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                    sum_ary: np.ndarray[int, np.dtype[np.float32]],
                    data_indices: np.ndarray[int, np.dtype[np.int32]],
                    feature_indices: np.ndarray[int, np.dtype[np.int32]]
                    ) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
    # Get the dimensions
    data_size = len(data_indices)
    feature_size = len(feature_indices)
    cme_mtx = np.zeros((data_size, feature_size), dtype=np.float32)

    # Populate the CME matrix
    for k in prange(data_size * feature_size):
        # Get index mapping
        i = k // feature_size
        j = k % feature_size

        data_idx = data_indices[i]
        feature_idx = feature_indices[j]

        cme_mtx[i, j] = compute_CME(X, data_idx, feature_idx, sum_ary)
            
    return cme_mtx