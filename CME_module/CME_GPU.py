import numpy as np
from numba import cuda as cuda
import cupy as cp
from . import config

def CME_pva_cuda(cme: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                    null_cme_ary: np.ndarray[tuple[int, int, int], np.dtype[np.float32]]
                    ) -> np.ndarray[int, np.dtype[np.float32]]:
    # Transfer data to GPU
    null_cme_ary_gpu = cp.asarray(null_cme_ary, dtype=cp.float32)
    cme_gpu = cp.asarray(cme, dtype=cp.float32)

    # Populate the null CME distribution by shuffling
    count_matrix_gpu = cp.sum(null_cme_ary_gpu < cme_gpu, axis=0)
    pval_mtx_gpu = count_matrix_gpu / float(null_cme_ary_gpu.shape[0])

    # Copy back the result and cleanup
    pval_mtx = cp.asnumpy(pval_mtx_gpu)
    del null_cme_ary_gpu, cme_gpu, count_matrix_gpu, pval_mtx_gpu

    return pval_mtx


def shuffle_and_normalize_cuda(X: np.ndarray[tuple[int, int], np.dtype[np.float32]]
                               ) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
    # Transfer data to GPU
    X_gpu = cp.asarray(X, dtype=cp.float32)

    # Shuffle and normalize
    rand_mat = cp.random.rand(*X_gpu.shape)
    shuffled_idx = rand_mat.argsort(axis=1)
    shuffled_X_gpu = cp.take_along_axis(X_gpu, shuffled_idx, axis=1)

    col_sum_ary = cp.sum(shuffled_X_gpu, axis=0, dtype=cp.float32)
    shuffled_X_gpu = shuffled_X_gpu / col_sum_ary[None, :] * 1e4

    # Copy back the result and cleanup
    shuffled_X = cp.asnumpy(shuffled_X_gpu)
    del X_gpu, rand_mat, shuffled_idx, shuffled_X_gpu, col_sum_ary

    return shuffled_X


def CME_cuda(X: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                normalize: bool = False,
                feature_indices: np.ndarray[int, np.dtype[np.int32]] = None
                ) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
    # Allocate X and sum_ary on GPU
    X_gpu = cp.asarray(X, dtype=cp.float32)
    sum_ary = cp.sum(X_gpu, axis=1, dtype=cp.float32)
    
    # Normalize the input matrix if needed
    if normalize:
        X_gpu_nan = cp.where(X_gpu == 0.0, cp.nan, X_gpu)
        medians = cp.nanmedian(X_gpu_nan, axis=1)
        X_gpu = X_gpu / medians[:, cp.newaxis]

        del X_gpu_nan, medians

    # Compute sum by gene
    sum_ary = cp.sum(X_gpu, axis=1, dtype=cp.float32)
    
    # Populate feature and data indices
    if feature_indices is None:
        feature_indices = np.arange(X.shape[0], dtype=np.int32)
    
    data_indices = np.arange(X.shape[0], dtype=np.int32)
    data_indices = data_indices[np.isin(data_indices, feature_indices, invert=True)]

    # Compute the symmetric part [upper triangle of a squire matrix]
    cme_sym = CME_sym_cuda_launcher(X_gpu, sum_ary, feature_indices)
    
    # # Compute the symmetric part [rectangular matrix]
    if len(data_indices) > 0:
        cme_asym = CME_asym_cuda_launcher(X_gpu, sum_ary, data_indices, feature_indices)
        cme = np.vstack((CME_sym, CME_asym))
    else:
        cme = cme_sym
    
    # Clean up memory
    del X_gpu, sum_ary

    return cme


@cuda.jit(device=True)
def compute_CME_cuda_device(X: cp.ndarray[tuple[int, int], cp.dtype[cp.float32]],
                            gene_i: int, gene_j: int,
                            sum_ary: cp.ndarray[int, cp.dtype[cp.float32]]
                            ) -> float:
    col_num = X.shape[1]
    min_sum = 0.0
  
    #Get the sum of the min between genes.
    for col_idx in range(col_num):
        gene_i_count = X[gene_i, col_idx]
        gene_j_count = X[gene_j, col_idx]
        min_sum += min(gene_i_count, gene_j_count)
    
    ratio_i = min_sum / sum_ary[gene_i]
    ratio_j = min_sum / sum_ary[gene_j]
    
    return 1.0 - max(ratio_i, ratio_j)


@cuda.jit
def CME_sym_cuda_kernel(X: cp.ndarray[tuple[int, int], cp.dtype[cp.float32]],
                        sum_ary: cp.ndarray[int, cp.dtype[cp.float32]],
                        feature_indices: cp.ndarray[int, cp.dtype[cp.int32]],
                        i_idx_ary: cp.ndarray[int, cp.dtype[cp.int32]],
                        j_idx_ary: cp.ndarray[int, cp.dtype[cp.int32]],
                        cme_mtx: cp.ndarray[tuple[int, int], cp.dtype[cp.float32]]
                        ) -> None:
    thread_idx = cuda.grid(1)
    
    feature_size = len(feature_indices)
    total_pairs = feature_size * (feature_size + 1) // 2
    
    if thread_idx < total_pairs:
        # Get index mapping
        feature_i = i_idx_ary[thread_idx]
        feature_j = j_idx_ary[thread_idx]
        
        gene_i = feature_indices[feature_i]
        gene_j = feature_indices[feature_j]
        
        cme_value = compute_CME_cuda_device(X, gene_i, gene_j, sum_ary)

        cme_mtx[feature_i, feature_j] = cme_mtx[feature_j, feature_i] = cme_value


def CME_sym_cuda_launcher(X: cp.ndarray[tuple[int, int], cp.dtype[cp.float32]],
                            sum_ary: cp.ndarray[int, cp.dtype[cp.float32]],
                            feature_indices: np.ndarray[int, np.dtype[np.int32]]
                            ) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
    # Get the dimensions
    feature_size = len(feature_indices)
    cme_mtx = np.zeros((feature_size, feature_size), dtype=np.float32)

    # Prepare the index mapping for the upper-right triangle
    (i_ind_ary, j_ind_ary) = np.triu_indices(feature_size)

    # Transfer meta data to GPU
    i_ind_ary_gpu = cp.asarray(i_ind_ary, dtype=cp.int32)
    j_ind_ary_gpu = cp.asarray(j_ind_ary, dtype=cp.int32)
    feature_indices_gpu = cp.asarray(feature_indices, dtype=cp.int32)
    cme_mtx_gpu = cp.zeros(cme_mtx.shape, dtype=cp.float32)

    # Set up cuda threadblock and grid sizes
    total_pairs = feature_size * (feature_size + 1) // 2
    threads_per_block = config.CUDA_config["threads_per_block"]
    blockspergrid = (total_pairs + threads_per_block - 1) // threads_per_block
    
    # Populate the CME matrix
    CME_sym_cuda_kernel[blockspergrid, threads_per_block](
        X, sum_ary, feature_indices_gpu, i_ind_ary_gpu, j_ind_ary_gpu, cme_mtx_gpu
    )
    
    # Copy back the result and clean up GPU memory
    cme_mtx = cp.asnumpy(cme_mtx_gpu)
    del i_ind_ary_gpu, j_ind_ary_gpu, feature_indices_gpu, cme_mtx_gpu

    return cme_mtx


@cuda.jit
def CME_asym_cuda_kernel(X: cp.ndarray[tuple[int, int], cp.dtype[cp.float32]],
                        sum_ary: cp.ndarray[int, cp.dtype[cp.float32]],
                        data_indices: cp.ndarray[int, cp.dtype[cp.int32]],
                        feature_indices: cp.ndarray[int, cp.dtype[cp.int32]],
                        cme_mtx: cp.ndarray[tuple[int, int], cp.dtype[cp.float32]]
                        ) -> None:
    # Get the dimensions
    data_size = len(data_indices)
    feature_size = len(feature_indices)
    
    # Obtain the cuda index
    thread_idx = cuda.grid(1)

    # Get index mapping
    row_idx = thread_idx // feature_size
    col_idx = thread_idx % feature_size
    
    if row_idx < data_size:
        data_idx = data_indices[row_idx]
        feature_idx = feature_indices[col_idx]
        
        cme_value = compute_CME_cuda_device(X, data_idx, feature_idx, sum_ary)
        cme_mtx[row_idx, col_idx] = cme_value


def CME_asym_cuda_launcher(X: cp.ndarray[tuple[int, int], cp.dtype[cp.float32]],
                            sum_ary: cp.ndarray[int, cp.dtype[cp.float32]],
                            data_indices: np.ndarray[int, np.dtype[np.int32]],
                            feature_indices: np.ndarray[int, np.dtype[np.int32]]
                            ) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
    # Get the dimensions
    data_size = len(data_indices)
    feature_size = len(feature_indices)
    cme_mtx = np.zeros((data_size, feature_size), dtype=np.float32)

    # Transfer meta data to GPU
    data_indices_gpu = cp.asarray(data_indices, dtype=cp.int32)
    feature_indices_gpu = cp.asarray(feature_indices, dtype=cp.int32)
    cme_mtx_gpu = cp.zeros((data_size, feature_size), dtype=cp.float32)

    # Set up cuda threadblock and grid sizes
    total_pairs = data_size * feature_size
    threads_per_block = config.CUDA_config["threads_per_block"]
    blockspergrid = (total_pairs + threads_per_block - 1) // threads_per_block
    
    # Populate the CME matrix
    CME_asym_cuda_kernel[blockspergrid, threads_per_block](
        X, sum_ary, data_indices_gpu, 
        feature_indices_gpu, cme_mtx_gpu
    )
    
    # Copy back the result and clean up GPU memory
    cme_mtx_gpu.copy_to_host(cme_mtx)
    del data_indices_gpu, feature_indices_gpu, cme_mtx_gpu

    return cme_mtx