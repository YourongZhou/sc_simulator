import numpy as np
from numba import cuda, njit, prange, float32, int32, int64
from multiprocessing import Pool
from multiprocessing import Process, Manager,cpu_count
import math
from numba.cuda.cudadrv.devicearray import DeviceNDArray
from sklearn.utils import compute_class_weight
from numba.cuda.random import create_xoroshiro128p_states, xoroshiro128p_uniform_float32

#numba.cuda.

threads_per_block = 256
threads_per_block_2D = (16, 16)


def CME_cuda(X, normalize=False, feature_indices=None):
    # Allocate X and sum_ary on GPU
    X_gpu = cuda.to_device(X)
    sum_ary = cuda.device_array(X.shape[0], dtype=np.float32)

    # Get the sum ary on GPU
    blockspergrid = (X.shape[0] + threads_per_block - 1) // threads_per_block
    row_sum_kernel[blockspergrid, threads_per_block](X_gpu, sum_ary)

    sum_ary_host = sum_ary.copy_to_host()
    print(sum_ary_host)
    
    # Normalize the input matrix if needed
    if normalize:
        normalize_by_row_kernel[blockspergrid, threads_per_block](X_gpu, sum_ary)
        row_sum_kernel[blockspergrid, threads_per_block](X_gpu, sum_ary)
    
    if feature_indices is None:
        feature_indices = np.arange(X.shape[0], dtype=np.int64)
    
    data_indices = np.arange(X.shape[0], dtype=np.int64)
    data_indices = data_indices[np.isin(data_indices, feature_indices, invert=True)]

    # 计算对称部分
    CME_sym = CME_sym_cuda_launcher(X_gpu, sum_ary, feature_indices)
    
    # 计算不对称部分（如果有的话）
    if len(data_indices) > 0:
        CME_asym = CME_asym_cuda_launcher(X_gpu, sum_ary, data_indices, feature_indices)
        # 合并结果
        CME_full = np.vstack((CME_sym, CME_asym))
    else:
        CME_full = CME_sym
    
    # Clean up memory
    del X_gpu, sum_ary

    return CME_full


@cuda.jit(device=True)
def compute_CME_cuda_device(X, gene_i, gene_j, sum_ary):
    col_num = X.shape[1]
    min_sum = 0.0
    
    # Get the sum of the min between genes
    for col_idx in range(col_num):
        gene_i_count = X[gene_i, col_idx]
        gene_j_count = X[gene_j, col_idx]
        min_sum += min(gene_i_count, gene_j_count)
    
    ratio_i = min_sum / sum_ary[gene_i]
    ratio_j = min_sum / sum_ary[gene_j]
    
    return 1.0 - max(ratio_i, ratio_j)


@cuda.jit
def CME_sym_cuda_kernel(X, sum_ary, feature_indices, i_idx_ary, j_idx_ary, cme_mtx):
    thread_idx = cuda.grid(1)
    
    feature_size = len(feature_indices)
    total_pairs = feature_size * (feature_size + 1) // 2
    
    if thread_idx < total_pairs:
        feature_i = i_idx_ary[thread_idx]
        feature_j = j_idx_ary[thread_idx]
        
        gene_i = feature_indices[feature_i]
        gene_j = feature_indices[feature_j]
        
        cme_value = compute_CME_cuda_device(X, gene_i, gene_j, sum_ary)
        cme_mtx[feature_i, feature_j] = cme_mtx[feature_j, feature_i] = cme_value


def CME_sym_cuda_launcher(X, sum_ary, feature_indices=None):
    if feature_indices is None:
        feature_indices = np.arange(X.shape[0], dtype=np.int64)
    
    feature_size = len(feature_indices)
    cme_mtx = np.zeros((feature_size, feature_size), dtype=np.float32)

    # Prepare the index mapping for the upper-right triangle
    (i_ind_ary, j_ind_ary) = np.triu_indices(feature_size)

    print(i_ind_ary)
    print(j_ind_ary)

    # Transfer meta data to GPU
    i_ind_ary_GPU = cuda.to_device(i_ind_ary.astype(np.int32))
    j_ind_ary_GPU = cuda.to_device(j_ind_ary.astype(np.int32))
    feature_indices_GPU = cuda.to_device(feature_indices.astype(np.int32))
    cme_mtx_GPU = cuda.to_device(cme_mtx)

    # Set up cuda threadblock and grid sizes
    total_pairs = feature_size * (feature_size + 1) // 2
    blockspergrid = (total_pairs + threads_per_block - 1) // threads_per_block
    
    # Launch the kernel
    CME_sym_cuda_kernel[blockspergrid, threads_per_block](
        X, sum_ary, feature_indices_GPU, i_ind_ary_GPU, j_ind_ary_GPU, cme_mtx_GPU
    )
    
    # Copy back the result
    cme_mtx_GPU.copy_to_host(cme_mtx)

    # Clean up GPU memory
    del i_ind_ary_GPU, j_ind_ary_GPU, feature_indices_GPU, cme_mtx_GPU

    return cme_mtx


@cuda.jit
def CME_asym_cuda_kernel(X, sum_ary, data_indices, feature_indices, cme_mtx):
    # Get the dimensions
    data_size = len(data_indices)
    feature_size = len(feature_indices)
    
    # Obtain the 2D grid indices
    thread_idx = cuda.grid(1)
    row_idx = thread_idx // feature_size
    col_idx = thread_idx % feature_size
    
    if row_idx < data_size:
        data_idx = data_indices[row_idx]
        feature_idx = feature_indices[col_idx]
        
        cme_value = compute_CME_cuda_device(X, data_idx, feature_idx, sum_ary)
        cme_mtx[row_idx, col_idx] = cme_value


def CME_asym_cuda_launcher(X, sum_ary, data_indices, feature_indices):
    data_size = len(data_indices)
    feature_size = len(feature_indices)
    cme_mtx = np.zeros((data_size, feature_size), dtype=np.float32)

    # Transfer meta data to GPU
    cme_mtx_GPU = cuda.to_device(cme_mtx)
    feature_indices_GPU = cuda.to_device(feature_indices)
    data_indices_GPU = cuda.to_device(data_indices)
    
    # Set up cuda threadblock and grid sizes
    total_pairs = data_size * feature_size
    blockspergrid = (total_pairs + threads_per_block - 1) // threads_per_block
    
    # 启动kernel
    CME_asym_cuda_kernel[blockspergrid, threads_per_block](
        X, sum_ary, data_indices_GPU, 
        feature_indices_GPU, cme_mtx_GPU
    )
    
    # Copy back the result
    cme_mtx_GPU.copy_to_host(cme_mtx)

    # Clean up GPU memory
    del data_indices_GPU, feature_indices_GPU, cme_mtx_GPU

    return cme_mtx


@cuda.jit
def row_sum_kernel(X, sum_ary):
    idx = cuda.grid(1)

    row_num = X.shape[0]
    col_num = X.shape[1]

    if idx < row_num: # Check if within bounds of rows
        current_row_sum = 0.0

        for i in range(col_num):
            current_row_sum += X[idx, i]

        sum_ary[idx] = current_row_sum


@cuda.jit
def normalize_by_row_kernel(X, sum_ary):
    row_idx = cuda.grid(1)

    row_num = X.shape[0]
    col_num = X.shape[1]

    if row_idx < row_num:
        row_sum = sum_ary[row_idx]
        if row_sum > 0:
            for col_idx in range(col_num):
                X[row_idx, col_idx] = X[row_idx, col_idx] * 1e4 / row_sum