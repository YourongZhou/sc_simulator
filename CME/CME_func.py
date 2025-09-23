import numpy as np
from numba import cuda, njit, prange, float32, int32, int64
from multiprocessing import Pool
from multiprocessing import Process, Manager,cpu_count
import math
from numba.cuda.cudadrv.devicearray import DeviceNDArray
from sklearn.utils import compute_class_weight
from numba.cuda.random import create_xoroshiro128p_states, xoroshiro128p_uniform_float32

def CME(X, normalize=False, feature_indices=None, gpu=False):

    if not isinstance(X, DeviceNDArray):
        X = X.astype(np.float32)

    if gpu:
        return CME_cuda(X, normalize, feature_indices)
    else:
        return CME_cpu(X, normalize, feature_indices)

### correction function via cpu ###
def CME_correction(X: np.ndarray[np.float32], cme: np.ndarray[np.float32], iterations=50, normalize=True, feature_indices=None, tp_cutoff=0.95, fp_cutoff=0.05, gpu=False) -> np.ndarray[np.float32]:
    # Populate the null CME distribution by shuffling
    null_CME = []

    X_shuffled = X.copy()

    for i in range(iterations):

        for row in X_shuffled:
            np.random.shuffle(row)

        X_shuffled = X_shuffled / X_shuffled.sum(axis=1)[None,:].T * 1e4
        
        cme_shuffled = CME(X_shuffled, normalize=normalize, feature_indices=feature_indices, gpu=gpu)
        null_CME.append(cme_shuffled)

    null_CME_ary = np.array(null_CME)

    cme_rank = np.zeros(cme.shape)
    for i in range(cme.shape[0]):
        for j in range(cme.shape[1]):

            cme_rank[i, j] = np.sum(null_CME_ary[:, i, j] < cme[i, j]) / float(iterations)

    cme_corrected = cme.copy()
    cme_corrected[cme_rank < fp_cutoff] = 0
    cme_corrected[(cme_rank >= fp_cutoff) &
                  (cme_rank <= tp_cutoff)] = 0.2
    
    return cme_corrected
### correction function via cpu ###




### correction function via cuda ###
@cuda.jit(device=True, inline=True)
def gpu_shuffle_row(row, n_cols, rng_states, thread_id):
    """
    Fisher-Yates shuffle for a single row on GPU.
    """
    for i in range(n_cols - 1, 0, -1):
        r = xoroshiro128p_uniform_float32(rng_states, thread_id)
        j = int(r * (i + 1))
        if j != i:
            tmp = row[i]
            row[i] = row[j]
            row[j] = tmp


@cuda.jit
def shuffle_and_normalize_kernel(d_X, rng_states, n_rows, n_cols):
    """
    每行一个线程：shuffle + normalize
    d_X: [n_rows, n_cols]
    """
    tid = cuda.blockIdx.x * cuda.blockDim.x + cuda.threadIdx.x
    if tid >= n_rows:
        return

    row = d_X[tid]

    # 1. Shuffle
    gpu_shuffle_row(row, n_cols, rng_states, tid)

    # 2. Normalize
    row_sum = 0.0
    for i in range(n_cols):
        row_sum += row[i]

    if row_sum > 0.0:
        scale = 1e4 / row_sum
        for i in range(n_cols):
            row[i] *= scale


def CME_correction(X: np.ndarray, cme: np.ndarray, iterations=50,
                   normalize=False, tp_cutoff=0.95, fp_cutoff=0.05, gpu=True):

    n_rows, n_cols = X.shape

    # CPU 输出缓存
    null_CME_list = []

    # 线程配置
    threads_per_block = 256
    blocks_per_grid = (n_rows + threads_per_block - 1) // threads_per_block

    # 为所有行分配 RNG 状态（只需 n_rows 个线程）
    rng_states = create_xoroshiro128p_states(n_rows, seed=42)

    # 循环每个 iteration
    for it in range(iterations):
        # --- Step 1: 把 X 复制到 GPU ---

        d_X = cuda.to_device(X.astype(np.float32))
        # --- Step 2: GPU Shuffle + Normalize ---

        shuffle_and_normalize_kernel[blocks_per_grid, threads_per_block](
            d_X, rng_states, n_rows, n_cols
        )

        cuda.synchronize()

        # --- Step 3: 拷贝回 CPU ---
        # X_shuffled = d_X.copy_to_host()

        # --- Step 4: 计算 CME ---
        cme_shuffled = CME(d_X, normalize=normalize, gpu=True)
        null_CME_list.append(cme_shuffled)

        # --- Step 5: 清理显存，防止累积 ---
        del d_X
        cuda.current_context().deallocations.clear()

    # Stack 成 numpy array
    null_CME = np.array(null_CME_list)

    # === CME rank correction ===
    cme_rank = np.mean(null_CME < cme[np.newaxis, :, :], axis=0)

    cme_corrected = cme.copy()
    cme_corrected[cme_rank < fp_cutoff] = 0
    cme_corrected[(cme_rank >= fp_cutoff) & (cme_rank <= tp_cutoff)] = 0.2

    return cme_corrected
### correction function via cuda ###



##################### CPU ver BEGIN ######################
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
        feature_indices = np.arange(X.shape[0], dtype=np.int32)

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
def CME_cpu(X: np.ndarray[np.float32], normalize=True, feature_indices=None) -> np.ndarray[np.float32]:

    #X = X.astype(X.shape, dtype=np.float32)

    if normalize:
        median_ary = np.apply_along_axis(lambda v: np.median(v[np.nonzero(v)]), 1, X)
        X_normed = X / median_ary[:, None]
    else:
        X_normed = X

    if feature_indices is None:
        feature_indices = np.arange(X.shape[0], dtype=np.int32)

    data_indices = np.arange(X.shape[0], dtype=np.int32)
    data_indices = data_indices[np.isin(data_indices, feature_indices, invert=True)]

    # Compute CME for symmetric part
    CME = CME_sym_numba(X_normed, feature_indices)

    # Compute CME for asymmetric part, if there is any.
    if len(data_indices) > 0:
        CME_asym = CME_asym_numba(X_normed, data_indices, feature_indices)

        CME = np.vstack((CME, CME_asym))

    return CME
##################### CPU ver END ######################


##################### CUDA BEGIN ######################
# CUDA版本的工具函数

@cuda.jit
def row_sum_kernel(X, sum_ary, n_rows, n_cols):
    """
    每行一个 block，每列一个 thread，按行求和
    """
    row = cuda.blockIdx.x
    tid = cuda.threadIdx.x

    # 正确：共享内存用 np.float32
    smem = cuda.shared.array(shape=1024, dtype=np.float32)

    # 每个线程累加一部分列
    temp = 0.0
    for col in range(tid, n_cols, cuda.blockDim.x):
        temp += X[row, col]

    smem[tid] = temp
    cuda.syncthreads()

    # 归约（Reduction）
    stride = cuda.blockDim.x // 2
    while stride > 0:
        if tid < stride:
            smem[tid] += smem[tid + stride]
        stride //= 2
        cuda.syncthreads()

    if tid == 0:
        sum_ary[row] = smem[0]


@cuda.jit(device=True)
def compute_CME_cuda_device(X, i, j, sum_ary):
    """
    CUDA设备端的CME计算函数
    """
    n_features = X.shape[1]
    min_sum = 0.0
    
    # 计算min(X[i,:], X[j,:])的和
    for k in range(n_features):
        val_i = X[i, k]
        val_j = X[j, k]
        min_sum += min(val_i, val_j)
    
    ratio_x = min_sum / sum_ary[i]
    ratio_y = min_sum / sum_ary[j]
    
    return 1.0 - max(ratio_x, ratio_y)

# 对称矩阵计算的CUDA版本
@cuda.jit
def CME_sym_cuda_kernel(X, sum_ary, feature_indices, cme_mtx):
    """
    对称CME矩阵计算的CUDA kernel
    """
    # 获取当前线程的全局索引
    idx = cuda.grid(1)
    
    n_features = X.shape[0]
    feature_size = len(feature_indices)
    total_pairs = feature_size * (feature_size + 1) // 2
    
    if idx < total_pairs:
        # 将线性索引转换为(i,j)坐标
        # 使用三角形索引计算
        i = 0
        j = idx
        while j > i:
            j -= i + 1
            i += 1
        
        x_idx = feature_indices[i]
        y_idx = feature_indices[j]
        
        cme_value = compute_CME_cuda_device(X, x_idx, y_idx, sum_ary)
        cme_mtx[i, j] = cme_value
        cme_mtx[j, i] = cme_value

def CME_sym_cuda(X, feature_indices=None):
    """
    CUDA版本的对称CME矩阵计算
    """
    if feature_indices is None:
        feature_indices = np.arange(X.shape[0], dtype=np.int64)
    
    feature_size = len(feature_indices)
    cme_mtx = np.zeros((feature_size, feature_size), dtype=np.float32)
    
    # 将数据转移到GPU
    if isinstance(X, DeviceNDArray):
        X_device = X
        n_rows = X.shape[0]
        n_cols = X.shape[1]
        # sum_ary = X.sum(axis=1).astype(np.float32)
        # sum_ary_device = cuda.to_device(sum_ary)
        sum_ary_device = cuda.device_array(n_rows, dtype=np.float32)
        threads_per_block = 32
        blocks_per_grid = n_rows  # 每行一个 block
        row_sum_kernel[blocks_per_grid, threads_per_block](X_device, sum_ary_device, n_rows, n_cols)
        feature_indices_device = cuda.to_device(feature_indices.astype(np.int32))
        cme_mtx_device = cuda.to_device(cme_mtx)
    else:
        X_device = cuda.to_device(X.astype(np.float32))
        sum_ary = X.sum(axis=1).astype(np.float32)
        sum_ary_device = cuda.to_device(sum_ary)
        feature_indices_device = cuda.to_device(feature_indices.astype(np.int32))
        cme_mtx_device = cuda.to_device(cme_mtx)
    
    # 计算总对数
    total_pairs = feature_size * (feature_size + 1) // 2
    
    # 配置CUDA线程
    threadsperblock = 256
    blockspergrid = (total_pairs + threadsperblock - 1) // threadsperblock
    
    # 启动kernel
    CME_sym_cuda_kernel[blockspergrid, threadsperblock](
        X_device, sum_ary_device, feature_indices_device, cme_mtx_device
    )
    
    # 将结果拷贝回主机
    cme_mtx_device.copy_to_host(cme_mtx)
    return cme_mtx

# 不对称矩阵计算的CUDA版本
@cuda.jit
def CME_asym_cuda_kernel(X, sum_ary, data_indices, feature_indices, cme_mtx):
    """
    不对称CME矩阵计算的CUDA kernel
    """
    # 获取2D网格中的位置
    i = cuda.blockIdx.x * cuda.blockDim.x + cuda.threadIdx.x
    j = cuda.blockIdx.y * cuda.blockDim.y + cuda.threadIdx.y
    
    data_size = len(data_indices)
    feature_size = len(feature_indices)
    
    if i < data_size and j < feature_size:
        data_idx = data_indices[i]
        feature_idx = feature_indices[j]
        
        cme_value = compute_CME_cuda_device(X, data_idx, feature_idx, sum_ary)
        cme_mtx[i, j] = cme_value

def CME_asym_cuda(X, data_indices, feature_indices):
    """
    CUDA版本的不对称CME矩阵计算
    """
    data_size = len(data_indices)
    feature_size = len(feature_indices)
    cme_mtx = np.zeros((data_size, feature_size), dtype=np.float32)
    
    if isinstance(X, DeviceNDArray):
        X_device = X
        n_rows = X.shape[0]
        n_cols = X.shape[1]
        # sum_ary = X.sum(axis=1).astype(np.float32)
        # sum_ary_device = cuda.to_device(sum_ary)
        sum_ary_device = cuda.device_array(n_rows, dtype=np.float32)
        threads_per_block = 32
        blocks_per_grid = n_rows  # 每行一个 block
        row_sum_kernel[blocks_per_grid, threads_per_block](X_device, sum_ary_device, n_rows, n_cols)
        feature_indices_device = cuda.to_device(feature_indices.astype(np.int32))
        cme_mtx_device = cuda.to_device(cme_mtx)

    else:
        # 将数据转移到GPU
        X_device = cuda.to_device(X.astype(np.float32))
        sum_ary = X.sum(axis=1).astype(np.float32)
        sum_ary_device = cuda.to_device(sum_ary)
        data_indices_device = cuda.to_device(data_indices.astype(np.int32))
        feature_indices_device = cuda.to_device(feature_indices.astype(np.int32))
        cme_mtx_device = cuda.to_device(cme_mtx)
    
    # 配置2D网格
    threads_per_block = (16, 16)  # 256 threads per block
    blocks_x = (data_size + threads_per_block[0] - 1) // threads_per_block[0]
    blocks_y = (feature_size + threads_per_block[1] - 1) // threads_per_block[1]
    
    # 启动kernel
    CME_asym_cuda_kernel[(blocks_x, blocks_y), threads_per_block](
        X_device, sum_ary_device, data_indices_device, 
        feature_indices_device, cme_mtx_device
    )
    
    # 将结果拷贝回主机
    cme_mtx_device.copy_to_host(cme_mtx)
    return cme_mtx

# 优化的对称矩阵计算（使用共享内存）
@cuda.jit
def CME_sym_cuda_optimized_kernel(X, sum_ary, feature_indices, cme_mtx):
    """
    使用共享内存优化的对称CME矩阵计算
    """
    idx = cuda.grid(1)
    feature_size = len(feature_indices)
    total_pairs = feature_size * (feature_size + 1) // 2
    
    if idx < total_pairs:
        # 计算三角形索引
        i = int((math.sqrt(8 * idx + 1) - 1) / 2)
        j = idx - i * (i + 1) // 2
        
        x_idx = feature_indices[i]
        y_idx = feature_indices[j]
        
        cme_value = compute_CME_cuda_device(X, x_idx, y_idx, sum_ary)
        cme_mtx[i, j] = cme_value
        cme_mtx[j, i] = cme_value

# 主函数 - CUDA版本
def CME_cuda(X, normalize=False, feature_indices=None):
    """
    CUDA版本的完整CME计算
    """
    
    if normalize:
        # normalize 仍然在CPU上做
        if isinstance(X, np.ndarray):
            median_ary = np.apply_along_axis(
                lambda v: np.median(v[np.nonzero(v)]), 1, X
            )
        else:
            # 如果是GPU数组，先拷贝回CPU做normalize
            X_cpu = X.copy_to_host()
            median_ary = np.apply_along_axis(
                lambda v: np.median(v[np.nonzero(v)]), 1, X_cpu
            )
        X_normed = X / median_ary[:, None]
    else:
        X_normed = X
    # if normalize:
    #     # 在CPU上计算中位数（这个操作不适合GPU）
    #     median_ary = np.apply_along_axis(lambda v: np.median(v[np.nonzero(v)]), 1, X)
    #     X_normed = X / median_ary[:, None]
    # else:
    #     X_normed = X
    
    if feature_indices is None:
        feature_indices = np.arange(X.shape[0], dtype=np.int64)
    
    data_indices = np.arange(X.shape[0], dtype=np.int64)
    data_indices = data_indices[np.isin(data_indices, feature_indices, invert=True)]
    
    # 计算对称部分
    CME_sym = CME_sym_cuda(X_normed, feature_indices)
    
    # 计算不对称部分（如果有的话）
    if len(data_indices) > 0:
        CME_asym = CME_asym_cuda(X_normed, data_indices, feature_indices)
        # 合并结果
        CME_full = np.vstack((CME_sym, CME_asym))
    else:
        CME_full = CME_sym
    
    return CME_full
##################### CUDA END ######################
