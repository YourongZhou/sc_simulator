import numpy as np
from numba import cuda, njit, prange, float32, int32, int64
from multiprocessing import Pool
from multiprocessing import Process, Manager,cpu_count
import math
from numba.cuda.cudadrv.devicearray import DeviceNDArray
from sklearn.utils import compute_class_weight
from numba.cuda.random import create_xoroshiro128p_states, xoroshiro128p_uniform_float32

# Global settings
threads_per_block = 256
threads_per_block_2D = (16, 16)


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
    blocks_x = (data_size + threads_per_block_2D[0] - 1) // threads_per_block_2D[0]
    blocks_y = (feature_size + threads_per_block_2D[1] - 1) // threads_per_block_2D[1]
    
    # 启动kernel
    CME_asym_cuda_kernel[(blocks_x, blocks_y), threads_per_block_2D](
        X_device, sum_ary_device, data_indices_device, 
        feature_indices_device, cme_mtx_device
    )
    
    # 将结果拷贝回主机
    cme_mtx_device.copy_to_host(cme_mtx)
    return cme_mtx


# Input: gene-cell matrix. feature_indices
# Return:
# * matrix on GPU
# * feature_indices on GPU
# * data_indices on GPU
def CPU_to_GPU_migrate(X, feature_indices = None):

    # Get the dimensions of the matrix
    n_rows = blocks_per_grid = X.shape[0]
    n_cols = X.shape[1]

    gpu_X = cuda.to_device(X.astype(np.float32))
    sum_ary_device = cuda.device_array(n_rows, dtype=np.float32)
    row_sum_kernel[blocks_per_grid, threads_per_block](gpu_X, sum_ary_device, n_rows, n_cols)

 

    gpu_data_indices = None

    if feature_indices is None:
        feature_indices = np.arange(X.shape[0], dtype=np.int64)
    
        data_indices = np.arange(X.shape[0], dtype=np.int64)
        data_indices = data_indices[np.isin(data_indices, feature_indices, invert=True)]

    data_size = len(data_indices)
    feature_size = len(feature_indices)

   asym_cme = cuda.device_array(（n_rows, n_col), dtype=np.float32)
    sym_cme 

    # compute summary
    row_sum_kernel[blocks_per_grid, threads_per_block](X_device, sum_ary_device, n_rows, n_cols)
    


def CME_correction(X: np.ndarray, cme: np.ndarray, iterations=50,
                   normalize=False, tp_cutoff=0.95, fp_cutoff=0.05, gpu=True):

    n_rows, n_cols = X.shape

    # CPU 输出缓存
    null_CME_list = []

    # 线程配置
    blocks_per_grid = (n_rows + threads_per_block - 1) // threads_per_block

    # 为所有行分配 RNG 状态（只需 n_rows 个线程）
    rng_states = create_xoroshiro128p_states(n_rows, seed=42)

    # 循环每个 iteration
    for it in range(iterations):
        # --- Step 1: 把 X 复制到 GPU ---

        # --- Step 1: 复制 X 到 GPU ---
        d_X = cuda.to_device(X.astype(np.float32))

        # --- Step 2: GPU Shuffle + 按列归一化 ---
        shuffle_kernel[blocks_per_grid, threads_per_block](
            d_X, rng_states, n_rows, n_cols
        )
        cuda.synchronize()

        n_rows, n_cols = d_X.shape
        # 每列一个 block，每列内部使用 256 个线程
        blocks_per_grid_norm = n_cols

        normalize_by_col_kernel[blocks_per_grid_norm, threads_per_block_norm](d_X, n_rows, n_cols)
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
def shuffle_kernel(d_X, rng_states, n_rows, n_cols):
    """
    每行一个线程块: shuffle + 按列归一化
    等价于 Python: X_shuffled = X_shuffled.T / X_shuffled.T.sum(axis=1)[None, :] * 1e4
    """
    tid = cuda.blockIdx.x * cuda.blockDim.x + cuda.threadIdx.x
    if tid >= n_rows:
        return

    # 当前行
    row = d_X[tid]

    # 1. Shuffle 每行
    gpu_shuffle_row(row, n_cols, rng_states, tid)

@cuda.jit
def normalize_by_col_kernel(d_X, n_rows, n_cols):
    """
    对 device 上的矩阵 d_X 进行按列归一化
    等价于: X = X.T / X.T.sum(axis=1)[None, :] * 1e4
    """
    # 每个线程负责一列中的一部分
    col = cuda.blockIdx.x
    thread_id = cuda.threadIdx.x
    stride = cuda.blockDim.x

    # Step 1. 每列先计算总和 (分块累加 + 原子操作)
    sum_val = 0.0
    for row in range(thread_id, n_rows, stride):
        sum_val += d_X[row, col]

    # 使用共享内存做块内归约
    smem = cuda.shared.array(256, dtype=np.float32)  # 假设线程数 ≤ 256
    smem[thread_id] = sum_val
    cuda.syncthreads()

    # 并行归约
    offset = stride // 2
    while offset > 0:
        if thread_id < offset:
            smem[thread_id] += smem[thread_id + offset]
        offset //= 2
        cuda.syncthreads()

    # 最终每列的和由线程 0 保留
    col_sum = smem[0]
    cuda.syncthreads()

    # Step 2. 归一化
    if col_sum > 0.0:
        scale = 1e4 / col_sum
        for row in range(thread_id, n_rows, stride):
            d_X[row, col] *= scale


@cuda.jit
def normalize_by_row_kernel(d_X, n_rows, n_cols):
    """
    对每行按和归一化: X[row, col] = X[row, col] / sum(row)
    """
    row = cuda.blockIdx.x  # 每个 block 处理一行
    tid = cuda.threadIdx.x

    # 1. 并行计算每行的和
    partial_sum = cuda.shared.array(256, dtype=np.float32)  # 假设线程数 <= 256
    local_sum = 0.0

    # 每个线程累加一部分
    for col in range(tid, n_cols, cuda.blockDim.x):
        local_sum += d_X[row, col]
    partial_sum[tid] = local_sum
    cuda.syncthreads()

    # 2. 线程间做归约求和（Reduction）
    step = cuda.blockDim.x // 2
    while step > 0:
        if tid < step:
            partial_sum[tid] += partial_sum[tid + step]
        step //= 2
        cuda.syncthreads()

    # 3. 按和归一化
    total_sum = partial_sum[0]
    total_sum /= 1e4 
    if total_sum > 0:
        for col in range(tid, n_cols, cuda.blockDim.x):
            d_X[row, col] /= total_sum



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
    blocks_x = (data_size + threads_per_block_2D[0] - 1) // threads_per_block_2D[0]
    blocks_y = (feature_size + threads_per_block_2D[1] - 1) // threads_per_block_2D[1]
    
    # 启动kernel
    CME_asym_cuda_kernel[(blocks_x, blocks_y), threads_per_block_2D](
        X_device, sum_ary_device, data_indices_device, 
        feature_indices_device, cme_mtx_device
    )
    
    # 将结果拷贝回主机
    cme_mtx_device.copy_to_host(cme_mtx)
    return cme_mtx


# 主函数 - CUDA版本
# Input: X is gene-cell matrix that is already transferred to GPU.
def CME_cuda(X, normalize=False, feature_indices=None):
    """
    CUDA版本的完整CME计算
    """
    
    if normalize:
        blocks_per_grid = n_rows = X.shape[0]
        n_cols = X.shape[1]

        X_normed = cuda.device_array_like(X)
        X_normed.copy_to_device(X)

        normalize_by_row_kernel[blocks_per_grid, threads_per_block](X_normed, n_rows, n_cols)
        cuda.synchronize()

        # X is no longer needed
        del X

        # normalize 仍然在CPU上做
        if isinstance(X, np.ndarray):
            # median_ary = np.apply_along_axis(
            #     lambda v: np.median(v[np.nonzero(v)]), 1, X
            # )
            sum_ary = np.apply_along_axis(
                lambda v: np.sum(v[np.nonzero(v)]), 1, X
            )   
            X_normed = X * 1e4 / col_sum_ary      
        else:
            # 如果是GPU数组，先拷贝回CPU做normalize
            # X_cpu = X.copy_to_host()
            # median_ary = np.apply_along_axis(
            #     lambda v: np.median(v[np.nonzero(v)]), 1, X_cpu
            # )
            n_rows = X.shape[0]
            n_cols = X.shape[1]
            X_normed = cuda.device_array_like(X)  # 分配同样大小的显存
            X_normed.copy_to_device(X)
            blocks_per_grid = n_rows  # 每行一个 block

            # 选择按和归一化
            normalize_by_row_kernel[blocks_per_grid, threads_per_block](X_normed, n_rows, n_cols)
            cuda.synchronize()
            res = X_normed.copy_to_host()
            print(res)

        # X_normed = X / median_ary[:, None]
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
