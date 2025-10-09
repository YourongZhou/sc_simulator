from .CME_CPU import CME_cpu, CME_pva_cpu, shuffle_and_normalize_cpu
import numpy as np
import typing
from . import config


# Only import CME_GPU if GPU is available
import GPUtil
try:
    gpus = GPUtil.getGPUs()
    if gpus:
        print("CUDA is available. CME imports GPU version.")
        for gpu in gpus:
            print(f"Found GPU: {gpu.name}")
            config.CME_config["GPU_available"] = True
    else:
        print("CUDA is not available. No NVIDIA GPUs found." \
        "CME will be imported without GPU support.")
except Exception as e:
    print(f"An error occurred: {e}")
    print("The system does not have a proper NVIDIA driver or nvidia-smi is not properly installed." \
    "CME will be imported without GPU support.")

if config.CME_config["GPU_available"]:
    from .CME_GPU import CME_cuda, CME_pva_cuda, shuffle_and_normalize_cuda


def CME(X: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                feature_indices: np.ndarray[int, np.dtype[np.int32]] = None
                ) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
    """Compute the CME matrix of a given gene-cell matrix.

    CME between two genes x and y is defined as:
        min_sum = sum(min(x_i, y_i))
        CME(x, y) = 1 - max(min_sum / sum(x), min_sum / sum(y) )
    
    Args:
        X (2D np.ndarray, type np.float32): The gene-cell matrix. Each row is a gene, each column is a cell.
        feature_indices (1D np.ndarray, type np.int32): The indices of the genes that serve as features. If None, all genes are used.

    Returns:
        cme (2D np.ndarray, type np.float32): The CME matrix bewtween genes. Could be a squire matrix or a rectangular matrix depending on the input feature_indices.

    Raises:
        N/A
    """
    # Get the config
    normalize = config.CME_config["normalize"]

    if config.CME_config["device"] == "CPU" or (not config.CME_config["GPU_available"]):
        print("Computing CME with CPU")
        return CME_cpu(X, normalize, feature_indices)
    elif config.CME_config["device"] == "GPU": 
        print("Computing CME with GPU")
        return CME_cuda(X, normalize, feature_indices)

    return None


def CME_corrected(X: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                    feature_indices: np.ndarray[tuple[int], np.dtype[np.int32]] = None):
    """Corrects the CME matrix according to the CME config, using the p-value associated with each CME value.

    Compute the p-value of each pair of genes. If the p-value is too small or too large, replace the CME value with a predefined value provided by the CME config.
    
    Args:
        X (2D np.ndarray, type np.float32): The gene-cell matrix. Each row is a gene, each column is a cell.
        feature_indices (1D np.ndarray, type np.int32): The indices of the genes that serve as features. If None, all genes are used.

    Returns:
        Corrected CME (2D np.ndarray, type np.float32): The corrected CME matrix bewtween genes, with false positive, phony CME and true positive CME values corrected according to CME config.

    Raises:
        N/A
    """

    # Get the configs
    iterations = config.CME_correction_config["iterations"]
    tp_cutoff = config.CME_correction_config["tp_cutoff"]
    fp_cutoff = config.CME_correction_config["fp_cutoff"]
    tp_val = config.CME_correction_config["tp_val"]
    fp_val = config.CME_correction_config["fp_val"]
    phony_val = config.CME_correction_config["phony_val"]
    normalize = config.CME_config["normalize"]

    cme = CME(X, feature_indices)
    null_cme_ary = rand_permute_CME(X, iterations, normalize, feature_indices)
    pval_mtx = CME_pva(cme, null_cme_ary)
    cme_corrected = CME_correction_by_pval(cme, pval_mtx, tp_cutoff, fp_cutoff, tp_val, fp_val, phony_val)

    return cme_corrected, pval_mtx


def CME_pva(cme: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                    null_cme_ary: np.ndarray[tuple[int, int, int], np.dtype[np.float32]]
                    ) -> typing.Optional[np.ndarray[int, np.dtype[np.float32]]]:
    """Compute the p-value of a CME matrix associated with a given gene-cell matrix.

    The idea is to check if the CME value between two genes increase or descrease after random shuffling of the gene-cell matrix by gene.
        * If the CME value increases significantly after shuffling, it indicates that the two genes are both markers of the same cell type.
        * If the CME value decreases significantly after shuffling, it indicates that the two genes are indeed being mutually exclusive.
        * If the CME value does not change significantly after shuffling, it indicates that the two genes are not mutually exclusive but we cannot determine if they are associated either.
    
    Args:
        cme (2D np.ndarray, type np.float32): The computed CME matrix.
        null_cme (3D np.ndarray, type np.float32): The CME matrices computed from shuffled gene-cell matrices.

    Returns:
        pval (2D np.ndarray, type np.float32): The p-value associated with each pair of genes.

    Raises:
        N/A
    """

    if config.CME_config["device"] == "CPU" or (not config.CME_config["GPU_available"]):
        print("Computing p-value with CPU")
        return CME_pva_cpu(cme, null_cme_ary)
    elif config.CME_config["device"] == "GPU":
        print("Computing p-value with GPU")
        return CME_pva_cuda(cme, null_cme_ary)
    
    return None


def shuffle_and_normalize(X: np.ndarray[tuple[int, int], np.dtype[np.float32]]
                            ) -> typing.Optional[np.ndarray[tuple[int, int], np.dtype[np.float32]]]:
    """Shuffle the gene-cell matrix by gene, and normalize the UMI in each cell to an irrelavent fixed value.
    
    Args:
        X (2D np.ndarray, type np.float32): The gene-cell matrix. Each row is a gene, each column is a cell.

    Returns:
        shuffled_X (2D np.ndarray, type np.float32): The gene-cell matrix after random shuffling.

    Raises:
        N/A
    """

    if config.CME_config["device"] == "CPU" or (not config.CME_config["GPU_available"]):
        print("Generating shuffled dummy matrices with CPU")
        return shuffle_and_normalize_cpu(X)
    elif config.CME_config["device"] == "GPU":
        print("Generating shuffled dummy matrices with GPU")
        return shuffle_and_normalize_cuda(X)
    
    return None


def CME_correction_by_pval(cme: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                            pval_mtx: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                            tp_cutoff: float = 0.95, fp_cutoff: float = 0.05,
                            tp_val: typing.Optional[float] = None,
                            fp_val: typing.Optional[float] = None,
                            phony_val: typing.Optional[float] = None
                            ) -> np.ndarray[tuple[int], np.dtype[np.float32]]:
    """Correct the CME values according to the p-value matrix and the CME config.

    Following the p-value matrix, correct the CME values according to the CME config.
    
    Args:
        cme (2D np.ndarray, type np.float32): The computed CME matrix.
        pval (2D np.ndarray, type np.float32): The p-value associated with each pair of genes.
        tp_cutoff (float): The cutoff to determine if a CME value is a true positive.
        fp_cutoff (float): The cutoff to determine if a CME value is a false positive
        tp_val (float or None): The value to replace true positive CME values. If None, true positive CME values are not corrected.
        fp_val (float or None): The value to replace false positive CME values. If None, false positive CME values are not corrected.
        phony_val (float or None): The value to replace phony CME values. If None, phony CME values are not corrected.

    Returns:
        null_cme_ary (3D np.ndarray, type np.float32): The CME matrices computed from shuffled gene-cell matrices.

    Raises:
        N/A
    """

    cme_corrected = cme.copy()
    
    # Correct the FP CME values
    if fp_val is not None:
        cme_corrected[pval_mtx < fp_cutoff] = fp_val

    # Correct the TP CME values
    if tp_val is not None:
        cme_corrected[pval_mtx > tp_cutoff] = tp_val
    
    # Correct the unrelated CME values
    if phony_val is not None:
        cme_corrected[(pval_mtx >= fp_cutoff) &
                  (pval_mtx <= tp_cutoff)] = phony_val
    
    return cme_corrected


def rand_permute_CME(X: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                        iterations: int = 50, normalize: bool = True,
                        feature_indices: np.ndarray[int, np.dtype[np.int32]] = None
                        ) -> np.ndarray[tuple[int, int, int], np.dtype[np.float32]]:
    """Gather the CME matrices after randomly shuffling the gene-cell matrix.

    Random shuffles the gene-cell matrix by gene, and computes the CME matrix for each shuffled matrix. The collection of these CME matrices is used later to compute the p-value associated with each CME value.
    
    Args:
        X (2D np.ndarray, type np.float32): The gene-cell matrix. Each row is a gene, each column is a cell.
        feature_indices (1D np.ndarray, type np.int32): The indices of the genes that serve as features. If None, all genes are used.
        iterations (int): The number of random shuffles to perform.

    Returns:
        null_cme_ary (3D np.ndarray, type np.float32): The CME matrices computed from shuffled gene-cell matrices.

    Raises:
        N/A
    """

    null_cme_ary = []

    for i in range(iterations):

        X_shuffled = shuffle_and_normalize(X)

        cme_shuffled = CME(X_shuffled, feature_indices=feature_indices)
        null_cme_ary.append(cme_shuffled)

    null_cme_ary = np.array(null_cme_ary)
    return null_cme_ary


__all__ = ["config", "CME", "CME_corrected", "CME_pva", "shuffle_and_normalize", "CME_correction_by_pval", "rand_permute_CME"]