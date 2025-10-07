from CME_CPU import CME_cpu, CME_pva_cpu, shuffle_and_normalize_cpu
from CME_GPU import CME_cuda, CME_pva_cuda, shuffle_and_normalize_cuda
import numpy as np
import typing

correction_config = {
    "iterations"    : 50,
    "tp_cutoff"     : 0.95,
    "fp_cutoff"     : 0.05,
    "fp_val"        : 0,
    "tp_val"        : None,
    "phony_val"  : 0.2
}

cme_config = {
    "device"            : "CPU",
    "normalize"         : True,
    "threads_per_block"   : 256
}

def CME(X: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                feature_indices: np.ndarray[int, np.dtype[np.int32]] = None
                ) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
    # Get the config
    normalize = cme_config["normalize"]

    if cme_config["device"] == "CPU":
        print("Computing CME with CPU")
        return CME_cpu(X, normalize, feature_indices)
    elif cme_config["device"] == "GPU": 
        print("Computing CME with GPU")
        return CME_cuda(X, normalize, feature_indices)

    return None


def CME_pva(cme: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                    null_cme_ary: np.ndarray[tuple[int, int, int], np.dtype[np.float32]]
                    ) -> np.ndarray[int, np.dtype[np.float32]]:
    if cme_config["device"] == "CPU":
        print("Computing p-value with CPU")
        return CME_pva_cpu(cme, null_cme_ary)
    elif cme_config["device"] == "GPU":
        print("Computing p-value with GPU")
        return CME_pva_cuda(cme, null_cme_ary)
    
    return None


def shuffle_and_normalize(X: np.ndarray[tuple[int, int], np.dtype[np.float32]]
                            ) -> np.ndarray[tuple[int, int], np.dtype[np.float32]]:
    if cme_config["device"] == "CPU":
        print("Generating shuffled dummy matrices with CPU")
        return shuffle_and_normalize_cpu(X)
    elif cme_config["device"] == "GPU":
        print("Generating shuffled dummy matrices with GPU")
        return shuffle_and_normalize_cuda(X)
    
    return None


def CME_corrected(X: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                    feature_indices: np.ndarray[tuple[int], np.dtype[np.int32]] = None):
    # Get the configs
    iterations = correction_config["iterations"]
    tp_cutoff = correction_config["tp_cutoff"]
    fp_cutoff = correction_config["fp_cutoff"]
    tp_val = correction_config["tp_val"]
    fp_val = correction_config["fp_val"]
    phony_val = correction_config["phony_val"]
    normalize = cme_config["normalize"]

    cme = CME(X, feature_indices)
    null_cme_ary = rand_permute_CME(X, iterations, normalize, feature_indices)
    pval_mtx = CME_pva(cme, null_cme_ary)
    cme_corrected = CME_correction_by_pval(cme, pval_mtx, tp_cutoff, fp_cutoff, tp_val, fp_val, phony_val)

    return cme_corrected, pval_mtx


def CME_correction_by_pval(cme: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                            pval_mtx: np.ndarray[tuple[int, int], np.dtype[np.float32]],
                            tp_cutoff: float = 0.95, fp_cutoff: float = 0.05,
                            tp_val: typing.Optional[float] = None,
                            fp_val: float = 0,
                            phony_val: typing.Optional[float] = None
                            ) -> np.ndarray[tuple[int], np.dtype[np.float32]]:
    cme_corrected = cme.copy()
    
    # Correct the FP CME values
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
    null_cme_ary = []

    for i in range(iterations):

        X_shuffled = shuffle_and_normalize(X)

        cme_shuffled = CME(X_shuffled, feature_indices=feature_indices)
        null_cme_ary.append(cme_shuffled)

    null_cme_ary = np.array(null_cme_ary)
    return null_cme_ary