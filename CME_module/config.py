# Configuration for CME and correction
CME_correction_config = {
    "iterations"    : 50,
    "tp_cutoff"     : 0.95,
    "fp_cutoff"     : 0.05,
    "fp_val"        : 0,
    "tp_val"        : None,
    "phony_val"  : 0.2
}

CME_config = {
    "GPU_available"     : False,
    "device"            : "CPU",
    "normalize"         : True,
}

CUDA_config = {
    "threads_per_block" : 256
}