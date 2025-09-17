# Return environment-specific directories based on workenv (0=local, 1=cluster)
import os
import sys

def set_paths(workenv: int):
    if workenv == 1:   # cluster
        jcm_root = "/sdcc/u/lsun1/JCM_2025"
        sys.path.append(os.path.join(jcm_root, 'ThirdPartySupport', 'Python'))
        jcm_optimizer_path = "/sdcc/u/lsun1/JCMoptimizer2"
    else:              # local
        jcm_root = "/home/sun2024/JCM_2024"
        sys.path.append(os.path.join(jcm_root, 'ThirdPartySupport', 'Python'))
        jcm_optimizer_path = "/home/sun2024/JCMoptimizer2"
    return jcm_root, jcm_optimizer_path