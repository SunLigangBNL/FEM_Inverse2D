# Return environment-specific directories based on workenv (0=local, 1=cluster)
def set_paths(workenv: int):
    if workenv == 1:   # cluster
        jcm_root = "/sdcc/u/lsun1/JCM_2025"
        jcm_optimizer_path = "/sdcc/u/lsun1/JCMoptimizer2"
        #tempworkdir = "/hpcgpfs01/scratch/lsun1/tempworkdir2025"
    else:              # local
        jcm_root = "/home/sun2024/JCM_2024"
        jcm_optimizer_path = "/home/sun2024/JCMoptimizer2"
        #tempworkdir = "/home/sun2024/snap/tempworkdir2025"
    return jcm_root, jcm_optimizer_path