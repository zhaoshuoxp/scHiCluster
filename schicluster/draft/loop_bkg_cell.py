# command time python /gale/ddn/snm3C/humanPFC/code/loop_bkg_cell.py --indir /gale/ddn/snm3C/humanPFC/smoothed_matrix/${res0}b_resolution/ --cell ${sample} --chrom ${c} --res ${res} --impute_mode pad2_std1_rp0.5_ws40

import sys
import h5py
import time
import cv2
cv2.useOptimized()
import argparse
import numpy as np
from scipy.sparse import save_npz, csr_matrix, coo_matrix
from sklearn.preprocessing import RobustScaler

def loop_bkg_cell(indir, cell, chrom, impute_mode, res,
            dist=10050000, cap=5, pad=5, gap=2, norm_mode='dist_trim'):

    if chrom[:3]=='chr':
        c = chrom
    else:
        c = 'chr' + chrom

    start_time = time.time()
    with h5py.File(f'{indir}{c}/{cell}_{c}_{impute_mode}.hdf5', 'r') as f:
        g = f['Matrix']
        E = csr_matrix((g['data'][()], g['indices'][()], g['indptr'][()]), g.attrs['shape']).tocoo()

    print('Load', time.time() - start_time)

    # TODO numba optimize
    start_time = time.time()
    ave, std, top, count = np.zeros((4, dist // res + pad + 1))
    for i in range(dist // res + pad + 1):
        tmp = E.diagonal(i)
        if len(tmp) > 0:
            top[i] = np.percentile(tmp, 99)
            tmp[tmp > top[i]] = top[i]
            ave[i] = np.mean(tmp)
            std[i] = np.std(tmp)
            count[i] = np.sum(tmp > 0)
        else:
            # Handle empty diagonal case
            top[i] = 0
            ave[i] = 0
            std[i] = 0 # Will be handled later to avoid div by zero
            count[i] = 0

    print('Curve', time.time() - start_time, '#Nonzero', np.sum(count))

    idx = np.triu_indices(E.shape[0], 0)
    idxfilter = ((idx[1] - idx[0]) < (dist // res + 1))
    idx = (idx[0][idxfilter], idx[1][idxfilter])
    mask = csr_matrix((np.ones(len(idx[0])), (idx[0], idx[1])), E.shape)

    start_time = time.time()

    # ================= FIX START: IndexError =================
    # 1. Calculate distances
    dist_indices = E.col - E.row
    
    # 2. Create mask for valid distances within top array bounds
    valid_mask = dist_indices < top.shape[0]
    
    # 3. Filter matrix if needed
    if not np.all(valid_mask):
        E = coo_matrix(
            (E.data[valid_mask], (E.row[valid_mask], E.col[valid_mask])),
            shape=E.shape
        )
    # ================= FIX END =================

    # Clip data by top percentile
    # E.col - E.row is now safe
    E.data = np.min([E.data, top[E.col - E.row]], axis=0)
    E = E.astype(np.float32).toarray()
    
    tmp = E[idx]
    dist_idx = idx[1] - idx[0]
    
    # ================= FIX START: Divide by Zero (Z-score) =================
    # 使用 np.errstate 忽略除以 0 的警告
    # 并在之后将结果为 NaN 或 Inf 的地方设为 0
    with np.errstate(divide='ignore', invalid='ignore'):
        # std[dist_idx] 可能包含 0
        denominator = std[dist_idx]
        numerator = tmp - ave[dist_idx]
        
        # 安全除法：如果分母为 0，结果设为 0
        z_scores = np.zeros_like(numerator)
        valid_div = denominator != 0
        z_scores[valid_div] = numerator[valid_div] / denominator[valid_div]
        
        tmp = z_scores

    # 原有的后处理逻辑保持不变
    tmp[count[dist_idx] < 100] = 0
    # tmp[std[dist_idx]==0] = 0 # 这一步已经在上面处理了
    tmp[tmp > cap] = cap
    tmp[tmp < -cap] = -cap
    E[idx] = tmp.copy()
    print('Norm', time.time() - start_time)

    start_time = time.time()
    w = pad * 2 + 1
    kernel = np.ones((w, w), np.float32)
    kernel[(pad - gap):(pad + gap + 1), (pad - gap):(pad + gap + 1)] = 0
    
    # ================= FIX START: Divide by Zero (Kernel) =================
    kernel_sum = np.sum(kernel)
    if kernel_sum > 0:
        kernel = kernel / kernel_sum
    else:
        # Fallback to avoid error, though logically shouldn't happen with default params
        kernel = kernel 
    # ================= FIX END =================

    N = cv2.filter2D(E, -1, kernel=kernel)
    
    # 再次使用 CSR 转换并保留小数位
    E = csr_matrix(np.around(E, decimals=6))
    N = csr_matrix(np.around(N, decimals=6)).multiply(mask)
    N = E - N
    print('Bkg', time.time() - start_time)

    save_npz(f'{indir}{c}/{cell}_{c}_{impute_mode}_{norm_mode}.E.npz', E)
    save_npz(f'{indir}{c}/{cell}_{c}_{impute_mode}_{norm_mode}.T.npz', N)
    return

'''
parser = argparse.ArgumentParser()
parser.add_argument('--indir', type=str, default=None, help='Directory of imputed matrix')
parser.add_argument('--cell', type=str, default=None, help='Specific identifier of a cell')
parser.add_argument('--chrom', type=str, default=None, help='Chromosome imputed')
parser.add_argument('--impute_mode', type=str, default=None, help='Suffix of imputed matrix file names')
parser.add_argument('--res', type=int, default=None, help='Bin size as integer to generate contact matrix')

parser.add_argument('--dist', type=int, default=10050000, help='Maximum distance threshold of contacts to use')
parser.add_argument('--cap', type=int, default=5, help='Trim Z-scores over the threshold')
parser.add_argument('--pad', type=int, default=5, help='One direction size of larger square for donut background')
parser.add_argument('--gap', type=int, default=2, help='One direction size of smaller square for donut background')
parser.add_argument('--norm_mode', type=str, default='dist_trim', help='Suffix of normalized file names')
opt = parser.parse_args()

loop_bkg_cell(opt.indir, opt.cell, opt.chrom, opt.impute_mode, opt.res,
        opt.dist, opt.cap, opt.pad, opt.gap, opt.norm_mode)
'''
