import argparse
import numpy as np
import pandas as pd
from multiprocessing import Pool
from scipy.sparse import csr_matrix, save_npz
import os

# 全局变量
celllist = []
rs = 0
n_bins_global = 0

def load_insulation(i):
    try:
        if os.path.getsize(celllist[i]) > 0:
            data = np.load(celllist[i])
            # 简单的长度检查
            if len(data) == 0: return [i, None]
            return [i, data]
    except:
        pass
    return [i, None]

def load_boundary(i):
    try:
        # 检查文件是否存在且不为空
        if os.path.exists(celllist[i]) and os.path.getsize(celllist[i]) > 0:
            tmp = pd.read_csv(celllist[i], sep='\t', header=None)
            # 确保最后一行有数据，且坐标合理
            max_bin = int(tmp.iloc[-1, 2] // rs)
            
            # 初始化数组，长度要足够容纳该染色体
            # 注意：这里如果每个细胞文件里的最大坐标不一样，可能会有问题
            # 建议使用传入的 n_bins_global 如果有的话，或者动态调整
            local_len = max_bin + 1
            data = np.zeros(local_len)
            
            domain_rows = tmp[tmp[3] == 'domain']
            if len(domain_rows) > 0:
                coords = domain_rows[[1, 2]].values // rs
                # 边界处 +1
                # 注意边界检查
                valid_idx = coords < local_len
                # 这里的逻辑是把 domain 的起点和终点标记为边界
                for row in coords:
                    if row[0] < local_len: data[row[0]] += 1
                    if row[1] < local_len: data[row[1]] += 1
            return [i, data]
    except Exception as e:
        # print(f"Error loading {celllist[i]}: {e}")
        pass
    return [i, None]

def domain_concatcell_chr(cell_list, outprefix, res, input_type='insulation', ncpus=10):
    global celllist, rs
    celllist = np.loadtxt(cell_list, dtype=str)
    rs = res

    print(f"Loading {len(celllist)} files for {input_type}...")
    
    p = Pool(ncpus)
    if input_type == 'insulation':
        result = p.map(load_insulation, np.arange(len(celllist)))
    elif input_type == 'boundary':
        result = p.map(load_boundary, np.arange(len(celllist)))
    p.close()
    p.join()

    # 1. 确定矩阵的列数 (Max Bins)
    # 找到第一个非空的 result 来确定长度，或者取最大长度
    max_len = 0
    valid_count = 0
    for i, data in result:
        if data is not None:
            max_len = max(max_len, len(data))
            valid_count += 1
    
    print(f"Valid files loaded: {valid_count}/{len(celllist)}. Matrix width: {max_len}")

    if max_len == 0:
        print("Error: No valid data found in any cell!")
        # 生成空文件防止流程报错
        if input_type == 'insulation':
            np.save(f'{outprefix}.{input_type}.npy', np.zeros((len(celllist), 0)))
        else:
            save_npz(f'{outprefix}.{input_type}.npz', csr_matrix((len(celllist), 0)))
        return

    # 2. 构建矩阵
    # 统一长度，不足补 0
    ins = np.zeros((len(celllist), max_len), dtype=np.float32)
    
    for i, data in result:
        if data is not None:
            # 如果当前细胞的数据长度小于 max_len，后面自动补0 (np.zeros已初始化)
            # 如果大于 (理论上不应该，前面取了max)，截断
            length = min(len(data), max_len)
            ins[i, :length] = data[:length]

    # 3. 保存
    print(f"Saving to {outprefix}.{input_type}...")
    if input_type == 'insulation':
        np.save(f'{outprefix}.{input_type}.npy', ins)
    elif input_type == 'boundary':
        # Boundary 矩阵通常很稀疏，存为稀疏矩阵
        save_npz(f'{outprefix}.{input_type}.npz', csr_matrix(ins))
    
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cell_list', type=str, required=True)
    parser.add_argument('--outprefix', type=str, required=True)
    parser.add_argument('--res', type=int, default=10000)
    parser.add_argument('--input_type', type=str, default='insulation')
    parser.add_argument('--ncpus', type=int, default=10)
    opt = parser.parse_args()

    domain_concatcell_chr(opt.cell_list, opt.outprefix, opt.res, opt.input_type, opt.ncpus)
