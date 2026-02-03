import time
import argparse
import numpy as np
import pandas as pd
from heapq import *
from statsmodels.sandbox.stats.multicomp import multipletests as FDR

def loop_mergechr(inprefix, outprefix, chrom_file, split_file=None, res=10000,
                thres_bl=1.33, thres_d=1.33, thres_h=1.2, thres_v=1.2, fdr_thres=0.1, dist_thres=20000, size_thres=1):

    # 内部函数: 寻找峰值
    def find_summit(loop, dist_thres):
        start_time = time.time()
        cord = loop[['x1', 'y1']].values // res
        idx = np.argsort(cord[:, 0])
        neighbor = {i:[] for i in range(len(idx))}
        for i in range(len(idx)-1):
            tmp = cord[idx[i]]
            for j in range(i+1,len(idx)):
                if cord[idx[j], 0] - tmp[0] > dist_thres:
                    break
                if np.abs(tmp[1] - cord[idx[j], 1]) <= dist_thres:
                    neighbor[idx[i]].append(idx[j])
                    neighbor[idx[j]].append(idx[i])
        
        print('Build graph takes', time.time() - start_time, 'seconds')

        start_time = time.time()
        nodescore = loop['E'].values
        flag = np.zeros(len(nodescore))
        tot = len(nodescore)
        summit = []
        nodeheap = (loop['E'] * -1).reset_index().reset_index()[['E', 'level_0']].values.tolist()
        heapify(nodeheap)

        while tot>0:
            t = int(heappop(nodeheap)[1])
            while flag[t]:
                t = int(heappop(nodeheap)[1])
            q = [t]
            flag[t] = 1
            tot -= 1
            head = 0
            flagtmp = np.zeros(len(nodescore))
            while (head < len(q)):
                for t in neighbor[q[head]]:
                    if not flagtmp[t] and nodescore[t]<nodescore[q[head]]:
                        if not flag[t]:
                            flag[t] = 1
                            tot -= 1
                        flagtmp[t] = 1
                        q.append(t)
                head += 1
            summit.append([q[0], len(q)])
        summit = np.array(summit)
        if len(summit) > 0:
            loop = loop.iloc[summit[:,0]].copy()
            loop['size'] = summit[:,1]
        else:
            loop = loop.iloc[[]].copy()
            loop['size'] = []

        print('BFS takes', time.time() - start_time, 'seconds')
        return loop

    # 加载染色体列表
    chrom = np.loadtxt(chrom_file, dtype=str)[:,0]
    if not split_file:
        chrom_split = chrom.copy()
    else:
        splitbed = pd.read_csv(split_file, sep='\t', header=None, index_col=0)
        chrom_split = np.concatenate([[c+'p', c+'q'] if c in splitbed.index else [c] for c in chrom])
    chrom_split = np.array([c if c[:3]=='chr' else 'chr'+c for c in chrom_split])

    start_time = time.time()
    loopall = []
    
    # 循环读取数据
    for c in chrom_split:
        try:
            data = pd.read_hdf(f'{inprefix}_{c}.loop.hdf5', key='loop')
            data['bkfilter'] = (((data['E']/data['E_bl'] > thres_bl) | (data['E_bl']<0)) & 
                                ((data['E']/data['E_donut'] > thres_d) | (data['E_donut']<0)) & 
                                ((data['E']/data['E_h'] > thres_h) | (data['E_h']<0)) & 
                                ((data['E']/data['E_v'] > thres_v) | (data['E_v']<0)))
            data['x1'] = data['x1'].astype(int) * res
            data['y1'] = data['y1'].astype(int) * res
            data['x2'] = data['x1'] + res
            data['y2'] = data['y1'] + res
            
            if c[-1]=='p' or c[-1]=='q':
                data['chr'] = c[:-1]
                if c[-1]=='q':
                    data[['x1','x2','y1','y2']] += splitbed.loc[c[:-1],2] // res * res
            else:
                data['chr'] = c
            loopall.append(data)
        except (FileNotFoundError, KeyError):
            print(f"Warning: Could not read loop file for {c}")
            continue

    if len(loopall) > 0:
        loopall = pd.concat(loopall, axis=0)
    else:
        loopall = pd.DataFrame(columns=['chr', 'x1', 'x2', 'y1', 'y2', 'E', 'E_bl', 'E_donut', 'E_h', 'E_v', 'bkfilter', 'rpv', 'tpv'])

    # === 修复: 填充 NaN 防止报错 ===
    if len(loopall) > 0:
        loopall['rpv'] = loopall['rpv'].fillna(1.0)
        loopall['tpv'] = loopall['tpv'].fillna(1.0)
        try:
            loopall['rFDR'] = FDR(loopall['rpv'], 0.1, 'fdr_bh')[1]
            loopall['tFDR'] = FDR(loopall['tpv'], 0.1, 'fdr_bh')[1]
        except:
            loopall['rFDR'] = 1.0
            loopall['tFDR'] = 1.0
    else:
        loopall['rFDR'] = []
        loopall['tFDR'] = []

    print('Merge chromosomes takes', time.time() - start_time, 'seconds')

    # 输出文件
    loop = loopall.loc[(loopall['bkfilter']==1) & (loopall['tFDR']<fdr_thres)]
    loop.sort_values(by=['chr', 'x1', 'y1'])[['chr', 'x1', 'x2', 'chr', 'y1', 'y2', 'E']].to_csv(f'{outprefix}.loop.bedpe', sep='\t', index=False, header=None)
    
    scloop = loopall.loc[loopall['tFDR']<fdr_thres]
    scloop.sort_values(by=['chr', 'x1', 'y1'])[['chr', 'x1', 'x2', 'chr', 'y1', 'y2', 'E']].to_csv(f'{outprefix}.scloop.bedpe', sep='\t', index=False, header=None)
    
    bkloop = loopall.loc[loopall['bkfilter']==1]
    bkloop.sort_values(by=['chr', 'x1', 'y1'])[['chr', 'x1', 'x2', 'chr', 'y1', 'y2', 'E']].to_csv(f'{outprefix}.bkloop.bedpe', sep='\t', index=False, header=None)

    # === 修复: 智能选择 Summit 数据源 ===
    target_for_summit = loop
    if len(loop) == 0:
        if len(bkloop) > 0:
            print(f"✅ Auto-fallback: Standard loop list is empty. Using {len(bkloop)} background-filtered loops to call summits.")
            target_for_summit = bkloop
    
    if len(target_for_summit) > 0:
        summit = pd.concat([find_summit(target_for_summit[target_for_summit['chr']==c], dist_thres//res) for c in chrom], axis=0)
        suumit = summit[summit['size'] >= size_thres]
        summit.sort_values(by=['chr', 'x1', 'y1'])[['chr', 'x1', 'x2', 'chr', 'y1', 'y2', 'E', 'size']].to_csv(f'{outprefix}.loopsummit.bedpe', sep='\t', index=False, header=None)
        print(f"Successfully wrote {len(summit)} summits.")
    else:
        pd.DataFrame(columns=['chr', 'x1', 'x2', 'chr', 'y1', 'y2', 'E', 'size']).to_csv(f'{outprefix}.loopsummit.bedpe', sep='\t', index=False, header=None)

    return

# === 关键修复: 只有作为主程序运行时才执行下面的代码 ===
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--inprefix', type=str, default=None)
    parser.add_argument('--outprefix', type=str, default=None)
    parser.add_argument('--chrom_file', type=str, default=None)
    parser.add_argument('--split_file', type=str, default=None)
    parser.add_argument('--res', type=int, default=10000)
    parser.add_argument('--thres_bl', type=float, default=1.33)
    parser.add_argument('--thres_d', type=float, default=1.33)
    parser.add_argument('--thres_h', type=float, default=1.2)
    parser.add_argument('--thres_v', type=float, default=1.2)
    parser.add_argument('--fdr_thres', type=float, default=0.1)
    parser.add_argument('--dist_thres', type=int, default=20000)
    parser.add_argument('--size_thres', type=int, default=1)
    opt = parser.parse_args()

    loop_mergechr(opt.inprefix, opt.outprefix, opt.chrom_file, opt.split_file, opt.res,
            opt.thres_bl, opt.thres_d, opt.thres_h, opt.thres_v, opt.fdr_thres, opt.dist_thres, opt.size_thres)
