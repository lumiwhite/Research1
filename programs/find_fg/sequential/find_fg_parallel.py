import time
import random
import multiprocessing
import concurrent.futures
import numpy as np
from functions import *

def search_task(worker_id):
    """
    1つのプロセスで実行される独立した探索タスク。
    """
    random.seed(int(time.time() * 1000) + worker_id)
    np.random.seed((int(time.time() * 1000) + worker_id) % (2**32 - 1))
    
    a_vec = [None] * L
    b_vec = [None] * L
    
    idx = 0
    cycles = gen_cycles(6)
    cx, cz = get_func_idx(cycles)

    max_reached = 0  # プロセスごとに、どこまで深く進めたかを記録する変数

    while idx < L:
        res = generate_random_apm(a_vec, b_vec, cx, cz, idx)
        
        if res is not None:
            a_vec[idx], b_vec[idx] = res
            idx += 1
            
            # 自分が到達したことのない深さまで進んだらログを出す
            if idx > max_reached:
                max_reached = idx
                print(f"[Worker {worker_id:02d}] 新記録: idx={idx} に到達 (a_vec={a_vec[:idx]})")
                
        else:
            if idx >= L_H:
                idx = 0
                max_reached = 0  # 最初からやり直すので新記録もリセット
            elif idx > 0:
                idx -= 1
            else:
                idx = 0
            
            for k in range(idx, L):
                a_vec[k] = b_vec[k] = None

    return a_vec, b_vec


def parallel_search(num_workers):
    """
    マルチプロセスで同時に探索を実行し、最初に見つかった解を返す。
    """
    print(f"{num_workers} 個のプロセスで並列探索を開始します...")
    start_time = time.time()
    
    # ProcessPoolExecutorによるマルチプロセス実行
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # 各ワーカーに一意のIDを渡してタスクを開始
        futures = [executor.submit(search_task, i) for i in range(num_workers)]
        
        # as_completed を使って、最初に完了した（解を見つけた）プロセスを取得
        for future in concurrent.futures.as_completed(futures):
            try:
                a_sol, b_sol = future.result()
                if a_sol is not None:
                    # 1つ解が見つかったら、実行中の他のプロセスをキャンセルする
                    for f in futures:
                        f.cancel()
                    
                    elapsed = time.time() - start_time
                    print(f"\n【探索成功！】計算時間: {elapsed:.3f} 秒")
                    return a_sol, b_sol
            except Exception as e:
                print(f"プロセス実行中にエラーが発生しました: {e}")

    return None, None

# ==========================================
# 実行部分
# ==========================================
if __name__ == "__main__":
    # PCの論理コア数を取得して、ワーカー数に設定する
    WORKERS = multiprocessing.cpu_count()
    
    # 並列探索の実行
    a_result, b_result = parallel_search(num_workers=WORKERS//2)  # CPUの半分をワーカーに割り当てる（必要に応じて調整）
    
    print("\nFinal Results:")
    print("a_vec:", a_result)
    print("b_vec:", b_result)

    