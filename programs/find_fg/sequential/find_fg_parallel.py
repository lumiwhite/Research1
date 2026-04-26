import time
import random
import sys
import threading
import multiprocessing
import concurrent.futures
import numpy as np
from functions import *

def search_task(worker_id, status_dict, start_time):
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

    max_reached = 0
    reset_count = 0
    
    latent_fail_count = 0
    LATENT_STALL_LIMIT = 50 
    
    # 追加：前半（idx < 6）用の失敗カウンタ
    active_fail_count = 0
    ACTIVE_STALL_LIMIT = 50 
    
    start_time = time.time()
    while idx < L:
        # 進捗ステータス文字列の作成
        vec_display = [str(x) if x is not None else "." for x in a_vec]
        display_str = "[" + " ".join(vec_display) + "]"
        elapsed = time.time() - start_time
        
        # printする代わりに共有辞書へ自分のステータスを保存する
        status_dict[worker_id] = f"[W{worker_id:02d}] idx:{idx:2d} | Max:{max_reached:2d} | Resets:{reset_count:3d} | Vec:{display_str} | {elapsed:.1f}s"
        
        res = generate_random_apm(a_vec, b_vec, cx, cz, idx)
        if res is not None:
            a_vec[idx], b_vec[idx] = res
            idx += 1
            if idx > max_reached:
                max_reached = idx
        else:
            if idx >= L_H:
                latent_fail_count += 1
                if latent_fail_count >= LATENT_STALL_LIMIT:
                    idx = 0
                    max_reached = 0
                    reset_count += 1
                    latent_fail_count = 0
                    active_fail_count = 0
                else:
                    idx -= 1
            else:
                # 変更：前半での失敗時にもリセット制限をかける
                active_fail_count += 1
                if active_fail_count >= ACTIVE_STALL_LIMIT:
                    idx = 0
                    max_reached = 0
                    reset_count += 1
                    active_fail_count = 0
                    latent_fail_count = 0
                else:
                    idx -= 1 if idx > 0 else 0
            
            for k in range(idx, L):
                a_vec[k] = b_vec[k] = None

    return a_vec, b_vec

def display_monitor(status_dict, num_workers, stop_event):
    """
    メインプロセスで動作し、全ワーカーの状況をターミナルに描画し続けるスレッド。
    """
    print("\n" * num_workers)
    
    while not stop_event.is_set():
        sys.stdout.write(f"\033[{num_workers}A")
        
        for i in range(num_workers):
            # --- 変更点：辞書へのアクセスをtry-exceptで囲み、安全にする ---
            try:
                text = status_dict.get(i, f"[W{i:02d}] 待機中...")
            except Exception:
                # 終了処理中で辞書が破棄されていた場合はエラーを無視する
                text = f"[W{i:02d}] 終了処理中..."
                
            sys.stdout.write(f"{text}\033[K\n")
            
        sys.stdout.flush()
        time.sleep(0.1)

def parallel_search(num_workers):
    print(f"{num_workers} 個のプロセスで並列探索を開始する...")
    start_time = time.time()
    
    manager = multiprocessing.Manager()
    status_dict = manager.dict()
    
    stop_event = threading.Event()
    monitor_thread = threading.Thread(target=display_monitor, args=(status_dict, num_workers, stop_event))
    
    # デーモンスレッドに設定（メインプロセス終了時に道連れで終了させる）
    monitor_thread.daemon = True 
    monitor_thread.start()
    
    try:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(search_task, i, status_dict, start_time) for i in range(num_workers)]
            
            for future in concurrent.futures.as_completed(futures):
                a_sol, b_sol = future.result()
                if a_sol is not None:
                    stop_event.set()
                    # 他のタスクをキャンセル
                    for f in futures:
                        f.cancel()
                    
                    elapsed = time.time() - start_time
                    print(f"\n【探索成功！】合計計算時間: {elapsed:.3f} 秒")
                    return a_sol, b_sol

    except KeyboardInterrupt:
        # Ctrl+C が押された場合の強制終了処理
        stop_event.set()
        # 画面の表示崩れを防ぐため、少し下にスクロールさせる
        print("\n" * (num_workers + 2))
        print("【強制終了】ユーザー操作により探索を中断した。")
        # ワーカープロセスの完全停止を待たずに強制終了する
        sys.exit(1)
        
    except Exception as e:
        stop_event.set()
        print(f"\nプロセス実行中にエラーが発生した: {e}")
        sys.exit(1)

    stop_event.set()
    return None, None

if __name__ == "__main__":
    # 12個のワーカーで固定したい場合はここで12を指定する
    WORKERS = 12 
    
    a_result, b_result = parallel_search(num_workers=WORKERS)
    
    print("\nFinal Results:")
    print(f"a_vec: {a_result}")
    print(f"b_vec: {b_result}")