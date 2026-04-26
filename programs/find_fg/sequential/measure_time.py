import subprocess
import re
import statistics

def measure_average_time(executable_path, p_value, num_runs=15):
    times = []
    # Cプログラムの出力から「合計時間: X.XXX 秒」の部分を抽出する正規表現
    pattern = re.compile(r"合計時間:\s*([0-9.]+)\s*秒")

    print(f"P={p_value} で {num_runs} 回のテストを実行中...")

    for i in range(num_runs):
        print(f"実行 {i+1}/{num_runs}...", end="", flush=True)
        
        # Cプログラムの実行
        result = subprocess.run([executable_path, str(p_value)], capture_output=True, text=True)

        # 出力から実行時間を抽出
        match = pattern.search(result.stdout)
        if match:
            exec_time = float(match.group(1))
            times.append(exec_time)
            print(f" {exec_time:.3f} 秒")
        else:
            print(" 失敗（出力のパースエラー、または解が見つからなかった）")

    # 結果の集計
    if times:
        avg_time = statistics.mean(times)
        max_time = max(times)
        min_time = min(times)
        stdev = statistics.stdev(times) if len(times) > 1 else 0.0

        print("\n" + "="*30)
        print("【実行時間測定結果】")
        print(f"実行回数: {len(times)} 回")
        print(f"平均実行時間: {avg_time:.3f} 秒")
        print(f"標準偏差:   {stdev:.3f} 秒")
        print(f"最速:       {min_time:.3f} 秒")
        print(f"最遅:       {max_time:.3f} 秒")
        print("="*30)
    else:
        print("\n有効な実行結果が得られなかった。")

if __name__ == "__main__":
    # 実行ファイル名（環境に合わせて変更：Windowsの場合は "./find_fg.exe" など）
    EXEC_PATH = "./find_fg" 
    
    # 探索するPの値
    P_VAL = 768
    
    # 測定する実行回数（最初は10回程度で様子を見ることを推奨）
    NUM_RUNS = 15
    
    measure_average_time(EXEC_PATH, P_VAL, NUM_RUNS)

# P = 768 で 15 回のテストを実行中...
# ==============================
# 【実行時間測定結果】
# 実行回数: 15 回
# 平均実行時間: 160.557 秒
# 標準偏差:   191.418 秒
# 最速:       11.047 秒
# 最遅:       797.792 秒
# ==============================