import subprocess
import re
import matplotlib.pyplot as plt
import os

# ==========================================
# 1. シミュレーション設定 (jointbp_llr_v2用)
# ==========================================
CPP_EXEC = "./jointbp_sim"  # コンパイル後の実行ファイル名
HX_FILE = "H_X.alist"
HZ_FILE = "H_Z.alist"

# テストする物理エラー率(p)のリスト
P_VALUES = [0.03, 0.04, 0.05, 0.06, 0.07, 0.08]

# 試行回数 (N=276万規模のため、まずは少なめでテストを推奨)
TRIALS = 100

# jointbp_llr_v2 の引数仕様に合わせたリスト
# 仕様想定: executable <hx> <hz> <p> <trials> <max_iter> <osd_order>
MAX_ITER = 50
OSD_ORDER = 10

results_p = []
results_fer = []

print(f"=== SC-QLDPCシミュレーション自動実行 (N=276万規模) ===")

for p in P_VALUES:
    print(f"\n[物理エラー率 p = {p}] 実行中...")
    
    # jointbp_llr_v2 のコマンドライン引数を構築
    cmd = [
        CPP_EXEC, 
        HX_FILE, 
        HZ_FILE, 
        str(p), 
        str(TRIALS), 
        str(MAX_ITER), 
        str(OSD_ORDER)
    ]
    
    try:
        # C++シミュレータの実行
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout_data = result.stdout
        
        # jointbp_llr_v2.cpp の最終出力行をパース
        # 例: "... fer=0.00530000 iters=..." という形式を抽出
        match = re.search(r'fer=([0-9.]+)', stdout_data)
        
        if match:
            fer = float(match.group(1))
            print(f" >> 完了: FER(論理エラー率) = {fer}")
            results_p.append(p)
            results_fer.append(fer)
        else:
            print(" !! エラー: C++の出力から 'fer=' を検出できませんでした。")
            print("--- 出力ログ ---")
            print(stdout_data[-500:]) # 最後の500文字を表示
            break
            
    except subprocess.CalledProcessError as e:
        print(f" !! シミュレータが異常終了しました (Exit Code: {e.returncode})")
        print(e.stderr)
        break

# ==========================================
# 2. ウォーターフォール曲線のプロット
# ==========================================
if results_fer:
    plt.figure(figsize=(10, 7))
    
    # 統計学的に0のエラー率は微小値にして対数グラフに載せる
    plot_fer = [f if f > 0 else 1e-6 for f in results_fer]
    
    plt.plot(results_p, plot_fer, 'ro-', linewidth=2, markersize=8, label='SC-QLDPC (JointBP-OSD)')
    
    plt.yscale('log')
    plt.xlabel('Physical Error Rate (p)', fontsize=12)
    plt.ylabel('Logical Frame Error Rate (FER)', fontsize=12)
    plt.title(f'Waterfall Curve (P=768, N=2,764,800)\nTrials={TRIALS}, OSD_Order={OSD_ORDER}', fontsize=14)
    plt.grid(True, which="both", linestyle="--", alpha=0.7)
    plt.legend()
    
    # 保存
    plt.savefig("waterfall_sc_qldpc.png", dpi=300)
    print("\n>> グラフを 'waterfall_sc_qldpc.png' に保存しました。")
    plt.show()