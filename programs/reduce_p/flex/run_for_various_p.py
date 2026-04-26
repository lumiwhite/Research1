import subprocess
import argparse
import re
import sys
import time

# python3 run_for_various_p.py --A 600 --B 800 --active 0 2 4 --nc 0,3 1,2 --timeout 100
# python3 run_for_various_p.py --A 500 --B 800 --active 0 2 4 --nc 0,3 0,1 --timeout 100

def get_unique_prime_factors_count(n):
    """ 値 n が持つ素因数の種類数を返す """
    count = 0
    d = 2
    while d * d <= n:
        if n % d == 0:
            count += 1
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        count += 1
    return count

def run_command(cmd, timeout=None):
    """ コマンドを実行し、出力をキャプチャして返す """
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if result.returncode != 0:
            return False, result.stderr
        return True, result.stdout
    except subprocess.TimeoutExpired:
        return False, "TIMEOUT"

def main():
    parser = argparse.ArgumentParser(description="指定範囲のPに対する一括探索パイプライン")
    parser.add_argument('--A', type=int, required=True, help="Pの下限 (A以上)")
    parser.add_argument('--B', type=int, required=True, help="Pの上限 (B未満)")
    parser.add_argument('--active', type=str, nargs='+', default=["0", "2", "4"], help="アクティブ行")
    parser.add_argument('--nc', type=str, nargs='*', default=[], help="非可換ペア (例: 0,3 1,2)")
    parser.add_argument('--timeout', type=int, default=60, help="1つのPに対する探索のタイムアウト(秒)")
    args = parser.parse_args()

    active_str = "".join(args.active)
    header_file = f"cycles_data_{active_str}.h"
    
    print("\n=== フェーズ 1: コンパイルと準備 ===")
    
    # 1. サイクル生成モジュールのコンパイルと実行
    print("サイクルデータ生成中...")
    run_command(["gcc", "-O3", "export_cycle.c", "-o", "export_cycle"])
    success, out = run_command(["./export_cycle", "--active"] + args.active + ["--out", header_file])
    if not success:
        print("サイクルデータの生成に失敗した。\n", out)
        sys.exit(1)
        
    # 2. 探索プログラムと検証プログラムのコンパイル
    print("Cプログラムのコンパイル中...")
    success, out = run_command([
        "gcc", "-O3", "-fopenmp", 
        f'-DCYCLES_HEADER="{header_file}"', 
        "find_fg.c", "-o", "find_fg", "-lm"
    ])
    if not success:
        print("find_fg.c のコンパイルに失敗した。\n", out)
        sys.exit(1)

    success, out = run_command([
        "gcc", "-O3", 
        f'-DCYCLES_HEADER="{header_file}"', 
        "condition_checker.c", "-o", "condition_checker", "-lm"
    ])
    if not success:
        print("condition_checker.c のコンパイルに失敗した。\n", out)
        sys.exit(1)

    print("準備完了。探索を開始する。\n")
    print("=" * 60)

    success_results = []

    # フェーズ 2: A以上B未満のPに対するループ
    for p in range(args.A, args.B):
        # 素因数がちょうど2種類(x^a * y^b)のものだけを対象とする
        if get_unique_prime_factors_count(p) != 2:
            continue
        
        print(f"[{p:4d}] 探索中... ", end="", flush=True)
        
        find_cmd = ["./find_fg", "--p", str(p), "--active"] + args.active + ["--save"]
        if args.nc:
            find_cmd += ["--nc"] + args.nc
        
        start_time = time.time()
        success, out = run_command(find_cmd, timeout=args.timeout)
        elapsed = time.time() - start_time

        if not success:
            if out == "TIMEOUT":
                print(f"❌ タイムアウト ({args.timeout}秒)")
            else:
                print("❌ 実行エラーまたは解なし")
            continue
        
        # 出力から a_vec と b_vec を正規表現で抽出
        a_match = re.search(r"a_vec:\s*\[(.*?)\]", out)
        b_match = re.search(r"b_vec:\s*\[(.*?)\]", out)

        if not a_match or not b_match:
            print(f"❌ 解が見つからなかった ({elapsed:.2f}秒)")
            continue

        a_vec = [x.strip() for x in a_match.group(1).split(',')]
        b_vec = [x.strip() for x in b_match.group(1).split(',')]

        # フェーズ 3: 見つかった解の検証
        check_cmd = ["./condition_checker", "--p", str(p), "--active"] + args.active
        if args.nc:
            check_cmd += ["--nc"] + args.nc
        check_cmd += ["--a"] + a_vec + ["--b"] + b_vec

        c_success, c_out = run_command(check_cmd)
        if c_success and "✅ 全ての条件を満たしています" in c_out:
            print(f"✅ 成功！ (探索時間: {elapsed:.2f}秒)")
            success_results.append({
                "P": p,
                "a_vec": a_vec,
                "b_vec": b_vec,
                "time": elapsed
            })
        else:
            print(f"⚠️ 解を発見したが検証に失敗した ({elapsed:.2f}秒)")
            
    # フェーズ 4: 結果の出力
    print("\n" + "=" * 60)
    print("【最終結果：すべての条件をクリアした P とパラメータ】")
    if not success_results:
        print("成功したパラメータは一つも見つからなかった。")
    else:
        for res in success_results:
            print(f"P = {res['P']} (探索時間: {res['time']:.2f}秒)")
            print(f"  a_vec: [{', '.join(res['a_vec'])}]")
            print(f"  b_vec: [{', '.join(res['b_vec'])}]")
            print("-" * 40)
    print("=" * 60)

if __name__ == "__main__":
    main()

# ============================================================
# 【最終結果：すべての条件をクリアした P とパラメータ】
# P = 144 (探索時間: 41.56秒)
#   a_vec: [85, 85, 97, 121, 121, 1, 1, 109, 91, 127, 37, 73]
#   b_vec: [126, 42, 112, 68, 92, 112, 84, 6, 87, 15, 126, 24]
# ----------------------------------------
# P = 160 (探索時間: 27.30秒)
#   a_vec: [153, 57, 33, 49, 129, 17, 41, 81, 101, 141, 1, 41]
#   b_vec: [26, 62, 80, 28, 152, 20, 90, 40, 115, 95, 80, 50]
# ----------------------------------------
# P = 192 (探索時間: 7.19秒)
#   a_vec: [113, 49, 65, 65, 1, 161, 49, 145, 25, 73, 49, 49]
#   b_vec: [102, 162, 152, 96, 0, 4, 54, 102, 87, 147, 186, 174]
# ----------------------------------------
# P = 216 (探索時間: 6.13秒)
#   a_vec: [133, 61, 37, 109, 73, 145, 109, 109, 19, 37, 163, 163]
#   b_vec: [214, 110, 30, 78, 120, 204, 54, 126, 165, 84, 27, 27]
# ----------------------------------------
# P = 288 (探索時間: 1.31秒)
#   a_vec: [269, 125, 133, 49, 169, 229, 253, 1, 169, 103, 271, 109]
#   b_vec: [222, 278, 90, 216, 84, 186, 54, 0, 108, 87, 27, 198]
# ----------------------------------------
# P = 324 (探索時間: 47.55秒)
#   a_vec: [109, 163, 1, 1, 163, 163, 1, 1, 289, 109, 1, 217]
#   b_vec: [228, 210, 261, 45, 171, 189, 90, 60, 308, 254, 30, 174]
# ----------------------------------------
# P = 384 (探索時間: 1.29秒)
#   a_vec: [265, 193, 241, 241, 97, 193, 289, 1, 49, 353, 353, 257]
#   b_vec: [279, 204, 114, 138, 252, 120, 332, 32, 106, 348, 132, 112]
# ----------------------------------------
# P = 432 (探索時間: 5.39秒)
#   a_vec: [97, 73, 1, 241, 145, 193, 73, 1, 145, 307, 109, 145]
#   b_vec: [124, 84, 24, 112, 360, 224, 264, 18, 123, 141, 234, 6]
# ----------------------------------------
# P = 448 (探索時間: 1.28秒)
#   a_vec: [435, 99, 1, 197, 29, 309, 57, 249, 293, 109, 225, 393]
#   b_vec: [273, 161, 224, 98, 14, 266, 220, 380, 162, 150, 176, 196]
# ----------------------------------------
# P = 500 (探索時間: 24.90秒)
#   a_vec: [261, 401, 201, 201, 301, 301, 1, 251, 351, 1, 251, 1]
#   b_vec: [398, 194, 340, 200, 160, 320, 275, 0, 80, 135, 225, 50]
# ----------------------------------------
# P = 576 (探索時間: 4.18秒)
#   a_vec: [415, 127, 145, 1, 145, 361, 193, 289, 337, 329, 1, 97]
#   b_vec: [375, 81, 432, 288, 0, 540, 0, 336, 104, 380, 448, 80]
# ----------------------------------------
# P = 640 (探索時間: 8.42秒)
#   a_vec: [637, 389, 129, 97, 273, 577, 41, 441, 331, 211, 441, 401]
#   b_vec: [326, 618, 384, 112, 168, 480, 580, 620, 465, 285, 460, 40]
# ----------------------------------------
# P = 648 (探索時間: 5.91秒)
#   a_vec: [37, 109, 217, 217, 433, 541, 1, 433, 271, 73, 1, 217]
#   b_vec: [250, 141, 516, 12, 528, 174, 126, 300, 273, 646, 540, 240]
# ----------------------------------------
# P = 768 (探索時間: 1.75秒)
#   a_vec: [433, 481, 65, 449, 1, 641, 433, 1, 25, 73, 1, 337]
#   b_vec: [694, 204, 680, 120, 320, 592, 630, 192, 675, 549, 144, 282]
# ----------------------------------------
# P = 784 (探索時間: 72.11秒)
#   a_vec: [393, 687, 393, 393, 1, 1, 225, 449, 225, 309, 561, 113]
#   b_vec: [427, 658, 644, 700, 532, 0, 48, 248, 142, 404, 368, 744]
# ----------------------------------------
# ============================================================