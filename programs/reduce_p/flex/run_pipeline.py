import subprocess
import argparse
import re
import sys
#  python3 -u run_pipeline.py --p 704 --active 0 2 4 --nc 0,3 1,2
def run_command(cmd, print_output=True):
    """コマンドを実行し、出力をキャプチャして返す"""
    if print_output:
        print(f"実行中: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("【エラー】コマンドが失敗した:\n" + result.stderr)
        sys.exit(1)
    return result.stdout

def main():
    parser = argparse.ArgumentParser(description="サイクル生成・探索・検証パイプライン")
    parser.add_argument('--p', type=str, required=True, help="ブロックサイズ P")
    parser.add_argument('--active', type=str, nargs='+', default=["0", "2", "4"], help="アクティブ行")
    parser.add_argument('--nc', type=str, nargs='*', default=[], help="非可換ペア (例: 0,3 1,2)")
    args = parser.parse_args()

    active_str = "".join(args.active)
    header_file = f"cycles_data_{active_str}.h"
    
    print("\n=== フェーズ 1: サイクルデータ生成 ===")
    run_command(["gcc", "-O3", "export_cycle.c", "-o", "export_cycle"])
    run_command(["./export_cycle", "--active"] + args.active + ["--out", header_file])
    print(f"{header_file} を生成した。")

    print("\n=== フェーズ 2: 探索の実行 (find_fg) ===")
    # 修正点: f'-DCYCLES_HEADER="{header_file}"' とシンプルに括る
    run_command([
        "gcc", "-O3", "-fopenmp", 
        f'-DCYCLES_HEADER="{header_file}"', 
        "find_fg.c", "-o", "find_fg", "-lm"
    ])
    
    find_cmd = ["./find_fg", "--p", args.p, "--active"] + args.active
    if args.nc:
        find_cmd += ["--nc"] + args.nc
        
    print("探索プログラムを実行中... (時間がかかる場合がある)")
    find_output = run_command(find_cmd, print_output=False)
    print(find_output)

    # 出力から a_vec と b_vec を正規表現で抽出
    a_match = re.search(r"a_vec:\s*\[(.*?)\]", find_output)
    b_match = re.search(r"b_vec:\s*\[(.*?)\]", find_output)

    if not a_match or not b_match:
        print("【結果】解が見つからなかったため、検証フェーズはスキップする。")
        sys.exit(0)

    # "1, 2, 3" を ["1", "2", "3"] のリストに変換
    a_vec = [x.strip() for x in a_match.group(1).split(',')]
    b_vec = [x.strip() for x in b_match.group(1).split(',')]

    print("\n=== フェーズ 3: 結果の検証 (condition_checker) ===")
    # 修正点: こちらも同様に修正
    run_command([
        "gcc", "-O3", 
        f'-DCYCLES_HEADER="{header_file}"', 
        "condition_checker.c", "-o", "condition_checker", "-lm"
    ])

    check_cmd = ["./condition_checker", "--p", args.p, "--active"] + args.active
    if args.nc:
        check_cmd += ["--nc"] + args.nc
    check_cmd += ["--a"] + a_vec + ["--b"] + b_vec

    check_output = run_command(check_cmd, print_output=False)
    print(check_output)

if __name__ == "__main__":
    main()