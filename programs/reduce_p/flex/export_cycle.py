import argparse
import itertools
import random

# ==========================================
# コマンドライン引数の設定
# ==========================================
parser = argparse.ArgumentParser(description="サイクルデータを生成します")
parser.add_argument('--active', type=int, nargs='+', default=[0, 2, 4], help="アクティブ部の行インデックスを指定 (例: 0 2 4)")
parser.add_argument('--out', type=str, default="", help="出力するC言語用ヘッダーファイル名 (省略時は自動生成)")
args = parser.parse_args()

# ==========================================
# ファイル名の自動決定
# ==========================================
if not args.out:
    # 例: [0, 2, 4] -> "024" に変換し、ファイル名に組み込む
    active_str = "".join(map(str, args.active))
    args.out = f"cycles_data_{active_str}.h"

# ==========================================
# パラメータ設定
# ==========================================
L = 12
C = args.active           # コマンドライン引数で指定された検査ノード（チェック行）
V = list(range(12))       # 変数ノード（カラム列）

def func_X(c, v):
    """ H_Xの(c, v)ブロックに対応する関数インデックスを計算する """
    if v < 6:
        return (v - c) % 6
    else:
        return 6 + ((v - 6 - c) % 6)

def func_Z(c, v):
    """ H_Zの(c, v)ブロックに対応する関数インデックスを計算する """
    if v < 6:
        return 6 + ((c - v) % 6)
    else:
        return (c - (v - 6)) % 6

def get_cycles():
    """ 完全2部グラフ K_{|C|,12} から長さ4および6のサイクルをすべて抽出する """
    raw_cycles = []
    
    # 長さ4のサイクル抽出
    for c0, c1 in itertools.combinations(C, 2):
        for v0, v1 in itertools.combinations(V, 2):
            raw_cycles.append([c0, v0, c1, v1])
            
    # 長さ6のサイクル抽出
    for v_comb in itertools.combinations(V, 3):
        v0, v1, v2 = v_comb
        for c_comb in itertools.combinations(C, 3):
            c0, c1, c2 = c_comb
            # K_3,3 部分グラフ内に存在する6つの独立な長さ6のサイクルをすべて網羅
            raw_cycles.append([c0, v0, c1, v1, c2, v2])
            raw_cycles.append([c0, v0, c1, v2, c2, v1])
            raw_cycles.append([c0, v1, c1, v0, c2, v2])
            raw_cycles.append([c0, v1, c1, v2, c2, v0])
            raw_cycles.append([c0, v2, c1, v0, c2, v1])
            raw_cycles.append([c0, v2, c1, v1, c2, v0])
        
    return raw_cycles

def gen_function_sequences(raw_cycles):
    """ 抽出したサイクルをF, Gの関数インデックスの巡回列に変換する """
    cx = []
    cz = []
    for cycle in raw_cycles:
        seq_x = []
        seq_z = []
        k = len(cycle) // 2
        for i in range(k):
            c_current = cycle[2*i]
            v_current = cycle[2*i + 1]
            c_next = cycle[(2*i + 2) % (2*k)]
            
            seq_x.append(func_X(c_current, v_current))
            seq_x.append(func_X(c_next, v_current))
            
            seq_z.append(func_Z(c_current, v_current))
            seq_z.append(func_Z(c_next, v_current))
        cx.append(seq_x)
        cz.append(seq_z)
    return cx, cz

def write_cycles_header(cx, cz, filename):
    """ C言語用のヘッダーファイルとしてエクスポートする """
    # ヘッダーガードのマクロ名もファイル名に合わせて動的に生成する
    guard_name = filename.upper().replace(".", "_")
    
    with open(filename, "w") as f:
        f.write(f"#ifndef {guard_name}\n#define {guard_name}\n\n")
        f.write("typedef struct {\n    int len;\n    int nodes[6];\n} Cycle;\n\n")
        
        for name, c_data in [("x", cx), ("z", cz)]:
            by_idx = {i: [] for i in range(L)}
            for cycle in c_data:
                by_idx[max(cycle)].append(cycle)
            
            max_c = max([len(by_idx[i]) for i in range(L)] + [1])
            
            f.write(f"const int num_cycles_{name}_by_idx[{L}] = {{")
            f.write(", ".join(str(len(by_idx[i])) for i in range(L)))
            f.write("};\n\n")
            
            f.write(f"const Cycle cycles_{name}_by_idx[{L}][{max_c}] = {{\n")
            for i in range(L):
                f.write("    {\n")
                if len(by_idx[i]) == 0:
                    f.write("        {0, {0,0,0,0,0,0}},\n")
                for cycle in by_idx[i]:
                    padded = cycle + [0] * (6 - len(cycle))
                    f.write(f"        {{{len(cycle)}, {{{', '.join(map(str, padded))}}}}},\n")
                f.write("    },\n")
            f.write("};\n\n")
            
        f.write("#endif\n")

if __name__ == "__main__":
    print(f"アクティブ行 {C} のサイクルデータを計算中...")
    
    raw_cycles = get_cycles()
    print(f"総サイクル数: {len(raw_cycles)}")
    
    cx, cz = gen_function_sequences(raw_cycles)
    write_cycles_header(cx, cz, args.out)
    print(f"--> {args.out} のエクスポートが完了した。")