# %%
from functions import *

# %%
# 探索の設定
a_vec = [None] * L2
b_vec = [None] * L2
functions = [None] * L
backtrack_counts = [0] * L  # 各インデックスでのバックトラック（失敗）回数を記録
RESET_LIMIT = 3

idx = 0
cycles = gen_cycles(6)
cx, cz = get_func_idx(cycles)

while idx < L:
    new_cx, new_cz = get_new_cycles(cx, cz, idx)

    res = generate_random_apm_debug(a_vec, b_vec, functions, new_cx, new_cz, idx)
    
    if res is not None:
        # 探索成功：配列全体を受け取って状態を更新する
        a_vec, b_vec, functions = res
        
        # 進展した場合は、その先の失敗カウントをリセットしておく
        if idx + 1 < L:
            backtrack_counts[idx + 1] = 0
        idx += 1
    else:
        # 探索失敗：バックトラッキング
        idx -= 1
        
        if idx < 0:
            print("\n全ての組み合わせを試行しましたが解が見つかりませんでした。")
            break
            
        # 戻った先のインデックスの失敗回数をカウント
        backtrack_counts[idx] += 1
        
        # リセット閾値に達したか判定
        if backtrack_counts[idx] >= RESET_LIMIT or idx >= L_H:
            # Index 6 以上で失敗した場合、局所的な修正では直らないため即座に全リセット
            print(f"\n【Reset】Index {idx} で致命的な手詰まり。Index 0 からやり直します。")
            idx = 0
            a_vec = [None] * L2
            b_vec = [None] * L2
            functions = [None] * L
            backtrack_counts = [0] * L
        else:
            print(f"\nIndex {idx+1} 失敗。Index {idx} へ戻ります (累積失敗: {backtrack_counts[idx]})")
            # 失敗した箇所以降のデータを正しくクリアする
            for k in range(idx + 1, L):
                a_vec[2*k] = a_vec[2*k+1] = None
                b_vec[2*k] = b_vec[2*k+1] = None
                functions[k] = None

print("\nFinal Results:")
print("a_vec:", a_vec)
print("b_vec:", b_vec)

# %%