# %%
from functions import *

# %%
# 探索の設定
a_vec = [None] * L
b_vec = [None] * L
backtrack_counts = [0] * L  # 各インデックスでのバックトラック（失敗）回数を記録
RESET_LIMIT = 3             # 同じインデックスにN回戻ってきたら0までリセット

idx = 0
cycles = gen_cycles(8)
cx, cz = get_func_idx(cycles)

while idx < L:
    print(f"\rCurrent Index: {idx} | History: {a_vec[:idx]}", end="")
    
    # 序盤のインデックス (idx < 3) は、後の自由度を確保するため
    # 1つの 'a' に固執せず、早めに切り替える（max_attempt1を小さくする）
    max_a = 50 if idx < 3 else 400
    
    res = generate_random_apm(a_vec, b_vec, cx, cz, idx, max_attempt1=max_a)
    
    if res is not None:
        # 探索成功：次のステップへ
        a_vec[idx], b_vec[idx] = res
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
        if backtrack_counts[idx] >= RESET_LIMIT:
            print(f"\n【Reset】Index {idx} で停滞（{backtrack_counts[idx]}回目）。Index 0 からやり直します。")
            idx = 0
            a_vec = [None] * L
            b_vec = [None] * L
            backtrack_counts = [0] * L # 全てのカウンタをリセット
        else:
            print(f"\nIndex {idx+1} 失敗。Index {idx} へ戻ります (累積失敗: {backtrack_counts[idx]})")
            # 失敗した箇所以降をクリア
            for k in range(idx + 1, L):
                a_vec[k] = b_vec[k] = None

print("\nFinal Results:")
print("a_vec:", a_vec)
print("b_vec:", b_vec)

# %%



