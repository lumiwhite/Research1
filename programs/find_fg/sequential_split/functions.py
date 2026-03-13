import random
import math
import numpy as np
import itertools

P = 1536
P_H = P//2
L = 12
L2 = 2*L
L_H = 12//2
J = 3

def gen_cycles(max_len, L=L, L_H=L_H):
    def is_valid(r_seq, c_seq):
        length = len(r_seq)
        for i in range(length):
            if r_seq[i] == r_seq[(i+1)%length]:
                return False
            if c_seq[i] == c_seq[(i+1)%length]:
                return False
        return True
    def gen_utcbc():
        utcbcs = set()
        r_seq = [0, 1, 2, 1]
        for c0 in range(L):
            c2 = (c0 + 1) % L_H + (L_H if c0 >= L_H else 0)
            for c1 in range(L):
                c3 = (c1 + 1) % L_H + (L_H if c1 >= L_H else 0)
                c_seq = [c0, c1, c2, c3]
                if is_valid(r_seq, c_seq):
                    positions = get_positions(r_seq, c_seq)
                    utcbcs.add(tuple(canonicalize(positions)))
        return utcbcs

    def get_positions(r_seq, c_seq):
        length = len(r_seq)
        pos = []
        for i in range(length):
            pos.append((r_seq[i], c_seq[i]))
            pos.append((r_seq[i], c_seq[(i+1)%length]))
        return pos

    def canonicalize(positions):
        symmetries = []
        curr = list(positions)
        curr_reverse = curr[::-1]
        for _ in range(len(curr)):
            curr = curr[1:] + curr[:1]
            symmetries.append(tuple(curr))
            curr_reverse = curr_reverse[1:] + curr_reverse[:1]
            symmetries.append(tuple(curr_reverse))
        return min(symmetries)
    cycles = set()
    for i in range(2, max_len//2+1):
        for r_seq in itertools.product(list(range(J)), repeat=i):
            for c_seq in itertools.product(list(range(L)), repeat=i):
                if is_valid(r_seq, c_seq):
                    positions = get_positions(r_seq, c_seq)
                    cycles.add(tuple(canonicalize(positions)))
    if max_len >= 8:
        cycles = cycles - gen_utcbc()
    return list(cycles)

def get_func_idx(cycles):
    cycles_x = []
    cycles_z = []
    for cycle in cycles:
        cycle_x = []
        cycle_z = []
        for r, c in cycle:
            if c < L_H:
                cycle_x.append((c-r) % L_H)
                cycle_z.append((r-c) % L_H + L_H)
            else:
                cycle_x.append((c-r) % L_H + L_H)
                cycle_z.append((r-c) % L_H)
        cycles_x.append(cycle_x)
        cycles_z.append(cycle_z)
    return cycles_x, cycles_z

def random_a(P):
    while True:
        a_val = random.randint(0, P-1)
        if math.gcd(a_val, P) == 1:
            return a_val

# functions.py 内の random_a を変更
# def random_a(P):
#     """
#     解空間を維持するため、a の候補を少数の「筋の良い値」に絞る。
#     (P_H = 384 の場合、a-1 が 2 や 4 などの小さな約数しか持たない値)
#     """
#     good_a_candidates = [7, 13, 25, 49, 97, 193, 257, 383] 
#     return random.choice(good_a_candidates)

def get_commute_b(a, b, c, P):
    '''
    cx+dがax+bと可換になるようなdを求める
    '''
    # 方程式を Ad ≡ B (mod P) の形に整理
    A = (a - 1) % P
    B = ((c - 1) * b) % P
    
    # 最大公約数 g = gcd(A, P) を求める
    g = math.gcd(A, P)
    # 解が存在するための必要十分条件は、B が g で割り切れることである
    if B % g != 0:
        return []
    # 方程式を g で割って簡約化する: A'd ≡ B' (mod P')
    # ここで A' = A/g, B' = B/g, P' = P/g であり、gcd(A', P') = 1 となる
    A_prime = A // g
    B_prime = B // g
    P_prime = P // g
    # A_prime と P_prime は互いに素であるため、モジュラ逆元が存在する
    # Python 3.8+ の pow(x, -1, m) を使用して特殊解 d0 を計算
    try:
        d_0 = (B_prime * pow(A_prime, -1, P_prime)) % P_prime
    except ValueError:
        # 理論上、解が存在する場合はここには到達しない
        return []
    # 全ての一般解は d = d0 + k * (P/g) (k = 0, 1, ..., g-1) となる
    # 法 P において解は正確に g 個存在する
    solutions = [(d_0 + k * P_prime) % P for k in range(g)]
    return sorted(solutions)

def get_non_commute_b(a, b, c, P):
    '''
    cx+dがax+bと非可換になるようなdを求める
    '''
    # 方程式を Ad ≡ B (mod P) の形に整理
    A = (a - 1) % P
    B = ((c - 1) * b) % P
    
    # 最大公約数 g = gcd(A, P) を求める
    g = math.gcd(A, P)
    # 解が存在するための必要十分条件は、B が g で割り切れることである
    if B % g != 0:
        return []
    # 方程式を g で割って簡約化する: A'd ≡ B' (mod P')
    # ここで A' = A/g, B' = B/g, P' = P/g であり、gcd(A', P') = 1 となる
    A_prime = A // g
    B_prime = B // g
    P_prime = P // g
    # A_prime と P_prime は互いに素であるため、モジュラ逆元が存在する
    # Python 3.8+ の pow(x, -1, m) を使用して特殊解 d0 を計算
    try:
        d_0 = (B_prime * pow(A_prime, -1, P_prime)) % P_prime
    except ValueError:
        # 理論上、解が存在する場合はここには到達しない
        return []
    # 全ての一般解は d = d0 + k * (P/g) (k = 0, 1, ..., g-1) となる
    # 法 P において解は正確に g 個存在する
    solutions = [(d_0 + k * P_prime) % P for k in range(g)]
    non_commute_solutions = set(range(P)) - set(solutions)
    # リスト形式に戻してソートする
    result_list = sorted(list(non_commute_solutions))
    return result_list

def get_new_cycles(cycles_x, cycles_z, idx):
    new_cycles_x = []
    new_cycles_z = []
    for cycle_x, cycle_z in zip(cycles_x, cycles_z):
        if idx in cycle_x and max(cycle_x) == idx:
            new_cycles_x.append(cycle_x)
        if idx in cycle_z and max(cycle_z) == idx:
            new_cycles_z.append(cycle_z)
    return new_cycles_x, new_cycles_z

def func_inv(func):
    """
    np.empty_like による未初期化メモリのバグを防ぐため、
    np.zeros で確実に0埋め初期化し、厳密な整数型を指定する。
    """
    func = np.asarray(func, dtype=int)
    inv = np.zeros(len(func), dtype=int)
    inv[func] = np.arange(len(func), dtype=int)
    return inv

def composite(left, right):
    # 両方が確実に整数配列として扱われるようにキャストして参照
    return np.asarray(left, dtype=int)[np.asarray(right, dtype=int)]

def is_closed(func, identity):
    return np.any(func == identity)

def joint_permutation(func1, func2, rule=0):
    if rule == 0:
        # 前半に func1、後半に func2 をそのまま繋げる
        return np.concatenate([func1, func2], axis=0)
        
    elif rule == 1:
        P_H = len(func1)
        res = np.empty(P_H * 2, dtype=int)
        
        # Domain 1 (偶数インデックス) を閉じる
        # func1の値 (0 から P_H-1) を偶数 (0, 2, 4...) にマッピングする
        res[0::2] = func1 * 2
        
        # Domain 2 (奇数インデックス) を閉じる
        # func2の値は (P_H から P-1) にオフセットされているため、
        # 一旦 P_H を引いて (0 から P_H-1) に戻してから、奇数 (1, 3, 5...) にマッピングする
        res[1::2] = (func2 - P_H) * 2 + 1
        
        return res

def exist_closed_cycle_x(cycles, functions):
    identity = list(range(P))
    for cycle in cycles:
        cycle_func = np.array(range(P))
        for i, idx in enumerate(cycle):
            if i%2 == 0:
                cycle_func = composite(cycle_func, functions[idx])
            else:
                cycle_func = composite(cycle_func, func_inv(functions[idx]))
        if is_closed(cycle_func, identity):
            return False
    return True
            
def exist_closed_cycle_z(cycles, functions):
    identity = list(range(P))
    for cycle in cycles:
        cycle_func = np.array(range(P))
        for i, idx in enumerate(cycle):
            if i%2 == 0:
                cycle_func = composite(cycle_func, func_inv(functions[idx]))
            else:
                cycle_func = composite(cycle_func, functions[idx])
        if is_closed(cycle_func, identity):
            return False
    return True

def generate_random_apm(a_vec, b_vec, funcs, new_cx, new_cz, idx, P=P, P_H=P_H, L=L, max_attempt=200):

    new_a_vec = a_vec.copy()
    new_b_vec = b_vec.copy()
    new_funcs = funcs.copy()
    for _ in range(max_attempt):
        # 可換性条件(条件A, B)を満たすようにAPMを生成
        a1 = random_a(P_H)
        a2 = random_a(P_H)
        valid_b1 = set(range(P_H))
        valid_b2 = set(range(P_H))
        ncom_reqs1 = []
        ncom_reqs2 = []
        if idx >= L_H:
            for i in range(L_H):
                if (idx == 8 and i == 1) or (idx == 9 and i == 0):
                    # 非可換ペアの条件は一旦保存しておく (ここでは絞り込まない)
                    ncom_reqs1.append(set(get_non_commute_b(new_a_vec[2*i], new_b_vec[2*i], a1, P_H)))
                    ncom_reqs2.append(set(get_non_commute_b(new_a_vec[2*i+1], new_b_vec[2*i+1], a2, P_H)))
                else:
                    # 可換ペアの条件はここで厳密に絞り込む (AND条件)
                    valid_b1 &= set(get_commute_b(new_a_vec[2*i], new_b_vec[2*i], a1, P_H))
                    valid_b2 &= set(get_commute_b(new_a_vec[2*i+1], new_b_vec[2*i+1], a2, P_H))
                    
            # 必須の可換条件を満たす b が存在しなければ、この a はスキップ
            if not valid_b1 or not valid_b2:
                continue
        b_cand1 = list(valid_b1)
        b_cand2 = list(valid_b2)

        valid_pairs = []
        for b1 in b_cand1:
            for b2 in b_cand2:
                is_valid_pair = True
                # 非可換の要求がある場合、b1かb2の少なくとも一方が要求リストに入っていればOK
                for req1, req2 in zip(ncom_reqs1, ncom_reqs2):
                    if b1 not in req1 and b2 not in req2:
                        is_valid_pair = False # 両方とも可換になってしまった場合は除外
                        break
                if is_valid_pair:
                    valid_pairs.append((b1, b2))
                    
        # 条件を満たすペアが1つも作れなかった場合はスキップ
        if not valid_pairs:
            continue

        random.shuffle(valid_pairs)
        test_limit = 500 if idx < L_H else 50
        for b1, b2 in valid_pairs[:test_limit]:
            # NumPy配列による関数の生成 (厳密な型指定でエラー回避)
            func1 = (int(a1) * np.arange(P_H, dtype=int) + int(b1)) % P_H
            func2 = (int(a2) * np.arange(P_H, dtype=int) + int(b2)) % P_H + P_H
            func = np.concatenate([func1, func2]).astype(int)
            
            new_funcs[idx] = func
            new_a_vec[2*idx] = a1
            new_a_vec[2*idx+1] = a2
            new_b_vec[2*idx] = b1
            new_b_vec[2*idx+1] = b2
            
            # ガース条件(C4, C6等の排除)のチェック
            if exist_closed_cycle_x(new_cx, new_funcs) and exist_closed_cycle_z(new_cz, new_funcs):
                print(f"a1={a1}, a2={a2}, b1={b1}, b2={b2} を見つけました。")
                return new_a_vec, new_b_vec, new_funcs

    return None

def generate_random_apm_debug(a_vec, b_vec, funcs, new_cx, new_cz, idx, P=P, P_H=P_H, L=L, max_attempt=200):
    new_a_vec = a_vec.copy()
    new_b_vec = b_vec.copy()
    new_funcs = funcs.copy()
    
    # --- デバッグ用のカウンター ---
    fail_commute = 0
    fail_non_commute = 0
    fail_cycle = 0
    
    for _ in range(max_attempt):
        a1 = random_a(P_H)
        a2 = random_a(P_H)
        
        valid_b1 = set(range(P_H))
        valid_b2 = set(range(P_H))
        
        ncom_reqs1 = []
        ncom_reqs2 = []
        
        if idx >= L_H:
            for i in range(L_H):
                if (idx == 8 and i == 1) or (idx == 9 and i == 0):
                    ncom_reqs1.append(set(get_non_commute_b(new_a_vec[2*i], new_b_vec[2*i], a1, P_H)))
                    ncom_reqs2.append(set(get_non_commute_b(new_a_vec[2*i+1], new_b_vec[2*i+1], a2, P_H)))
                else:
                    valid_b1 &= set(get_commute_b(new_a_vec[2*i], new_b_vec[2*i], a1, P_H))
                    valid_b2 &= set(get_commute_b(new_a_vec[2*i+1], new_b_vec[2*i+1], a2, P_H))
                    
            if not valid_b1 or not valid_b2:
                fail_commute += 1
                continue

        b_cand1 = list(valid_b1)
        b_cand2 = list(valid_b2)
        
        valid_pairs = []
        for b1 in b_cand1:
            for b2 in b_cand2:
                is_valid_pair = True
                for req1, req2 in zip(ncom_reqs1, ncom_reqs2):
                    if b1 not in req1 and b2 not in req2:
                        is_valid_pair = False 
                        break
                if is_valid_pair:
                    valid_pairs.append((b1, b2))
                    
        if not valid_pairs:
            fail_non_commute += 1
            continue

        random.shuffle(valid_pairs)
        
        # サイクル検査を通過できたかどうかのフラグ
        passed_cycle_check = False
        test_limit = 500 if idx < L_H else 50
        for b1, b2 in valid_pairs[:test_limit]:
            func1 = (int(a1) * np.arange(P_H, dtype=int) + int(b1)) % P_H
            func2 = (int(a2) * np.arange(P_H, dtype=int) + int(b2)) % P_H + P_H
            func = joint_permutation(func1, func2, 1)
            
            new_funcs[idx] = func
            new_a_vec[2*idx] = a1
            new_a_vec[2*idx+1] = a2
            new_b_vec[2*idx] = b1
            new_b_vec[2*idx+1] = b2
            
            if exist_closed_cycle_x(new_cx, new_funcs) and exist_closed_cycle_z(new_cz, new_funcs):
                print(f"a1={a1}, a2={a2}, b1={b1}, b2={b2} を見つけました。")
                return new_a_vec, new_b_vec, new_funcs
            else:
                pass # サイクル形成してしまった場合は次のペアを試す

        # テストしたペアすべてがサイクルを形成してしまった場合
        fail_cycle += 1

    # 失敗内訳の出力
    print(f" [解析] 失敗内訳: 可換枯渇={fail_commute}, 非可換失敗={fail_non_commute}, サイクル形成={fail_cycle}")
    return None