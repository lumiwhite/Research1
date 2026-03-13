import random
import math
import numpy as np
import itertools

P = 768
L = 12
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
# def random_a(P):
#     while True:
#         a_val = random.randint(0, P-1)
#         if math.gcd(a_val, P) == 1:
#             return a_val
def random_a(P):
    """
    gcd(a-1, P) が小さい a を優先的に選択する。
    P=768の場合、gcd(a-1, P) の最小値は 2。
    """
    for _ in range(1000):
        a_val = random.randint(2, P - 2)
        # a は P と互いに素、かつ a-1 が 2 以外の大きな約数を P と共有しないものを優先
        if math.gcd(a_val, P) == 1:
            if math.gcd(a_val - 1, P) <= 2: # 最適な a
                return a_val
    # 見つからない場合は通常のランダム
    while True:
        a_val = random.randint(2, P - 2)
        if math.gcd(a_val, P) == 1: return a_val

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
def func_inv(func, P=P):
    a, b = func
    a_inv = pow(a, -1, P)
    return [a_inv, (-a_inv * b) % P]
def composite(left, right, P=P):
    return [(left[0] * right[0]) % P, (left[0] * right[1] + left[1]) % P]
def is_closed(func, P=P):
    """不動点が存在するか判定"""
    a, b = func
    if a == 1: return b == 0
    d = math.gcd(a - 1, P)
    return b % d == 0
def find_closed_cycle_x(cycles, a_vec, b_vec):
    for cycle in cycles:
        cycle_func = [1, 0]
        for i, idx in enumerate(cycle):
            func = [a_vec[idx], b_vec[idx]]
            if i%2 == 0:
                cycle_func = composite(cycle_func, func)
            else:
                cycle_func = composite(cycle_func, func_inv(func))
        if is_closed(cycle_func):
            return False
    return True
            
def find_closed_cycle_z(cycles, a_vec, b_vec):
    for cycle in cycles:
        cycle_func = [1, 0]
        for i, idx in enumerate(cycle):
            func = [a_vec[idx], b_vec[idx]]
            if i % 2 == 1:
                cycle_func = composite(cycle_func, func)
            else:
                cycle_func = composite(cycle_func, func_inv(func))
        if is_closed(cycle_func):
            return False
    return True
def generate_random_apm(a_vec, b_vec, cycles_x, cycles_z, idx, P=P, L=L, max_attempt1=400, max_attempt2=200):
    new_a_vec = a_vec.copy()
    new_b_vec = b_vec.copy()
    new_cycles_x, new_cycles_z = get_new_cycles(cycles_x, cycles_z, idx)
    for _ in range(max_attempt1):
        # 可換性条件(条件A, B)を満たすようにAPMを生成
        b_cand_list = []
        a_val = random_a(P)
        if idx < L:
            b_cand = list(range(P))
        else:
            b_sets = []
            if idx >= L:
                for i in range(L):
                # 条件B: 非可換ペア (g2-f1, g3-f0)
                    if (idx == 8 and i == 1) or (idx == 9 and i == 0):
                        b_sets.append(set(get_non_commute_b(a_vec[i], b_vec[i], a_val, P)))
                    else:
                        b_sets.append(set(get_commute_b(a_vec[i], b_vec[i], a_val, P)))
            common = b_sets[0].intersection(*b_sets[1:])
            if not common: continue
            b_cand = list(common)
        for _ in range(len(b_cand)+5):
            b_val = random.choice(b_cand)
            new_a_vec[idx] = a_val
            new_b_vec[idx] = b_val
            # ガースが8になる条件(条件A, B)を満たすようにAPMを生成
            if  find_closed_cycle_x(new_cycles_x, new_a_vec, new_b_vec)and find_closed_cycle_z(new_cycles_z, new_a_vec, new_b_vec):
                print(f"a={a_val}, b={b_val}を見つけました。")
                return a_val, b_val
        # print(f"{max_attempt2}回の試行で解を見つけられませんでした。新しいaを計算します。")
    print(f"{max_attempt1}回の試行で解を見つけられませんでした。ひとつ前のインデックス{idx-1}の探索に戻ります")
    return None
        
def base_matrix(idx, L, J):
    mat_x = np.zeros((J, L), dtype=int)
    mat_z = np.zeros((J, L), dtype=int)
    for k in range(idx+1):
        for i in range(J):
            j_x = (k + i) % L_H
            j_z = (i - k) % L_H
            if k >= L_H:
                j_x += L_H
                j_z += L_H
            mat_x[i][j_x] = 1
            mat_z[i][j_z] = 1
    return mat_x, mat_z
