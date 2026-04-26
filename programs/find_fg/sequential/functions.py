import random
import math
import numpy as np
import itertools

P = 768
L = 12
L_H = 12 // 2
J = 3
ACTIVE_S = [0, 2, 4]

# 【追加】aの候補を事前計算し、グローバルリストとして保持する
VALID_A_LIST = [a for a in range(1, P) if math.gcd(a, P) == 1]

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

def get_commute_b(a, b, c, P):
    A = (a - 1) % P
    B = ((c - 1) * b) % P
    g = math.gcd(A, P)
    if B % g != 0:
        return []
    A_prime = A // g
    B_prime = B // g
    P_prime = P // g
    try:
        d_0 = (B_prime * pow(A_prime, -1, P_prime)) % P_prime
    except ValueError:
        return []
    solutions = [(d_0 + k * P_prime) % P for k in range(g)]
    return sorted(solutions)

# 探索の順番を変えたら変更する必要がある
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
    a, b = func
    if a == 1: return b == 0
    d = math.gcd(a - 1, P)
    return b % d == 0

def is_commute(func1, func2, P=P):
    left = composite(func1, func2, P)
    right = composite(func2, func1, P)
    return left == right

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
            if i%2 == 1:
                cycle_func = composite(cycle_func, func)
            else:
                cycle_func = composite(cycle_func, func_inv(func))
        if is_closed(cycle_func):
            return False
    return True

def generate_random_apm(a_vec, b_vec, cycles_x, cycles_z, idx, P=P, L=L):
    new_a_vec = a_vec.copy()
    new_b_vec = b_vec.copy()
    new_cycles_x, new_cycles_z = get_new_cycles(cycles_x, cycles_z, idx)
    a_count = 0
    # 事前計算したaのリストをコピーしてシャッフル（ランダム性を確保しつつ全探索）
    a_candidates = VALID_A_LIST.copy()
    random.shuffle(a_candidates)
    
    for a_val in a_candidates:
        a_count += 1
        b_cand = []
        
        if idx < L_H:
            b_cand = list(range(P))
          
        else:
            if idx == 8 or idx == 9:
                base_i = 9 - idx
                initial_cands = list(set(range(P)) - set(get_commute_b(a_vec[base_i], b_vec[base_i], a_val, P)))
                if not initial_cands:
                    continue
            else:
                initial_cands = list(range(P))

            for b_val in initial_cands:
                is_valid = True
                for i in range(L_H):
                    if (idx == 8 and i == 1) or (idx == 9 and i == 0):
                        if is_commute([a_val, b_val], [a_vec[i], b_vec[i]], P):
                            is_valid = False
                            break
                    else:
                        if not is_commute([a_val, b_val], [a_vec[i], b_vec[i]], P):
                            is_valid = False
                            break
                            
                if is_valid:
                    b_cand.append(b_val)

        if not b_cand:
            continue
            
        random.shuffle(b_cand)

        # 見つかった候補の中からサイクルチェックを実行
        for b_val in b_cand: 
            new_a_vec[idx] = a_val
            new_b_vec[idx] = b_val
            if find_closed_cycle_x(new_cycles_x, new_a_vec, new_b_vec) and find_closed_cycle_z(new_cycles_z, new_a_vec, new_b_vec):
                return a_val, b_val
        # if a_count % 50 == 0:
        #     print(f"  (idx {idx}: {a_count}個の 'a' を試行中...)", end="\r")
    # 全てのaの候補を試しても見つからなければ、明確な手詰まりとしてNoneを返す
    return None
    