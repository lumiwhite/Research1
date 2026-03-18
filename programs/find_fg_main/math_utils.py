# math_utils.py
import math
import random
import numpy as np
import sympy as sp
from config import P, L

def random_a(p_val=P):
    """互いに素な a をランダムに生成する"""
    while True:
        a_val = random.randint(0, p_val - 1)
        if math.gcd(a_val, p_val) == 1:
            return a_val

def gen_coprime_array():
    """解空間を広げるための「エリート候補」を含めた a の配列を生成する"""
    normal_candidates = [i for i in range(2, P) if math.gcd(i, P) == 1]
    
    # a - 1 が P と大きな公約数を持つ「エリート候補」を抽出
    elite_candidates = []
    for a in normal_candidates:
        if math.gcd(a - 1, P) >= 64:
            elite_candidates.append(a)
            
    result = []
    for _ in range(L):
        if random.random() < 0.8 and elite_candidates:
            result.append(random.choice(elite_candidates))
        else:
            result.append(random.choice(normal_candidates))
            
    return result

def composite_affine(left, right):
    """2つのアフィン関数の合成を行う"""
    a_new = (left[0] * right[0]) % P
    b_new = (left[0] * right[1] + left[1]) % P
    return [a_new, b_new]

def func_inv(input_func):
    """アフィン関数の逆関数を求める"""
    a, b = input_func
    try:
        a_inv = pow(a, -1, P)
    except ValueError:
        raise ValueError("Inverse does not exist")
    b_new = (-1 * a_inv * b) % P
    return [a_inv, b_new]

def solve_snf(Ga_list, p_val=P):
    """
    スミス標準形を用いて、Ga * b ≡ 0 (mod P) を満たす解空間を特定する。
    """
    if len(Ga_list) == 0 or len(Ga_list[0]) == 0:
        return [], []

    n = len(Ga_list)
    m = len(Ga_list[0])
    
    if isinstance(Ga_list, np.ndarray):
        A = Ga_list.tolist()
    else:
        A = [list(row) for row in Ga_list]
        
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    V = [[1 if i == j else 0 for j in range(m)] for i in range(m)]

    def row_add(i, j, q):
        for k in range(m): A[i][k] -= q * A[j][k]
        for k in range(n): U[i][k] -= q * U[j][k]
    def col_add(i, j, q):
        for k in range(n): A[k][i] -= q * A[k][j]
        for k in range(m): V[k][i] -= q * V[k][j]
    def row_swap(i, j):
        A[i], A[j] = A[j], A[i]
        U[i], U[j] = U[j], U[i]
    def col_swap(i, j):
        for k in range(n): A[k][i], A[k][j] = A[k][j], A[k][i]
        for k in range(m): V[k][i], V[k][j] = V[k][j], V[k][i]
    def col_mult(i, q):
        for k in range(n): A[k][i] *= q
        for k in range(m): V[k][i] *= q

    for t in range(min(n, m)):
        while True:
            min_val = float('inf')
            pi, pj = -1, -1
            for i in range(t, n):
                for j in range(t, m):
                    if A[i][j] != 0 and abs(A[i][j]) < min_val:
                        min_val = abs(A[i][j])
                        pi, pj = i, j
            
            if pi == -1:
                break
                
            if pi != t: row_swap(t, pi)
            if pj != t: col_swap(t, pj)
            if A[t][t] < 0: col_mult(t, -1)
            
            changed = False
            for j in range(t + 1, m):
                if A[t][j] != 0:
                    q = A[t][j] // A[t][t]
                    col_add(j, t, q)
                    changed = True
            for i in range(t + 1, n):
                if A[i][t] != 0:
                    q = A[i][t] // A[t][t]
                    row_add(i, t, q)
                    changed = True
            
            if not changed:
                divides_all = True
                for i in range(t + 1, n):
                    for j in range(t + 1, m):
                        if A[i][j] % A[t][t] != 0:
                            row_add(t, i, -1)
                            changed = True
                            divides_all = False
                            break
                    if not divides_all: break
                if not changed:
                    break

    for i in range(min(n, m)):
        if A[i][i] < 0:
            col_mult(i, -1)

    y_ranges = []
    for i in range(m):
        d_i = A[i][i] if i < min(n, m) else 0
        if d_i == 0:
            step = 1
            num_vals = p_val
        else:
            g = math.gcd(d_i, p_val)
            step = p_val // g
            num_vals = g
        y_ranges.append((step, num_vals))
        
    return y_ranges, V

def solve_hnf(Ga_list, p_val=P):
    """
    エルミート標準形のアプローチを用いて、Ga * b ≡ 0 (mod P) を満たす解の基底を求める。
    """
    if len(Ga_list) == 0 or len(Ga_list[0]) == 0:
        return []
        
    rows = len(Ga_list)
    cols = len(Ga_list[0])
    
    Ga = sp.Matrix(Ga_list)
    P_I = p_val * sp.eye(rows)
    M = Ga.row_join(P_I)
    
    null_basis_rat = M.nullspace()
    
    basis_vectors = []
    for vec in null_basis_rat:
        lcm_val = 1
        for val in vec:
            lcm_val = sp.lcm(lcm_val, val.q)
        
        int_vec = vec * lcm_val
        b_vec = [int(int_vec[i]) % p_val for i in range(cols)]
        
        if any(v != 0 for v in b_vec):
            basis_vectors.append(b_vec)
            
    unique_basis = []
    for b in basis_vectors:
        if b not in unique_basis:
            unique_basis.append(b)
            
    return unique_basis