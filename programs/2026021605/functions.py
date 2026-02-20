import math
import itertools
import numpy as np
import random
from dataclasses import dataclass
import matplotlib.pyplot as plt
from collections import Counter
import sympy as sp
from sympy import gcd, gcdex, Matrix, list2numpy, eye
from itertools import product
from z3 import *

L=12
l_h = L // 2
J = 3
P=768

def generate_coprime_array():
    candidates = [i for i in range(2, P) if math.gcd(i, P) == 1]
    
    result = random.choices(candidates, k=L)
    return result

def generate_g(a_vec):
    G = np.zeros((36, 12), dtype=int)
    for i in range(6):
        for j in range(6):
            G[6*i+j, i] = (1 - a_vec[6+j])
            G[6*i+j, 6+j] = (a_vec[i]-1)
    return G

def print_g_matrix(matrix, split_idx=6):
    for row in matrix:
        left = row[:split_idx]
        right = row[split_idx:]
        left_str = " ".join(f"{val:3}" for val in left)
        right_str = " ".join(f"{val:3}" for val in right)
        print(f"{left_str} | {right_str}")

def print_g_smart(matrix, split_idx=6):
    for row in matrix:
        elements = []
        for i, val in enumerate(row):
            if i == split_idx:
                elements.append("|")
            if val == 0:
                elements.append(f"   ·")
            else:
                elements.append(f"{val:4}")
        
        # 行番号を付けて表示
        print(f"{' '.join(elements)}")

def solve_with_z3(cond_a, cond_b, cond_c):
    solver = Solver()
    b = [Int(f'b_{i}') for i in range(L)]
    for x in b:
        solver.add(x >= 0, x < P)

    for row in cond_a:
        expr = Sum([row[i] * b[i] for i in range(L)])
        solver.add(expr % P == 0)

    for row in cond_b:
        expr = Sum([row[j] * b[j] for j in range(L)])
        solver.add(expr % P != 0)
        
    for row, a_c in cond_c:
        expr = Sum([row[j] * b[j] for j in range(L)])
        solver.add(expr % a_c != 0)

    if solver.check() == sat:
        model = solver.model()
        return [model[x].as_long() for x in b]
    else:
        return None

def is_commute(f, g):
    a_1 = f[0]
    b_1 = f[1]
    a_2 = g[0]
    b_2 = g[1]
    if ((a_1-1)*b_2 - (a_2-1)*b_1)%P == 0:
        return 1
    else:
        return 0
    
def commute_matrix(a_vec, b_vec):
    result = []
    for i in range(l_h):
        row = []
        for j in range(l_h):
            row.append(is_commute([a_vec[i], b_vec[i]], [a_vec[6+j], b_vec[6+j]]))
        result.append(row)
    return result 

def is_valid(r_seq, c_seq):
    length = len(r_seq)
    for i in range(length):
        if r_seq[i] == r_seq[(i+1)%length]:
            return False
        if c_seq[i] == c_seq[(i+1)%length]:
            return False
    return True
    

def generate_cycles(max_len):
    def generate_utcbc():
        utcbcs = set()
        r_seq = [0, 1, 2, 1]
        for c0 in range(L):
            c2 = (c0 + 1) % l_h + (l_h if c0 >= l_h else 0)
            for c1 in range(L):
                c3 = (c1 + 1) % l_h + (l_h if c1 >= l_h else 0)
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
    cycles = cycles - generate_utcbc()
    return list(cycles)

def composite_affine(left, right):
    a_new = (left[0] * right[0]) % P
    b_new = (left[0] * right[1] + left[1]) % P
    return [a_new, b_new]

def func_inv(input):
        a,b = input
        try:
            a_inv = pow(a, -1, P)
        except ValueError:
            raise ValueError("Inverse does not exist")
        b_new = (-1 * a_inv * b) % P
        return [a_inv, b_new]

def generate_h_xz():
    h_x = []
    h_z = []
    for i in range(l_h):
        row_x = []
        row_z = []
        for j in range(l_h):
            idx = (j-i)%l_h
            row_x.append(idx)
            idx = (i-j)%l_h
            row_z.append(6+idx)
        for j in range(l_h):
            idx = (j-i)%l_h
            row_x.append(6+idx)
            idx = (i-j)%l_h
            row_z.append(idx)
        h_x.append(row_x)
        h_z.append(row_z)
    return h_x, h_z

def inv(val):
    try:
        return pow(val, -1, P)
    except ValueError:
        raise ValueError("Inverse does not exist")
    
   
def generate_constraints(cycles, a_vec, h_x, h_z):
    constraints = []
    for cycle in cycles:
        N = len(cycle)
        row_x = [0] * L
        row_z = [0] * L
        idx_x = [h_x[r][c] for r, c in cycle]
        idx_z = [h_z[r][c] for r, c in cycle]
        a_x = [a_vec[idx] for idx in idx_x]
        a_z = [a_vec[idx] for idx in idx_z]

        term_x = inv(a_x[0])
        for i in range(0, N, 2):
            row_x[idx_x[i]] = (row_x[idx_x[i]] - term_x) % P
            row_x[idx_x[i+1]] = (row_x[idx_x[i+1]] + term_x) % P
            if i + 2 < N:
                term_x = (term_x * a_x[i+1] * inv(a_x[i+2])) % P
        a_c_x = 1
        for i in range(0, N, 2):
            a_c_x = (a_c_x * inv(a_x[i]) * a_x[i+1]) % P
        gcd_x = math.gcd(a_c_x-1, P)
        row_x = [val * (P//gcd_x) % P for val in row_x]
        constraints.append(row_x)

        row_z[idx_z[0]] = (row_z[idx_z[0]] + 1) % P
        term_z = 1
        for i in range(0, N - 2, 2):
            term_z = (term_z * a_z[i] * inv(a_z[i+1])) % P
            row_z[idx_z[i+1]] = (row_z[idx_z[i+1]] - term_z) % P
            row_z[idx_z[i+2]] = (row_z[idx_z[i+2]] + term_z) % P
        term_z = (term_z * a_z[N-2] * inv(a_z[N-1])) % P
        row_z[idx_z[N-1]] = (row_z[idx_z[N-1]] - term_z) % P
        a_c_z = 1
        for i in range(0, N, 2):
            a_c_z = (a_c_z * a_z[i] * inv(a_z[i+1])) % P
        gcd_z = math.gcd(a_c_z-1, P)
        row_z = [val * (P//gcd_z) % P for val in row_z]
        constraints.append(row_z)
    return constraints

def is_equiv(func1, func2):
    if (func1 != func2):
        print(False)
    return 0

def generate_functions(cycles, a_vec, b_vec, h_x, h_z):
    functions = []
    for cycle in cycles:
        idx_x = [h_x[r][c] for r, c in cycle]
        idx_z = [h_z[r][c] for r, c in cycle]
        a_x = [a_vec[idx] for idx in idx_x]
        b_x = [b_vec[idx] for idx in idx_x]
        a_z = [a_vec[idx] for idx in idx_z]
        b_z = [b_vec[idx] for idx in idx_z]
        function_x = [1, 0]
        function_z = [1, 0]
        for i in range(len(cycle) // 2):
            function_x = composite_affine(function_x, [a_x[2*i], b_x[2*i]])
            function_x = composite_affine(function_x, func_inv([a_x[2*i+1], b_x[2*i+1]]))
            function_z = composite_affine(function_z, func_inv([a_z[2*i], b_z[2*i]]))
            function_z = composite_affine(function_z, [a_z[2*i+1], b_z[2*i+1]])
            
        functions.append(function_x)
        functions.append(function_z)
    return functions

def is_closed(input):
    a = input[0]
    b = input[1]
    if a == 1 and b == 0:
        return True
    d = math.gcd(a-1, P)
    if b % d == 0:
        return True
    else:
        return False
    
def make_matrix(func):
    a,b = func
    matrix = []
    for i in range(P):
        row = []
        for j in range(P):
            if (a*i + b) % P == j:
                row.append(1)
            else:
                row.append(0)
        matrix.append(row)
    return matrix
def create_h_x(a_vec, b_vec, h_x):
    rows = []
    for i in range(l_h):
        current_row_matrices = []
        for j in range(L):
            idx = h_x[i][j]
            params = [a_vec[idx], b_vec[idx]]
            mat = make_matrix(params)
            current_row_matrices.append(mat)
        rows.append(current_row_matrices)
    h_x = np.block(rows)
    return h_x
def create_h_z(a_vec, b_vec, h_z):
    rows = []
    for i in range(l_h):
        current_row_matrices = []
        for j in range(L):
            idx = h_z[i][j]
            params = [a_vec[idx], b_vec[idx]]
            mat = make_matrix(params)
            current_row_matrices.append(mat.T)
        rows.append(current_row_matrices)
    h_x = np.block(rows)
    return h_x

def generate_g(a_vec):
    G = Matrix.zeros(l_h**2, L)
    for i in range(l_h):
        for j in range(l_h):
            G[l_h*i+j, i] = (1 - a_vec[l_h+j])
            G[l_h*i+j, l_h+j] = (a_vec[i]-1)
    return G

def check(a_vec, b_vec):
    commute_mat = commute_matrix(a_vec, b_vec)
    if commute_mat != [[1, 1, 1, 0, 1, 1],
    [1, 1, 0, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1]]:
        return False
    cycles = generate_cycles(6)
    h_x, h_z = generate_h_xz()
    functions = generate_functions(cycles, a_vec, b_vec, h_x, h_z)
    for function in functions:
        if is_closed(function):
            return False
    return True

def crt_combine_val(v3, v256):
    """
    v = v3 (mod 3) かつ v = v256 (mod 256) となる v (mod 768) を求める
    """
    m1, m2 = 3, 256
    # Garner's Algorithm の簡略版: x = v3 + m1 * ((v256 - v3) * inv(m1, m2) % m2)
    m1_inv = pow(m1, -1, m2)
    h = ((v256 - v3) * m1_inv) % m2
    return (v3 + m1 * h) % 768

def crt_combine_matrix(M3, M256):
    """
    Matrix の各要素に対して CRT 合成を行う
    """
    res = M3.copy()
    for i in range(res.rows):
        for j in range(res.cols):
            res[i, j] = crt_combine_val(M3[i, j], M256[i, j])
    return res

def crt_combine_val(v3, v256):
    """
    v = v3 (mod 3) かつ v = v256 (mod 256) となる v (mod 768) を求める
    """
    m1, m2 = 3, 256
    # Garner's Algorithm の簡略版: x = v3 + m1 * ((v256 - v3) * inv(m1, m2) % m2)
    m1_inv = pow(m1, -1, m2)
    h = ((v256 - v3) * m1_inv) % m2
    return (v3 + m1 * h) % 768

def crt_combine_matrix(M3, M256):
    """
    Matrix の各要素に対して CRT 合成を行う
    """
    res = M3.copy()
    for i in range(res.rows):
        for j in range(res.cols):
            res[i, j] = crt_combine_val(M3[i, j], M256[i, j])
    return res

def crt_lift(x_mod3, x_mod256):
    """
    x = x_mod3 (mod 3) かつ x = x_mod256 (mod 256) となる x (mod 768) を合成する。
    中国剰余定理の公式: x = a1*M1*y1 + a2*M2*y2 (mod M)
    """
    m1, m2 = 3, 256
    M = m1 * m2 # 768
    
    # M1 = 256, M2 = 3
    # y1 = inv(256, 3) = 1, y2 = inv(3, 256) = 171
    y1 = pow(m2, -1, m1)
    y2 = pow(m1, -1, m2)
    
    res = (x_mod3 * m2 * y1 + x_mod256 * m1 * y2) % M
    return res

# 行列/ベクトル全体に適用するラッパー
def crt_combine_vectors(v3, v256):
    res = v3.copy()
    for i in range(res.rows):
        res[i, 0] = crt_lift(v3[i, 0], v256[i, 0])
    return res

def solve_modular_kernel(A, P):
    """
    剰余環 Z_P 上における行列 A の核 (Ax = 0 mod P) の基底を求める。
    
    引数:
        A (sympy.Matrix): 対象となる行列
        P (int): 法 (Modulus)
        
    戻り値:
        list: 法 P における解空間の基底ベクトルのリスト
    """
    rows, cols = A.shape
    
    # 整数上の方程式 Ax + Pk = 0 を解くために、[A | P*I] という行列を作る。
    # この行列の整数核 (nullspace) を求めれば、その最初の 'cols' 成分が Ax = 0 mod P の解になる。
    M = Matrix.hstack(A, P * eye(rows))
    
    # 整数体上での核を計算
    null_basis = M.nullspace()
    
    basis_vectors = []
    for v in null_basis:
        # 最初の cols 行が変数 x に対応する部分
        u = v[:cols, :]
        
        # SymPyの計算過程で有理数（分数）が含まれる場合があるため、
        # 分母を払って整数に変換してから mod P を適用する
        denoms = [val.q for val in u if val != 0]
        lcm = 1
        for d in denoms:
            import math
            lcm = (lcm * d) // math.gcd(lcm, d)
        
        u_int = (u * lcm).applyfunc(lambda x: x % P)
        
        # 零ベクトルでなければ基底として採用
        if not u_int.is_zero_matrix:
            # すでに保存されている基底と重複（スカラー倍など）していないか簡易チェック
            is_new = True
            for existing in basis_vectors:
                if u_int == existing:
                    is_new = False
                    break
            if is_new:
                basis_vectors.append(u_int)
                
    return basis_vectors

from z3 import *
from sympy import Matrix

def solve_decomposition_z3(V, b, p):
    """
    V * c = b (mod p) となる係数 c を Z3 (整数計画) で求める。
    V: カーネル基底行列
    b: ターゲットベクトル
    p: 法
    """
    solver = Solver()
    k = V.cols # 係数 c の数（基底の数）
    rows = V.rows
    
    # 変数 c_0, ..., c_{k-1} の定義
    c = [Int(f'c_{i}') for i in range(k)]
    
    # 範囲制約 (任意ですが、0以上P未満の解を探すと綺麗です)
    for var in c:
        solver.add(var >= 0, var < p)
        
    # 方程式の制約: V * c ≡ b (mod p)
    # 各行について: sum(V_ij * c_j) - b_i = m_i * p
    # Z3の % 演算子を使って記述します
    for i in range(rows):
        # 行列の成分は sympy の Integer なので int() でキャスト
        expr = Sum([int(V[i, j]) * c[j] for j in range(k)])
        target = int(b[i, 0])
        solver.add((expr - target) % p == 0)
        
    # 解の探索
    if solver.check() == sat:
        model = solver.model()
        # 結果を SymPy Matrix として返す
        return Matrix([model[c[i]].as_long() for i in range(k)])
    else:
        return None
    
    
def is_in_general_solution_strict(x_vec, V, forbidden_vectors, a_vec, p):
    if x_vec is None: return False
    
    # 1. b_vec の復元
    b_mat = (V * x_vec).applyfunc(lambda x: x % p)
    b_list = [int(val) for val in b_mat]
    
    # 2. 条件A (特定の可換性) のチェック
    # インデックス (1, 3) が 0 になっているものは、条件Aを満たさないので排除
    # commute_matrix を使わずに直接計算することで高速化
    # f1 = [a_vec[1], b_list[1]], g3 = [a_vec[6+3], b_list[6+3]]
    # (a_1 - 1)*b_9 - (a_9 - 1)*b_1 == 0 (mod P)
    a1, b1 = a_vec[1], b_list[1]
    a9, b9 = a_vec[9], b_list[9] # index 6+3 = 9
    if ((a1 - 1) * b9 - (a9 - 1) * b1) % p != 0:
        return False # 可換でなくなった場合は不合格

    # 3. 条件B/C (禁止ベクトル/ガース) のチェック
    x_mat = Matrix(x_vec)
    for r in forbidden_vectors:
        if (r.T * x_mat)[0] % p == 0:
            return False # サイクルが閉じる(または非可換部が可換になる)場合は不合格
            
    return True


def reconstruct_b_vectors(x_solutions, V_matrix, p):
    """
    da 次元の解 x から L 次元の係数ベクトル b を復元する
    b = V * x (mod p)
    """
    b_vectors = []
    for x_sol in x_solutions:
        # 基底行列 V との積を計算し mod P を適用
        b_mat = (V_matrix * x_sol).applyfunc(lambda val: val % p)
        # Matrix 型から Python のリスト形式に変換
        b_list = [int(val) for val in b_mat]
        b_vectors.append(b_list)
    return b_vectors