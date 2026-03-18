# matrix_builder.py
import math
import numpy as np
from config import P, L, l_h

def gen_h_xz():
    """LDPCのベースとなるパリティ検査行列の配置を生成する"""
    h_x = []
    h_z = []
    for i in range(l_h):
        row_x = []
        row_z = []
        for j in range(l_h):
            idx = (j - i) % l_h
            row_x.append(idx)
            idx = (i - j) % l_h
            row_z.append(6 + idx)
        for j in range(l_h):
            idx = (j - i) % l_h
            row_x.append(6 + idx)
            idx = (i - j) % l_h
            row_z.append(idx)
        h_x.append(row_x)
        h_z.append(row_z)
    return h_x, h_z

def gen_g(a_vec):
    """aのベクトルから、可換性を表現するためのベース行列Gを生成する"""
    G = np.zeros((36, 12), dtype=int)
    for i in range(6):
        for j in range(6):
            G[6*i+j, i] = (1 - a_vec[6+j])
            G[6*i+j, 6+j] = (a_vec[i] - 1)
    return G

def gen_g_ab(G):
    """行列Gを、可換条件(Ga)と非可換条件(Gb)に分割する"""
    all_indices = [i for i in range(l_h**2)]
    gb_indices = [3, 8]
    ga_indices = [i for i in all_indices if i not in gb_indices]
    Ga = G[ga_indices]
    Gb = G[gb_indices]
    return Ga, Gb

def gen_g_mat(a_vec):
    """a_vecから直接GaとGbを生成するヘルパー関数"""
    G = gen_g(a_vec)
    return gen_g_ab(G)

def gen_c_constraints(cycles, a_vec, h_x, h_z):
    """ガース制約(条件C)を行列形式で生成する"""
    def inv(val):
        try:
            return pow(val, -1, P)
        except ValueError:
            raise ValueError("Inverse does not exist")
    
    constraints = []
    for cycle in cycles:
        N = len(cycle)
        row_x = [0] * L
        row_z = [0] * L
        idx_x = [h_x[r][c] for r, c in cycle]
        idx_z = [h_z[r][c] for r, c in cycle]
        a_x = [a_vec[idx] for idx in idx_x]
        a_z = [a_vec[idx] for idx in idx_z]

        # X側の制約構築
        term_x = inv(a_x[0]) * a_x[1]
        row_x[idx_x[0]] -= 1
        row_x = [(val * term_x) % P for val in row_x]
        
        for i in range(1, N - 1, 2):
            term_x = (inv(a_x[i+1]) * a_x[i+2]) % P
            row_x[idx_x[i]] += 1
            row_x[idx_x[i+1]] -= 1
            row_x = [(val * term_x) % P for val in row_x]
            
        row_x[idx_x[N-1]] += 1
        
        a_c_x = 1
        for i in range(0, N, 2):
            a_c_x = (a_c_x * inv(a_x[i]) * a_x[i+1]) % P
        mul = P // math.gcd(a_c_x - 1, P)
        row_x = [(val * mul) % P for val in row_x]
        constraints.append(row_x)

        # Z側の制約構築
        for i in range(0, N - 2, 2):
            term_z = (inv(a_z[i+1]) * a_z[i+2]) % P
            row_z[idx_z[i]] += 1
            row_z[idx_z[i+1]] -= 1
            row_z = [(val * term_z) % P for val in row_z]

        term_z = inv(a_z[N-1]) % P
        row_z[idx_z[N-2]] += 1
        row_z[idx_z[N-1]] -= 1
        row_z = [(val * term_z) % P for val in row_z]
        
        a_c_z = 1
        for i in range(0, N, 2):
            a_c_z = (a_c_z * a_z[i] * inv(a_z[i+1])) % P
        mul = P // math.gcd(a_c_z - 1, P)
        row_z = [(val * mul) % P for val in row_z]
        constraints.append(row_z)
        
    return constraints

def gen_constraints(Gb, a_vec, cycles, h_x, h_z):
    """条件Bと条件Cの「禁止ベクトル(0になってはいけないベクトル)」を統合する"""
    constraints = gen_c_constraints(cycles, a_vec, h_x, h_z)

    unique_forbidden_vectors = []
    seen_vectors = set()

    # 1. 条件B (潜在部の非可換性) からの制約
    for row in Gb:
        row_tuple = tuple(row)
        if row_tuple not in seen_vectors:
            unique_forbidden_vectors.append(row)
            seen_vectors.add(row_tuple)

    # 2. 条件C (短いサイクルの回避) からの制約
    for row in constraints:
        row_tuple = tuple(row)
        if row_tuple not in seen_vectors:
            unique_forbidden_vectors.append(row)
            seen_vectors.add(row_tuple)
            
    return unique_forbidden_vectors