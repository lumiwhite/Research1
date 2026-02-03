import math
import itertools
import numpy as np
import random
from dataclasses import dataclass
import matplotlib.pyplot as plt
from collections import Counter
import sympy as sp
from sympy import gcd, gcdex
from itertools import product
from z3 import *

L=12
l_h = L // 2
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
    """
    P: 法
    L: ベクトルの次元 (b_0, ..., b_{L-1})
    M_comm: 可換性条件の係数行列 (リストのリスト)
    forbidden_constraints: (係数ベクトル w, 法 g) のリスト。 w・b % g != 0 を満たす必要がある。
    """
    solver = Solver()
    
    # 変数の定義 (0 <= b_i < P)
    b = [Int(f'b_{i}') for i in range(L)]
    for x in b:
        solver.add(x >= 0, x < P)

    # 1. 等式制約: M_comm * b == 0 (mod P)
    for row in cond_a:
        # sum(row[i] * b[i]) % P == 0
        expr = Sum([row[i] * b[i] for i in range(L)])
        solver.add(expr % P == 0)

    for row in cond_b:
        expr = Sum([row[j] * b[j] for j in range(L)])
        solver.add(expr % P != 0)
        
    for row, a_c in cond_c:
        expr = Sum([row[j] * b[j] for j in range(L)])
        solver.add(expr % a_c != 0)

    # 解の探索
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

def generate_cycles():
    def generate_utcbc():
        utcbcs = set()
        for row_0 in range(l_h):
            for row_1 in range(l_h):
                for row_2 in range(l_h):
                        for col_0 in range(L):
                             for col_1 in range(L):
                                rows = [row_0, row_1, row_2]
                                cols = [col_0, col_1]
                                if  cols[0] < l_h:
                                    col_2 = (rows[1]+cols[0]-rows[0]) % l_h
                                else:
                                    col_2 = l_h + (rows[1]+cols[0]-rows[0]) % l_h
                                if  cols[1] < l_h:
                                    col_3 = (rows[2]+cols[1]-rows[1]) % l_h
                                else:
                                    col_3 = l_h + (rows[2]+cols[1]-rows[1]) % l_h
                                row_3 = (col_3+rows[0]-cols[1]) % l_h
                                rows.extend([row_3])
                                cols.extend([col_2, col_3])
                                is_valid = True
                                for i in range(4):
                                    if rows[i] == rows[(i+1)%4] or cols[i] == cols[(i+1)%4]:
                                        is_valid = False
                                        break
                                
                                if is_valid:
                                    utcbcs.add(canonicalize(rows, cols))
        return utcbcs

    def get_valid_sequences(elements, length=4):
        valid = []
        for seq in product(elements, repeat=length):
            # 隣接チェック (r1!=r2, r2!=r3, r3!=r4, r4!=r1)
            if all(seq[i] != seq[(i+1)%length] for i in range(length)):
                valid.append(seq)
        return valid

    def canonicalize(r_seq, c_seq):
        positions = tuple((r_seq[i], c_seq[i]) for i in range(4))
        symmetries = []
        curr = list(positions)
        curr_equiv = get_equiv_pos(curr)
        for _ in range(4):
            curr = curr[1:] + curr[:1]
            symmetries.append(tuple(curr))
            curr_equiv = curr_equiv[1:] + curr_equiv[:1]
            symmetries.append(tuple(curr_equiv))
        return min(symmetries)

    def get_equiv_pos(pos):
        pos1 = tuple((pos[0][0], pos[1][1]))
        pos2 = tuple((pos[3][0], pos[0][1]))
        pos3 = tuple((pos[2][0], pos[3][1]))
        pos4 = tuple((pos[1][0], pos[2][1]))
        return [pos1, pos2, pos3, pos4]
    c_4 = []
    for r1 in range(l_h):
        for r2 in range(r1 + 1, l_h): # r1より大きいインデックスのみ選ぶ
            for c1 in range(L):
                for c2 in range(c1 + 1, L): # c1より大きいインデックスのみ選ぶ
                    pos1 = [r1, c1]
                    pos2 = [r2, c2]
                    c_4.append([pos1, pos2])
    c_6 = []
    for r_set in itertools.combinations(list(range(l_h)), 3):
        r1, r2, r3 = r_set
        for c_set in itertools.combinations(list(range(L)), 3):
            for c_perm in itertools.permutations(c_set):
                c1, c2, c3 = c_perm
                cycle_positions = [(r1, c1), (r2, c2), (r3, c3)]
                c_6.append(cycle_positions)

    row_seqs = get_valid_sequences(range(l_h))  # 630通り
    col_seqs = get_valid_sequences(range(L))  # 14652通り
    c_8 = []
    c_8 = set()
    for r in row_seqs:
        for c in col_seqs:
            c_8.add(canonicalize(r, c))
    c_8 = list(c_8 - generate_utcbc())
    return c_4, c_6, c_8

def get_cycle(pos):
    length = len(pos)
    cycle = []
    for i in range(length):
        cycle.append([pos[i][0], pos[i][1]])
        cycle.append([pos[i][0], pos[(i+1)%length][1]])
    return cycle

def composite_affine(left, right):
    a_new = (left[0] * right[0]) % P
    b_new = (left[0] * right[1] + left[1]) % P
    return [a_new, b_new]

def func_inv(input):
        try:
            a_inv = pow(input[0], -1, P)
        except ValueError:
             raise ValueError("Inverse does not exist")
        b_new = (-1 * a_inv * input[1]) % P
        return [a_inv, b_new]

def generate_h_xz():
    h_x = []
    for i in range(l_h):
        row = []
        for j in range(l_h):
            idx = (j-i)%l_h
            row.append(idx)
        for j in range(l_h):
            idx = (j-i)%l_h
            row.append(6+idx)
        h_x.append(row)
    h_z = []
    for i in range(l_h):
        row = []
        for j in range(l_h):
            idx = (i-j)%l_h
            row.append(6+idx)
        for j in range(l_h):
            idx = (i-j)%l_h
            row.append(idx)
        h_x
        h_z.append(row)
    return h_x, h_z

def inv(val):
    try:
        return pow(val, -1, P)
    except ValueError:
        raise ValueError("Inverse does not exist")
    
def generate_constraints(c_4, c_6, c_8, a_vec, h_x, h_z):
    constraints = []
    cycle_4 = [get_cycle(c) for c in c_4]
    for cycle in cycle_4:
        row_x = [0]*L
        row_z = [0]*L
        functions_x = []
        functions_z = []
        for i in range(4):
            functions_x.append(h_x[cycle[i][0]][cycle[i][1]])
            functions_z.append(h_z[cycle[i][0]][cycle[i][1]])
        idx_x = [functions_x[i] for i in range(4)]
        idx_z = [functions_z[i] for i in range(4)]
        a_x = [a_vec[idx] for idx in idx_x]
        a_z = [a_vec[idx] for idx in idx_z]
        row_x[idx_x[0]] -= inv(a_x[0])
        row_x[idx_x[1]] += inv(a_x[0])
        row_x[idx_x[2]] -= inv(a_x[0]) * a_x[1] * inv(a_x[2])
        row_x[idx_x[3]] += inv(a_x[0]) * a_x[1] * inv(a_x[2])
        a_c_x = inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] % P
        row_x = [val % P for val in row_x]
        constraints.append([row_x, a_c_x])
        row_z[idx_z[0]] += 1
        row_z[idx_z[1]] -= a_z[0] * inv(a_z[1])
        row_z[idx_z[2]] += a_z[0] * inv(a_z[1])
        row_z[idx_z[3]] -= a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3])
        a_c_z = a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3]) % P
        row_z = [val % P for val in row_z]
        constraints.append([row_z, a_c_z]) 

    cycle_6 = [get_cycle(c) for c in c_6]
    for cycle in cycle_6:
        row_x = [0]*L
        row_z = [0]*L
        functions_x = []
        functions_z = []
        for i in range(6):
            functions_x.append(h_x[cycle[i][0]][cycle[i][1]])
            functions_z.append(h_z[cycle[i][0]][cycle[i][1]])
        idx_x = [functions_x[i] for i in range(6)]
        idx_z = [functions_z[i] for i in range(6)]
        a_x = [a_vec[idx] for idx in idx_x]
        a_z = [a_vec[idx] for idx in idx_z]
        row_x[idx_x[0]] -= inv(a_x[0])
        row_x[idx_x[1]] += inv(a_x[0])
        row_x[idx_x[2]] -= inv(a_x[0]) * a_x[1] * inv(a_x[2])
        row_x[idx_x[3]] += inv(a_x[0]) * a_x[1] * inv(a_x[2])
        row_x[idx_x[4]] -= inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4])
        row_x[idx_x[5]] += inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4])
        a_c_x = inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4]) * a_x[5] % P
        row_x = [val % P for val in row_x]
        constraints.append([row_x, a_c_x])
        row_z[idx_z[0]] += 1
        row_z[idx_z[1]] -= a_z[0] * inv(a_z[1])
        row_z[idx_z[2]] += a_z[0] * inv(a_z[1])
        row_z[idx_z[3]] -= a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3])
        row_z[idx_z[4]] += a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3])
        row_z[idx_z[5]] -= a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3]) * a_z[4] * inv(a_z[5])
        a_c_z = a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3]) * a_z[4] * inv(a_z[5]) % P
        row_z = [val % P for val in row_z]
        constraints.append([row_z, a_c_z])

    cycle_8 = [get_cycle(c) for c in c_8]
    for cycle in cycle_8:
        row_x = [0]*L
        row_z = [0]*L
        functions_x = []
        functions_z = []
        for i in range(8):
            functions_x.append(h_x[cycle[i][0]][cycle[i][1]])
            functions_z.append(h_z[cycle[i][0]][cycle[i][1]])
        idx_x = [functions_x[i] for i in range(8)]
        idx_z = [functions_z[i] for i in range(8)]
        a_x = [a_vec[idx] for idx in idx_x]
        a_z = [a_vec[idx] for idx in idx_z]
        row_x[idx_x[0]] -= inv(a_x[0])
        row_x[idx_x[1]] += inv(a_x[0])
        row_x[idx_x[2]] -= inv(a_x[0]) * a_x[1] * inv(a_x[2])
        row_x[idx_x[3]] += inv(a_x[0]) * a_x[1] * inv(a_x[2])
        row_x[idx_x[4]] -= inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4])
        row_x[idx_x[5]] += inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4])
        row_x[idx_x[4]] -= inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4]) * a_x[5] * inv(a_x[6])
        row_x[idx_x[5]] += inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4]) * a_x[5] * inv(a_x[6])
        a_c_x = inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4]) * a_x[5] % P
        row_x = [val % P for val in row_x]
        constraints.append([row_x, a_c_x])
        row_z[idx_z[0]] += 1
        row_z[idx_z[1]] -= a_z[0] * inv(a_z[1])
        row_z[idx_z[2]] += a_z[0] * inv(a_z[1])
        row_z[idx_z[3]] -= a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3])
        row_z[idx_z[4]] += a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3])
        row_z[idx_z[5]] -= a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3]) * a_z[4] * inv(a_z[5])
        row_z[idx_z[6]] += a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3]) * a_z[4] * inv(a_z[5])
        row_z[idx_z[7]] -= a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3]) * a_z[4] * inv(a_z[5]) * a_z[6] * inv(a_z[7])
        a_c_z = a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3]) * a_z[4] * inv(a_z[5]) % P
        row_z = [val % P for val in row_z]
        constraints.append([row_z, a_c_z])
    return constraints

def is_equiv(func1, func2):
    if (func1 != func2):
        print(False)
    return 0

def generate_functions(c_4, c_6, c_8, a_vec, b_vec, h_x, h_z):
    functions = [] 
    cycle_4 = [get_cycle(c) for c in c_4]
    for cycle in cycle_4:
        functions_x = []
        functions_z = []
        for i in range(4):
            functions_x.append(h_x[cycle[i][0]][cycle[i][1]])
            functions_z.append(h_z[cycle[i][0]][cycle[i][1]])
        idx_x = [functions_x[i] for i in range(4)]
        idx_z = [functions_z[i] for i in range(4)]
        a_x = [a_vec[idx] for idx in idx_x]
        b_x = [b_vec[idx] for idx in idx_x]
        a_z = [a_vec[idx] for idx in idx_z]
        b_z = [b_vec[idx] for idx in idx_z]
        function_x = [1, 0]
        function_z = [1, 0]
        for i in range(2):
            function_x = composite_affine(function_x, func_inv([a_x[2*i], b_x[2*i]]))
            function_x = composite_affine(function_x, [a_x[2*i+1], b_x[2*i+1]])
            function_z = composite_affine(function_z, [a_z[2*i], b_z[2*i]])
            function_z = composite_affine(function_z, func_inv([a_z[2*i+1], b_z[2*i+1]]))
        # a_c_x = inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] % P
        # b_c_x = (inv(a_x[0]) * a_x[1] * inv(a_x[2])*(b_x[3]-b_x[2]) + inv(a_x[0])*(b_x[1]-b_x[0])) % P
        # a_c_z = a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3]) % P
        # b_c_z = (b_z[0] + a_z[0]*inv(a_z[1])*(b_z[2]-b_z[1])-a_z[0]*inv(a_z[1])*a_z[2]*inv(a_z[3])*b_z[3]) % P
        # is_equiv([a_c_x, b_c_x], function_x)
        # is_equiv([a_c_z, b_c_z], function_z)
        functions.append(function_x)
        functions.append(function_z)

    cycle_6 = [get_cycle(c) for c in c_6]
    for cycle in cycle_6:
        functions_x = []
        functions_z = []
        for i in range(6):
            functions_x.append(h_x[cycle[i][0]][cycle[i][1]])
            functions_z.append(h_z[cycle[i][0]][cycle[i][1]])
        idx_x = [functions_x[i] for i in range(6)]
        idx_z = [functions_z[i] for i in range(6)]
        a_x = [a_vec[idx] for idx in idx_x]
        b_x = [b_vec[idx] for idx in idx_x]
        a_z = [a_vec[idx] for idx in idx_z]
        b_z = [b_vec[idx] for idx in idx_z]
        function_x = [1, 0]
        function_z = [1, 0]
        for i in range(3):
            function_x = composite_affine(function_x, func_inv([a_x[2*i], b_x[2*i]]))
            function_x = composite_affine(function_x, [a_x[2*i+1], b_x[2*i+1]])
            function_z = composite_affine(function_z, [a_z[2*i], b_z[2*i]])
            function_z = composite_affine(function_z, func_inv([a_z[2*i+1], b_z[2*i+1]]))
        # a_c_x = inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4]) * a_x[5] % P
        # b_c_x = (inv(a_x[0])*a_x[1]*inv(a_x[2])*a_x[3]*inv(a_x[4])*(b_x[5]-b_x[4]) + inv(a_x[0])*a_x[1]*inv(a_x[2])*(b_x[3]-b_x[2]) + inv(a_x[0])*(b_x[1]-b_x[0])) % P
        # a_c_z = a_z[0]*inv(a_z[1])*a_z[2]*inv(a_z[3])*a_z[4]*inv(a_z[5]) % P
        # b_c_z = (b_z[0]+a_z[0]*inv(a_z[1])*(b_z[2]-b_z[1])+a_z[0]*inv(a_z[1])*a_z[2]*inv(a_z[3])*(b_z[4]-b_z[3])-a_z[0]*inv(a_z[1])*a_z[2]*inv(a_z[3])*a_z[4]*inv(a_z[5])*b_z[5]) % P
        # is_equiv([a_c_x, b_c_x], function_x)
        # is_equiv([a_c_z, b_c_z], function_z)
        functions.append(function_x)
        functions.append(function_z)

    # cycle_8 = [get_cycle(c) for c in c_8]
    # for cycle in cycle_8:
    #     functions_x = []
    #     functions_z = []
    #     for i in range(8):
    #         functions_x.append(h_x[cycle[i][0]][cycle[i][1]])
    #         functions_z.append(h_z[cycle[i][0]][cycle[i][1]])
    #     idx_x = [functions_x[i] for i in range(8)]
    #     idx_z = [functions_z[i] for i in range(8)]
    #     a_x = [a_vec[idx] for idx in idx_x]
    #     b_x = [b_vec[idx] for idx in idx_x]
    #     a_z = [a_vec[idx] for idx in idx_z]
    #     b_z = [b_vec[idx] for idx in idx_z]
    #     function_x = [1, 0]
    #     function_z = [1, 0]
    #     function_x = [1, 0]
    #     function_z = [1, 0]
    #     for i in range(4):
    #         function_x = composite_affine(function_x, func_inv([a_x[2*i], b_x[2*i]]))
    #         function_x = composite_affine(function_x, [a_x[2*i+1], b_x[2*i+1]])
    #         function_z = composite_affine(function_z, [a_z[2*i], b_z[2*i]])
    #         function_z = composite_affine(function_z, func_inv([a_z[2*i+1], b_z[2*i+1]]))
    #     # a_c_x = inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4]) * a_x[5] * inv(a_x[6]) * a_x[7] % P
    #     # b_c_x = (inv(a_x[0])*a_x[1]*inv(a_x[2])*a_x[3]*inv(a_x[4])*a_x[5]*inv(a_x[6])*(b_x[7]-b_x[6]) + inv(a_x[0])*a_x[1]*inv(a_x[2])*a_x[3]*inv(a_x[4])*(b_x[5]-b_x[4]) + inv(a_x[0])*a_x[1]*inv(a_x[2])*(b_x[3]-b_x[2]) + inv(a_x[0])*(b_x[1]-b_x[0])) % P
    #     # a_c_z = a_z[0]*inv(a_z[1])*a_z[2]*inv(a_z[3])*a_z[4]*inv(a_z[5])*a_z[6]*inv(a_z[7]) % P
    #     # b_c_z = (b_z[0]+a_z[0]*inv(a_z[1])*(b_z[2]-b_z[1])+a_z[0]*inv(a_z[1])*a_z[2]*inv(a_z[3])*(b_z[4]-b_z[3])+a_z[0]*inv(a_z[1])*a_z[2]*inv(a_z[3])*a_z[4]*inv(a_z[5])*(b_z[6]-b_z[5])-a_z[0]*inv(a_z[1])*a_z[2]*inv(a_z[3])*a_z[4]*inv(a_z[5])*a_z[6]*inv(a_z[7])*b_z[7]) % P
    #     # is_equiv([a_c_x, b_c_x], function_x)
    #     # is_equiv([a_c_z, b_c_z], function_z)
    #     functions.append(function_x)
    #     functions.append(function_z)
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
    



