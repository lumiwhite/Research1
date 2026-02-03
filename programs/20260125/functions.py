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

def is_valid(positions):
    length = len(positions)
    for i in range(length):
        if positions[i][0] == positions[(i+1)%length][0] and positions[i][1] == positions[(i+1)%length][1]:
            return False
        if positions[i][0] != positions[(i+1)%length][0] and positions[i][1] != positions[(i+1)%length][1]:
            return False
    return True
    

def generate_cycles():
    def generate_utcbc():
        utcbcs = set()
        for r_seq in itertools.combinations(list(range(l_h)), 3):
            r0, r1, r2 = r_seq
            for c_seq in itertools.combinations(list(range(L)), 2):
                c0, c1 = c_seq
                if  c0 < l_h:
                    c2 = (r1+c0-r0) % l_h
                else:
                    c2 = l_h + (r1+c0-r0) % l_h
                if  c1 < l_h:
                    c3 = (r2+c1-r1) % l_h
                else:
                    c3 = l_h + (r2+c1-r1) % l_h
                r3 = (c3+r0-c1) % l_h
                rows=[r0,r1,r2,r3]
                cols=[c0,c1,c2,c3]
                positions = get_positions(rows, cols)
                if is_valid(positions):
                    utcbcs.add(tuple(canonicalize(positions)))
        return utcbcs

    def get_positions(r_seq, c_seq):
        length = len(r_seq)
        pos = []
        for i in range(length):
            # [] ではなく () を使うことで Hashable にする
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
    for r_seq in itertools.combinations(list(range(l_h)), 2):
        for c_seq in itertools.combinations(list(range(L)), 2):
            positions = get_positions(r_seq, c_seq)
            if is_valid(positions):
                cycles.add(tuple(canonicalize(positions)))
                
    for r_seq in itertools.combinations(list(range(l_h)), 3):
        for c_seq in itertools.combinations(list(range(L)), 3):
            positions = get_positions(r_seq, c_seq)
            if is_valid(positions):
                cycles.add(tuple(canonicalize(positions)))

    for r_seq in itertools.combinations(list(range(l_h)), 4):
        for c_seq in itertools.combinations(list(range(L)), 4):
            positions = get_positions(r_seq, c_seq)
            if is_valid(positions):
                cycles.add(tuple(canonicalize(positions)))
    cycles = cycles - generate_utcbc()
    return list(cycles)

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
        row_x = [0]*L
        row_z = [0]*L
        idx_x = []
        idx_z = []
        if len(cycle) == 4:
            for i in range(4):
                idx_x.append(h_x[cycle[i][0]][cycle[i][1]])
                idx_z.append(h_z[cycle[i][0]][cycle[i][1]])
            a_x = [a_vec[idx] for idx in idx_x]
            a_z = [a_vec[idx] for idx in idx_z]
            row_x[idx_x[0]] -= inv(a_x[0])
            row_x[idx_x[1]] += inv(a_x[0])
            row_x[idx_x[2]] -= inv(a_x[0]) * a_x[1] * inv(a_x[2])
            row_x[idx_x[3]] += inv(a_x[0]) * a_x[1] * inv(a_x[2])
            a_c_x = inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] % P
            row_x = [val % P for val in row_x]
            row_z[idx_z[0]] += 1
            row_z[idx_z[1]] -= a_z[0] * inv(a_z[1])
            row_z[idx_z[2]] += a_z[0] * inv(a_z[1])
            row_z[idx_z[3]] -= a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3])
            a_c_z = a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3]) % P
            row_z = [val % P for val in row_z]
            constraints.append([row_x, a_c_x])
            constraints.append([row_z, a_c_z]) 
        elif len(cycle) == 6:
            for i in range(6):
                idx_x.append(h_x[cycle[i][0]][cycle[i][1]])
                idx_z.append(h_z[cycle[i][0]][cycle[i][1]])
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
        elif len(cycle) == 8:
            for i in range(8):
                idx_x.append(h_x[cycle[i][0]][cycle[i][1]])
                idx_z.append(h_z[cycle[i][0]][cycle[i][1]])
            a_x = [a_vec[idx] for idx in idx_x]
            a_z = [a_vec[idx] for idx in idx_z]
            row_x[idx_x[0]] -= inv(a_x[0])
            row_x[idx_x[1]] += inv(a_x[0])
            row_x[idx_x[2]] -= inv(a_x[0]) * a_x[1] * inv(a_x[2])
            row_x[idx_x[3]] += inv(a_x[0]) * a_x[1] * inv(a_x[2])
            row_x[idx_x[4]] -= inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4])
            row_x[idx_x[5]] += inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4])
            row_x[idx_x[6]] -= inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4]) * a_x[5] * inv(a_x[6])
            row_x[idx_x[7]] += inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4]) * a_x[5] * inv(a_x[6])  
            a_c_x = inv(a_x[0]) * a_x[1] * inv(a_x[2]) * a_x[3] * inv(a_x[4]) * a_x[5] * inv(a_x[6]) * a_x[7] % P
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
            a_c_z = a_z[0] * inv(a_z[1]) * a_z[2] * inv(a_z[3]) * a_z[4] * inv(a_z[5]) * a_z[6] * inv(a_z[7]) % P
            row_z = [val % P for val in row_z]
            constraints.append([row_z, a_c_z])
    return constraints

def is_equiv(func1, func2):
    if (func1 != func2):
        print(False)
    return 0

def generate_functions(cycles, a_vec, b_vec, h_x, h_z):
    functions = []
    for cycle in cycles:
        idx_x = []
        idx_z = []
        function_x = [1, 0]
        function_z = [1, 0]
        if len(cycle) == 4:
            for i in range(4):
                idx_x.append(h_x[cycle[i][0]][cycle[i][1]])
                idx_z.append(h_z[cycle[i][0]][cycle[i][1]])
            a_x = [a_vec[idx] for idx in idx_x]
            b_x = [b_vec[idx] for idx in idx_x]
            a_z = [a_vec[idx] for idx in idx_z]
            b_z = [b_vec[idx] for idx in idx_z]
            for i in range(2):
                function_x = composite_affine(function_x, func_inv([a_x[2*i], b_x[2*i]]))
                function_x = composite_affine(function_x, [a_x[2*i+1], b_x[2*i+1]])
                function_z = composite_affine(function_z, [a_z[2*i], b_z[2*i]])
                function_z = composite_affine(function_z, func_inv([a_z[2*i+1], b_z[2*i+1]]))
            functions.append(function_x)
            functions.append(function_z)
        elif len(cycle) == 6:
            for i in range(6):
                idx_x.append(h_x[cycle[i][0]][cycle[i][1]])
                idx_z.append(h_z[cycle[i][0]][cycle[i][1]])
            a_x = [a_vec[idx] for idx in idx_x]
            b_x = [b_vec[idx] for idx in idx_x]
            a_z = [a_vec[idx] for idx in idx_z]
            b_z = [b_vec[idx] for idx in idx_z]
            for i in range(3):
                function_x = composite_affine(function_x, func_inv([a_x[2*i], b_x[2*i]]))
                function_x = composite_affine(function_x, [a_x[2*i+1], b_x[2*i+1]])
                function_z = composite_affine(function_z, [a_z[2*i], b_z[2*i]])
                function_z = composite_affine(function_z, func_inv([a_z[2*i+1], b_z[2*i+1]]))
            functions.append(function_x)
            functions.append(function_z)
        elif len(cycle) == 8:
            for i in range(8):
                idx_x.append(h_x[cycle[i][0]][cycle[i][1]])
                idx_z.append(h_z[cycle[i][0]][cycle[i][1]])
            a_x = [a_vec[idx] for idx in idx_x]
            b_x = [b_vec[idx] for idx in idx_x]
            a_z = [a_vec[idx] for idx in idx_z]
            b_z = [b_vec[idx] for idx in idx_z]
            for i in range(4):
                function_x = composite_affine(function_x, func_inv([a_x[2*i], b_x[2*i]]))
                function_x = composite_affine(function_x, [a_x[2*i+1], b_x[2*i+1]])
                function_z = composite_affine(function_z, [a_z[2*i], b_z[2*i]])
                function_z = composite_affine(function_z, func_inv([a_z[2*i+1], b_z[2*i+1]]))
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
    



