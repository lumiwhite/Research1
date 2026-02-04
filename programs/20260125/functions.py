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
        for r_seq in itertools.permutations(list(range(l_h//2)), 3):
            r0, r1, r2 = r_seq
            for c_seq in itertools.permutations(list(range(L)), 2):
                c0, c1 = c_seq
                if  c0 < l_h:
                    c2 = (r1+c0-r0) % l_h
                else:
                    c2 = l_h + (r1+c0-r0) % l_h
                if  c1 < l_h:
                    c3 = (r2+c1-r1) % l_h
                else:
                    c3 = l_h + (r2+c1-r1) % l_h
                r3 = (c3+r0-c1) % (l_h//2)
                rows=[r0,r1,r2,r3]
                cols=[c0,c1,c2,c3]
                if is_valid(rows, cols):
                    positions = get_positions(rows, cols)
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
        for r_seq in itertools.permutations(list(range(l_h//2)), i):
            for c_seq in itertools.permutations(list(range(L)), i):
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
        constraints.append([row_x, a_c_x])

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
        constraints.append([row_z, a_c_z])
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


