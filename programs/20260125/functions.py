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


def inv(val):
    try:
        return pow(val, -1, P)
    except ValueError:
        raise ValueError("Inverse does not exist")


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


