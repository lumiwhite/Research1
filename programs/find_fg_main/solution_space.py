# solution_space.py
import numpy as np
import random
import itertools
from config import P
from math_utils import solve_snf
from matrix_builder import gen_g, gen_g_ab

def extract_snf_basis(y_ranges, V, p_val=P):
    """SNFの出力結果から、解空間を生成するための基底ベクトル群を抽出する。"""
    basis_info = []
    cols = len(V)
    m = len(y_ranges)
    
    for j in range(m):
        step, num_vals = y_ranges[j]
        if num_vals <= 1:
            continue
            
        basis_vec = [0] * cols
        for i in range(cols):
            basis_vec[i] = (V[i][j] * step) % p_val
            
        if any(val != 0 for val in basis_vec):
            basis_info.append({
                'vector': basis_vec,
                'num_vals': num_vals
            })
    return basis_info

def analyze_diff_space_fast(basis_info_Ga, Gb, p_val=P):
    """条件判定に影響する基底を分離し、差空間の有効なパラメータパターンを抽出する。"""
    free_bases = []
    constrained_bases = []
    
    for info in basis_info_Ga:
        vec = np.array(info['vector']) 
        W = np.dot(Gb, vec) % p_val
        
        if np.all(W == 0):
            free_bases.append(info)
        else:
            info['W'] = W
            constrained_bases.append(info)
            
    valid_constrained_patterns = []
    if not constrained_bases:
        return free_bases, constrained_bases, valid_constrained_patterns
        
    c_ranges = [range(info['num_vals']) for info in constrained_bases]
    rows_gb = Gb.shape[0]
    
    for coeffs in itertools.product(*c_ranges):
        current_W = np.zeros(rows_gb, dtype=int)
        for c, info in zip(coeffs, constrained_bases):
            current_W = (current_W + c * info['W']) % p_val
            
        if np.all(current_W != 0):
            valid_constrained_patterns.append(coeffs)
            
    return free_bases, constrained_bases, valid_constrained_patterns

def extract_b_solution_space(a_vec, p_val=P):
    """a_vecから条件A・条件Bを満たすbの解空間を抽出・分離する。"""
    G = gen_g(a_vec)
    Ga, Gb = gen_g_ab(G)
    
    y_ranges_Ga, V_Ga = solve_snf(Ga, p_val)
    ker_Ga_s = extract_snf_basis(y_ranges_Ga, V_Ga, p_val)
    
    free, const, valid_pats = analyze_diff_space_fast(ker_Ga_s, Gb, p_val)
    return free, const, valid_pats

def count_total_solutions(free_bases, valid_constrained_patterns):
    """自由基底と有効な制約パターンから、条件を満たす解の総数を計算する。"""
    total_solutions = 1
    for info in free_bases:
        total_solutions *= info['num_vals']
    total_solutions *= len(valid_constrained_patterns)
    return total_solutions

def generate_solution_from_patterns(free_bases, constrained_bases, valid_patterns, p_val=P, cols=12):
    """特定された自由基底と有効なパターンから、差空間の解をランダムに1つ生成する。"""
    if not valid_patterns:
        return None
        
    b = np.zeros(cols, dtype=int)
    for info in free_bases:
        c = random.randint(0, info['num_vals'] - 1)
        vec = np.array(info['vector'])
        b = (b + c * vec) % p_val
            
    pattern = random.choice(valid_patterns)
    for c, info in zip(pattern, constrained_bases):
        vec = np.array(info['vector'])
        b = (b + c * vec) % p_val
            
    return b