# cycle_utils.py
import itertools
import math
from config import P, J, L, l_h
from math_utils import composite_affine, func_inv

def gen_cycles(lengths):
    """指定された長さのリスト(例: [4, 6])に対応するサイクル(閉路)を生成する"""
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
    for cycle_len in lengths:
        # 2部グラフの閉路は必ず偶数長になるため、奇数や4未満はスキップする
        if cycle_len % 2 != 0 or cycle_len < 4:
            continue
            
        i = cycle_len // 2
        for r_seq in itertools.product(list(range(J)), repeat=i):
            for c_seq in itertools.product(list(range(L)), repeat=i):
                if is_valid(r_seq, c_seq):
                    positions = get_positions(r_seq, c_seq)
                    cycles.add(tuple(canonicalize(positions)))
                    
    # UTCBCは長さ8の特殊なサイクルであるため、8が含まれる場合のみ生成・除外する
    if 8 in lengths:
        cycles = cycles - gen_utcbc()
        
    return list(cycles)
def get_function_x(cycle, a_vec, b_vec, h_x):
    idx_x = [h_x[r][c] for r, c in cycle]
    a_x = [a_vec[idx] for idx in idx_x]
    b_x = [b_vec[idx] for idx in idx_x]
    function_x = [1, 0]
    for i in range(len(cycle) // 2):
        function_x = composite_affine(function_x, [a_x[2*i], b_x[2*i]])
        function_x = composite_affine(function_x, func_inv([a_x[2*i+1], b_x[2*i+1]]))
    return function_x

def get_function_z(cycle, a_vec, b_vec, h_z):
    idx_z = [h_z[r][c] for r, c in cycle]
    a_z = [a_vec[idx] for idx in idx_z]
    b_z = [b_vec[idx] for idx in idx_z]
    function_z = [1, 0]
    for i in range(len(cycle) // 2):
        function_z = composite_affine(function_z, [a_z[2*i], b_z[2*i]])
        function_z = composite_affine(function_z, func_inv([a_z[2*i+1], b_z[2*i+1]]))
    return function_z

def get_functions(cycle, a_vec, b_vec, h_x, h_z):
    f_x = get_function_x(cycle, a_vec, b_vec, h_x)
    f_z = get_function_z(cycle, a_vec, b_vec, h_z)
    return f_x, f_z

def is_closed(input_func):
    """置換関数が閉路を形成しているか(恒等写像に近いか)を判定する"""
    a = input_func[0]
    b = input_func[1]
    if a == 1 and b == 0:
        return True
    d = math.gcd(a - 1, P)
    if b % d == 0:
        return True
    return False

def count_cycles(a_vec, b_vec):
    """与えられたa, bのベクトルにおいて、形成されるサイクル数をカウントする"""
    from matrix_builder import gen_h_xz  # 循環インポート回避のためのローカルインポート
    
    h_x, h_z = gen_h_xz()
    cycles = gen_cycles(8)
    result = {}
    for cycle in cycles:
        length = len(cycle)
        f_x, f_z = get_functions(cycle, a_vec, b_vec, h_x, h_z)
        if is_closed(f_x) or is_closed(f_z):
            if length not in result:
                result[length] = 1
            else:
                result[length] += 1
    return result