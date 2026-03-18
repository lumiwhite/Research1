# z3_solver.py
import numpy as np
from z3 import *
from config import P, L
from matrix_builder import gen_g_mat, gen_constraints

def solve_z3_for_pattern(free_bases, constrained_bases, pattern, constraints, p_val=P, cols=L):
    """固定された制約パターンのもとで、Z3を用いて条件Cを満たすbベクトルを見つける。"""
    solver = Solver()
    
    P_bv = BitVecVal(p_val, 32)
    ZERO = BitVecVal(0, 32)
    
    # 1. 探索変数: 自由基底の係数 c_i
    c_vars = []
    for i, info in enumerate(free_bases):
        c = BitVec(f'c_{i}', 32)
        solver.add(UGE(c, 0), ULT(c, info['num_vals']))
        c_vars.append(c)
        
    # 2. 中間変数: 最終的なシフトベクトル b_k
    b_vars = [BitVec(f'b_{k}', 32) for k in range(cols)]
    for b in b_vars:
        solver.add(UGE(b, 0), ULT(b, p_val))
        
    # 3. 制約基底から生成される定数ベクトルを計算
    const_b = np.zeros(cols, dtype=int)
    for val, info in zip(pattern, constrained_bases):
        const_b += val * np.array(info['vector'])
        
    # 4. 基底の線形結合と b_vars を結びつける制約式
    for k in range(cols):
        expr_terms = [BitVecVal(int(const_b[k]) % p_val, 32)]
        
        for i, info in enumerate(free_bases):
            v_val = int(info['vector'][k]) % p_val
            if v_val != 0:
                expr_terms.append(BitVecVal(v_val, 32) * c_vars[i])
                
        linear_expr = Sum(expr_terms) if len(expr_terms) > 1 else expr_terms[0]
        solver.add(b_vars[k] == URem(linear_expr, P_bv))
        
    # 5. 条件C (ガース制約)
    for row in constraints:
        row_terms = []
        for j in range(cols):
            val = int(row[j]) % p_val
            if val != 0:
                row_terms.append(BitVecVal(val, 32) * b_vars[j])
                
        if row_terms:
            row_expr = Sum(row_terms) if len(row_terms) > 1 else row_terms[0]
            solver.add(URem(row_expr, P_bv) != ZERO)
            
    # 6. ソルバの実行
    res = solver.check()
    if res == sat:
        model = solver.model()
        return np.array([model[b].as_long() for b in b_vars])
        
    return None

def find_b_with_z3(free_bases, constrained_bases, valid_patterns, constraints, p_val=P, cols=L):
    """有効な制約パターンをループしてZ3ソルバを回す。"""
    for i, pattern in enumerate(valid_patterns):
        b_sol = solve_z3_for_pattern(free_bases, constrained_bases, pattern, constraints, p_val, cols)
        
        if b_sol is not None:
            return b_sol
            
        print(f"\rパターン {i+1}/{len(valid_patterns)} は条件Cを満たす解なし...", end="")
        
    print("\n全てのパターンを探索したが、ガース条件を満たす解は存在しなかった。")
    return None

def find_b_from_a(a_vec, cycles, h_x, h_z, p_val=P):
    """(旧機能) a_vecからZ3を用いて直接bを探索するフォールバック関数。"""
    Ga, Gb = gen_g_mat(a_vec)
    constraints = gen_constraints(Gb, a_vec, cycles, h_x, h_z)
    solver = Solver()
    
    b = [BitVec(f'b_{i}', 32) for i in range(L)]
    P_bv = BitVecVal(p_val, 32)
    ZERO = BitVecVal(0, 32)

    for x in b:
        solver.add(UGE(x, 0), ULT(x, P_bv))

    for row in Ga:
        expr = Sum([BitVecVal(int(row[i])%p_val, 32) * b[i] for i in range(L)])
        solver.add(expr % P_bv == ZERO)

    for row in constraints:
        expr = Sum([BitVecVal(int(row[j])%p_val, 32) * b[j] for j in range(L)])
        solver.add(expr % P_bv != ZERO)
        
    res = solver.check()
    if res == sat:
        model = solver.model()
        return a_vec, [model[x].as_long() for x in b]
    else:
        print(f"Solver check result: {res}")
        return None
    

def solve_z3_for_pattern(free_bases, constrained_bases, pattern, constraints, P, cols=12):
    """
    固定された制約パターン(条件B充足)のもとで、
    Z3(BitVec)を用いて自由基底の係数を探索し、ガース条件を満たす b ベクトルを見つける。
    """
    solver = Solver()
    
    P_bv = BitVecVal(P, 32)
    ZERO = BitVecVal(0, 32)
    
    # 1. 探索変数: 自由基底の係数 c_i (0 から num_vals - 1 の整数)
    c_vars = []
    for i, info in enumerate(free_bases):
        c = BitVec(f'c_{i}', 32)
        # BitVecの符号なし比較演算子を使用
        solver.add(UGE(c, 0), ULT(c, info['num_vals']))
        c_vars.append(c)
        
    # 2. 中間変数: 最終的なシフトベクトル b_k (0 から P - 1 の整数)
    b_vars = [BitVec(f'b_{k}', 32) for k in range(cols)]
    for b in b_vars:
        solver.add(UGE(b, 0), ULT(b, P))
        
    # 3. 制約基底から生成される定数ベクトルを計算 (条件Bを担保)
    const_b = np.zeros(cols, dtype=int)
    for val, info in zip(pattern, constrained_bases):
        const_b += val * np.array(info['vector'])
        
    # 4. 基底の線形結合と b_vars を mod P で結びつける制約式
    for k in range(cols):
        # 定数項 const_b[k]
        expr_terms = [BitVecVal(int(const_b[k]) % P, 32)]
        
        # 変数項 c_i * v_{ik}
        for i, info in enumerate(free_bases):
            v_val = int(info['vector'][k]) % P
            if v_val != 0: # 係数が0の項はZ3のパース時間を節約するため省く
                expr_terms.append(BitVecVal(v_val, 32) * c_vars[i])
                
        # sum() で足し合わせ、符号なしモジュロ(URem)で b_k と等値化する
        linear_expr = Sum(expr_terms) if len(expr_terms) > 1 else expr_terms[0]
        solver.add(b_vars[k] == URem(linear_expr, P_bv))
        
    # 5. 条件C (ガース制約) の不等式制約を追加 (Cx != 0 mod P)
    for row in constraints:
        row_terms = []
        for j in range(cols):
            val = int(row[j]) % P
            if val != 0:
                row_terms.append(BitVecVal(val, 32) * b_vars[j])
                
        if row_terms:
            row_expr = Sum(row_terms) if len(row_terms) > 1 else row_terms[0]
            solver.add(URem(row_expr, P_bv) != ZERO)
            
    # 6. ソルバの実行
    res = solver.check()
    if res == sat:
        model = solver.model()
        # 解が見つかった場合、b ベクトルを抽出してNumPy配列で返す
        return np.array([model[b].as_long() for b in b_vars])
    elif res == unknown:
        # タイムアウト等が発生した場合は通知
        # print("Solver check result: unknown (timeout)")
        pass
        
    return None

def find_b_with_z3(free_bases, constrained_bases, valid_patterns, constraints, P, cols=12):
    """
    有効な制約パターンをループしてZ3ソルバを回し、すべての条件を満たす b を見つける。
    """
    # print(f"探索開始: {len(valid_patterns)} 個の有効パターンについてZ3を実行する。")
    
    for i, pattern in enumerate(valid_patterns):
        b_sol = solve_z3_for_pattern(free_bases, constrained_bases, pattern, constraints, P, cols)
        
        if b_sol is not None:
            # print(f"\nパターン {i+1}/{len(valid_patterns)} で解を発見した。")
            return b_sol
            
        # 実行進捗を同じ行に上書き表示
        print(f"\rパターン {i+1}/{len(valid_patterns)} は条件Cを満たす解なし...", end="")
        
    print("\n全てのパターンを探索したが、ガース条件を満たす解は存在しなかった。")
    return None