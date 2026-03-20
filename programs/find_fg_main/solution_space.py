# solution_space.py
import numpy as np
import random
import itertools
from config import P
from math_utils import solve_snf
from matrix_builder import gen_g, gen_g_ab
import math

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

def analyze_diff_space_fast(basis_info_Ga, Gb_combined, p_val=P):
    """
    条件判定に影響する基底を分離し、差空間の有効なパラメータパターンを抽出する。
    Gb_combined には、条件Bと条件C(長さ4など)の禁止制約がすべて含まれる。
    """
    free_bases = []
    constrained_bases = []
    
    for info in basis_info_Ga:
        vec = np.array(info['vector']) 
        W = np.dot(Gb_combined, vec) % p_val
        
        # すべての制約に対して 0 になるなら、完全に自由な基底
        if np.all(W == 0):
            free_bases.append(info)
        else:
            info['W'] = W
            constrained_bases.append(info)
            
    valid_constrained_patterns = []
    if not constrained_bases:
        return free_bases, constrained_bases, valid_constrained_patterns
        
    c_ranges = [range(info['num_vals']) for info in constrained_bases]
    rows_gb = Gb_combined.shape[0]
    
    for coeffs in itertools.product(*c_ranges):
        current_W = np.zeros(rows_gb, dtype=int)
        for c, info in zip(coeffs, constrained_bases):
            current_W = (current_W + c * info['W']) % p_val
            
        # 全ての禁止制約について、計算結果が 0 でない (not 0 mod P) ことを確認
        if np.all(current_W != 0):
            valid_constrained_patterns.append(coeffs)
            
    return free_bases, constrained_bases, valid_constrained_patterns

def extract_b_solution_space(a_vec, constraints_4=None, p_val=P):
    """
    a_vecから条件Aを満たし、かつ条件Bと条件C(長さ4)をクリアするbの解空間を抽出・分離する。
    """
    G = gen_g(a_vec)
    Ga, Gb = gen_g_ab(G)
    
    # 長さ4の制約が渡された場合、Gbと縦に結合して1つの厳しい制約行列にする
    if constraints_4 is not None and len(constraints_4) > 0:
        Gb_combined = np.vstack((Gb, np.array(constraints_4)))
    else:
        Gb_combined = Gb
    
    y_ranges_Ga, V_Ga = solve_snf(Ga, p_val)
    ker_Ga_s = extract_snf_basis(y_ranges_Ga, V_Ga, p_val)
    
    # 結合した制約行列で空間を解析する
    free, const, valid_pats = analyze_diff_space_fast(ker_Ga_s, Gb_combined, p_val)
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

def find_b_with_numpy(free_bases, constrained_bases, valid_patterns, constraints_all, P, cols=12):
    """
    NumPyの超並列行列計算を用いて、数百万の解空間をバッチ処理でテストする。
    メモリ爆発を防ぐため、候補を分割して計算する。
    """
    import numpy as np
    import itertools
    import random

    # 制約行列の準備 (メモリ節約と高速化のため int32 を使用)
    C_all = np.array(constraints_all, dtype=np.int32)
    if len(C_all) == 0:
        C_all = np.zeros((1, cols), dtype=np.int32)

    # 1. 自由基底の全係数パターン(c_0, c_1...)を生成
    if free_bases:
        c_ranges = [range(info['num_vals']) for info in free_bases]
        all_c_coeffs = list(itertools.product(*c_ranges))
    else:
        all_c_coeffs = [tuple()]

    # 候補が多すぎる場合はランダムにサンプリング (今回は上限を広めに設定)
    if len(all_c_coeffs) > 500000:
        all_c_coeffs = random.sample(all_c_coeffs, 500000)

    # 2. 自由基底によるベース空間を構築
    if free_bases:
        free_vecs = np.array([info['vector'] for info in free_bases], dtype=np.int32)
        c_matrix = np.array(all_c_coeffs, dtype=np.int32)
        B_free = (c_matrix @ free_vecs) % P
    else:
        B_free = np.zeros((1, cols), dtype=np.int32)

    # 1回の行列積で処理する候補の数 (メモリと速度のバランスを取る)
    BATCH_SIZE = 10000

    # 3. 制約パターンのループ
    for pattern in valid_patterns:
        const_b = np.zeros(cols, dtype=np.int32)
        for val, info in zip(pattern, constrained_bases):
            const_b = (const_b + val * np.array(info['vector'], dtype=np.int32)) % P

        # 候補となる全ての b ベクトルを生成
        B_candidates = (B_free + const_b) % P
        num_candidates = B_candidates.shape[0]

        # 4. メモリ爆発を防ぐためのバッチ処理
        for i in range(0, num_candidates, BATCH_SIZE):
            B_batch = B_candidates[i:i + BATCH_SIZE]

            # W_batch の shape = (BATCH_SIZE, 制約の数)
            W_batch = (B_batch @ C_all.T) % P

            # すべての制約に対して計算結果が 0 でないかを一括判定
            valid_mask = np.all(W_batch != 0, axis=1)

            # 1つでも条件をクリアした候補があれば、最初のものを返す
            if np.any(valid_mask):
                first_valid_idx = np.argmax(valid_mask)
                return B_batch[first_valid_idx]

    return None

def get_snf_basis(a_vec, p_val):
    """条件A(可換性)のみを満たす解空間の基底を抽出する"""
    from matrix_builder import gen_g, gen_g_ab
    from math_utils import solve_snf
    
    G = gen_g(a_vec)
    Ga, Gb = gen_g_ab(G)
    
    y_ranges_Ga, V_Ga = solve_snf(Ga, p_val)
    ker_Ga_s = extract_snf_basis(y_ranges_Ga, V_Ga, p_val)
    
    return ker_Ga_s, Gb

def find_b_incremental_numpy(ker_Ga_s, Gb, constraints_4, constraints_6, p_val, cols=12):
    """
    ジェネレータとNumPyを用いて、メモリ爆発を防ぎながら数百万の解空間を「段階的」に全探索する。
    """
    if not ker_Ga_s:
        return None
        
    V = np.array([info['vector'] for info in ker_Ga_s], dtype=np.int32)
    num_vals = [info['num_vals'] for info in ker_Ga_s]
    
    # 第1段階の制約: 条件B (Gb) と 長さ4の制約 (C4) を結合
    C4_all = np.array(Gb.tolist() + constraints_4, dtype=np.int32)
    
    # 第2段階の制約: 長さ6の制約 (C6)
    C6_all = np.array(constraints_6, dtype=np.int32)
    
    c_ranges = [range(nv) for nv in num_vals]
    total_space = math.prod(num_vals)
    
    # メモリを節約しつつGPU/CPU並列を活かせる最適なバッチサイズ (5万件)
    BATCH_SIZE = 50000 
    
    # イテレータを作成（全組み合わせをメモリに展開せず、順番に生成する準備）
    coeff_iterator = itertools.product(*c_ranges)
    
    # 空間の広さが何百万あっても、5万件ずつ切り出して順次処理する
    for _ in range(0, total_space, BATCH_SIZE):
        # ジェネレータから BATCH_SIZE 分だけ係数の組み合わせを取り出す
        batch_tuples = list(itertools.islice(coeff_iterator, BATCH_SIZE))
        if not batch_tuples:
            break
            
        c_batch = np.array(batch_tuples, dtype=np.int32)
        B_batch = (c_batch @ V) % p_val
        
        # ==================================================
        # 【第1段階】条件B と 長さ4の制約 でバッチを一括ふるい落とし
        # ==================================================
        W4 = (B_batch @ C4_all.T) % p_val
        valid_4_mask = np.all(W4 != 0, axis=1)
        
        B_passed_4 = B_batch[valid_4_mask]
        
        # 5万件のうち、長さ4をクリアしたものが1つもなければ次の5万件へ
        if len(B_passed_4) == 0:
            continue  
            
        # ==================================================
        # 【第2段階】長さ4をクリアした少数精鋭だけを長さ6でテスト
        # ==================================================
        if len(C6_all) > 0:
            W6 = (B_passed_4 @ C6_all.T) % p_val
            valid_6_mask = np.all(W6 != 0, axis=1)
        else:
            # C6の制約がない場合はすべてクリア扱い
            valid_6_mask = np.ones(len(B_passed_4), dtype=bool)
        
        # 全ての条件をクリアした解が1つでも見つかれば、即座に終了して返す
        if np.any(valid_6_mask):
            first_valid_idx = np.argmax(valid_6_mask)
            return B_passed_4[first_valid_idx]
            
    # 全空間を探索しても見つからなかった場合
    return None