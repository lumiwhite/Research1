from constants import *
from z3 import Solver, Int, Sum, sat, Or

def gen_coprime_array():
    # 通常の互いに素な候補（奇数かつ3の倍数でない）
    normal_candidates = [i for i in range(2, P) if math.gcd(i, P) == 1]
    
    # a - 1 が P=768 と大きな公約数を持つ「エリート候補」を抽出
    # （公約数が64以上のものをエリートとする）
    elite_candidates = []
    for a in normal_candidates:
        if math.gcd(a - 1, P) >= 64:
            elite_candidates.append(a)
            
    result = []
    for _ in range(L): # L = 12
        # 例: 80%の確率でエリート候補から、20%の確率で通常候補から選ぶ
        # これにより、解空間を広げつつ、ランダム性（探索の多様性）も維持する
        if random.random() < 0.8 and elite_candidates:
            result.append(random.choice(elite_candidates))
        else:
            result.append(random.choice(normal_candidates))
            
    return result

def gen_cycles(max_len):
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
    for i in range(2, max_len//2+1):
        for r_seq in itertools.product(list(range(J)), repeat=i):
            for c_seq in itertools.product(list(range(L)), repeat=i):
                if is_valid(r_seq, c_seq):
                    positions = get_positions(r_seq, c_seq)
                    cycles.add(tuple(canonicalize(positions)))
    cycles = cycles - gen_utcbc()
    return list(cycles)

def gen_h_xz():
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

def gen_c_constraints(cycles, a_vec, h_x, h_z):
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

        term_x = inv(a_x[0])*a_x[1]
        row_x[idx_x[0]] -= 1
        row_x = [(val * term_x) % P for val in row_x]
        
        for i in range(1, N-1, 2):
            term_x = (inv(a_x[i+1]) * a_x[i+2]) % P
            row_x[idx_x[i]] += 1
            row_x[idx_x[i+1]] -= 1
            row_x = [(val * term_x) % P for val in row_x]
            
        row_x[idx_x[N-1]] += 1
        
        a_c_x = 1
        for i in range(0, N, 2):
            a_c_x = (a_c_x * inv(a_x[i]) * a_x[i+1]) % P
        mul = P // math.gcd(a_c_x-1, P)
        row_x = [(val * mul) % P for val in row_x]
        constraints.append(row_x)

        for i in range(0, N-2, 2):
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
        mul = P // math.gcd(a_c_z-1, P)
        row_z = [(val * mul) % P for val in row_z]
        constraints.append(row_z)
    return constraints

def gen_constraints(Gb, a_vec, cycles, h_x, h_z):
    
    constraints = gen_c_constraints(cycles, a_vec, h_x, h_z)

    # 全ての禁止ベクトル（法ベクトル）を個別にリスト化する
    unique_forbidden_vectors = []
    seen_vectors = set()

    # 1. 条件B (潜在部の非可換性) からの制約 r_i
    for row in Gb:
        row_tuple = tuple(row)
        if row_tuple not in seen_vectors:
            unique_forbidden_vectors.append(row) # 列ベクトルとして保存
            seen_vectors.add(row_tuple)

    # 2. 条件C (短いサイクルの回避) からの制約 c_prime
    for row in constraints:
        row_tuple = tuple(row)
        if row_tuple not in seen_vectors:
            unique_forbidden_vectors.append(row) # 列ベクトルとして保存
            seen_vectors.add(row_tuple)
    return unique_forbidden_vectors


def gen_g_mat(a_vec):
    def gen_g(a_vec):
        G = np.zeros((36, 12), dtype=int)
        for i in range(6):
            for j in range(6):
                G[6*i+j, i] = (1 - a_vec[6+j])
                G[6*i+j, 6+j] = (a_vec[i]-1)
        return G
    G = gen_g(a_vec)
    all_indices = [i for i in range(l_h**2)]
    gb_indices = [3, 8]
    ga_indices = [i for i in all_indices if i not in gb_indices]
    Ga = G[ga_indices]
    Gb = G[gb_indices]
    return Ga, Gb


def find_b_from_random_a(cycles, h_x, h_z, p_val=P):
    attempt = 0
    # seen_avec = []
    while True:
        attempt += 1
        print(f"試行回数: {attempt}")
        # a_vec = gen_coprime_array()
        a_vec = [763, 679, 397, 61, 697, 373, 289, 257, 625, 41, 193, 449]
        # if a_vec in seen_avec:
        #     continue
        # seen_avec.append(a_vec)
        
        Ga, Gb = gen_g_mat(a_vec)
        constraints = gen_constraints(Gb, a_vec, cycles, h_x, h_z)
        solver = Solver()
        solver.set("timeout", 30000)  # タイムアウトを1000msに設定
        # 32bit BitVector を使用
        b = [BitVec(f'b_{i}', 32) for i in range(L)]
        P = BitVecVal(p_val, 32)
        ZERO = BitVecVal(0, 32)

        # 範囲制約
        for x in b:
            solver.add(UGE(x, 0), ULT(x, P))

        # 等式制約 (Ax = 0 mod P)
        for row in Ga:
            # row[i] も BitVecVal に変換して計算
            expr = Sum([BitVecVal(int(row[i]), 32) * b[i] for i in range(L)])
            solver.add(expr % P == ZERO)

        # 不等式制約 (Cx != 0 mod P)
        for row in constraints:
            expr = Sum([BitVecVal(int(row[j]), 32) * b[j] for j in range(L)])
            solver.add(expr % P != ZERO)
            
        if solver.check() == sat:
            model = solver.model()
            return a_vec, [model[x].as_long() for x in b]

def find_b_from_a(a_vec, cycles, h_x, h_z, p_val=P):
        Ga, Gb = gen_g_mat(a_vec)
        constraints = gen_constraints(Gb, a_vec, cycles, h_x, h_z)
        solver = Solver()
        # 32bit BitVector を使用
        b = [BitVec(f'b_{i}', 32) for i in range(L)]
        P = BitVecVal(p_val, 32)
        ZERO = BitVecVal(0, 32)

        # 範囲制約
        for x in b:
            solver.add(UGE(x, 0), ULT(x, P))

        # 等式制約 (Ax = 0 mod P)
        for row in Ga:
            # row[i] も BitVecVal に変換して計算
            expr = Sum([BitVecVal(int(row[i])%p_val, 32) * b[i] for i in range(L)])
            solver.add(expr % P == ZERO)

        # 不等式制約 (Cx != 0 mod P)
        for row in constraints:
            expr = Sum([BitVecVal(int(row[j])%p_val, 32) * b[j] for j in range(L)])
            solver.add(expr % P != ZERO)
            
        res = solver.check()
        if res == sat:
            model = solver.model()
            return a_vec, [model[x].as_long() for x in b]
        else:
            # unsat なのか timeout(unknown) なのかを出力して原因を切り分ける
            print(f"Solver check result: {res}")
            return None





