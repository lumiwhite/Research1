from constants import *

def gen_coprime_array():
    candidates = [i for i in range(2, P) if math.gcd(i, P) == 1]
    
    result = random.choices(candidates, k=L)
    return result

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

        term_x = inv(a_x[0])
        for i in range(0, N, 2):
            row_x[idx_x[i]] = (row_x[idx_x[i]] - term_x) % P
            row_x[idx_x[i+1]] = (row_x[idx_x[i+1]] + term_x) % P
            if i + 2 < N:
                term_x = (term_x * a_x[i+1] * inv(a_x[i+2])) % P
        a_c_x = 1
        for i in range(0, N, 2):
            a_c_x = (a_c_x * inv(a_x[i]) * a_x[i+1]) % P
        mul = P // math.gcd(a_c_x-1, P)
        row_x = [(val * mul) % P for val in row_x]
        constraints.append(row_x)

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
        mul = P // math.gcd(a_c_z-1, P)
        row_z = [(val * mul) % P for val in row_z]
        constraints.append(row_z)
    return constraints

def gen_constraints(Gb, a_vec, girth):
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
    cycles = gen_cycles(girth)
    h_x, h_z = gen_h_xz()
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
    print(f"個別に回避すべき禁止制約（超平面）の数: {len(unique_forbidden_vectors)}")
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

def solve_with_bitvec(Ga, constraints, count):
    solver = Solver()
    
    # 32ビットのビットベクトルを定義
    b = [BitVec(f'b_{i}', 32) for i in range(L)]
    
    # 範囲制約: 0 <= b_i < P
    for x in b:
        solver.add(UGE(x, 0))
        solver.add(ULT(x, P))

    # 等式制約 (Ga * b = 0 mod P)
    for row in Ga:
        expr = Sum([BitVecVal(int(row[i]), 32) * b[i] for i in range(L)])
        solver.add(expr % P == 0)

    # 不等式制約 (C * b != 0 mod P)
    for row in constraints:
        expr = Sum([BitVecVal(int(row[j]), 32) * b[j] for j in range(L)])
        solver.add(expr % P != 0)
    
    solutions = []
    
    for _ in range(count):
        if solver.check() == sat:
            model = solver.model()
            # 現在のモデルから解を抽出
            res = [model[x].as_long() for x in b]
            solutions.append(res)
            
            # 少なくとも1つの変数が今の解と異なることを保証する
            solver.add(Or([b[i] != res[i] for i in range(L)]))
        else:
            # これ以上解が存在しない場合はループを抜ける
            break
            
    return solutions












