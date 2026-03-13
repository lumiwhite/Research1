#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <z3.h>

#define L 12
#define L_H 6
#define P 768
#define MAX_SEEN 10000

// ユーティリティ: 最大公約数
int gcd(int a, int b) {
    while (b != 0) {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

// 候補を保持するグローバル配列
int normal_candidates[P];
int num_normal = 0;
int elite_candidates[P];
int num_elite = 0;

// P=768 と互いに素な候補の初期化
void init_candidates() {
    for (int i = 2; i < P; i++) {
        if (gcd(i, P) == 1) {
            normal_candidates[num_normal++] = i;
            // a - 1 が P と大きな公約数を持つ「エリート候補」
            if (gcd(i - 1, P) >= 64) {
                elite_candidates[num_elite++] = i;
            }
        }
    }
}

// 互いに素な配列の生成
void gen_coprime_array(int a_vec[L]) {
    for (int i = 0; i < L; i++) {
        double r = (double)rand() / RAND_MAX;
        if (r < 0.8 && num_elite > 0) {
            a_vec[i] = elite_candidates[rand() % num_elite];
        } else {
            a_vec[i] = normal_candidates[rand() % num_normal];
        }
    }
}

// 行列Gの生成
void gen_g_mat(int a_vec[L], int Ga[34][L], int Gb[2][L]) {
    int G[36][L] = {0};
    for (int i = 0; i < L_H; i++) {
        for (int j = 0; j < L_H; j++) {
            G[L_H * i + j][i] = 1 - a_vec[L_H + j];
            G[L_H * i + j][L_H + j] = a_vec[i] - 1;
        }
    }
    
    int ga_idx = 0;
    int gb_idx = 0;
    for (int i = 0; i < 36; i++) {
        // インデックス 3, 8 を Gb に振り分け
        if (i == 3 || i == 8) {
            for (int j = 0; j < L; j++) Gb[gb_idx][j] = G[i][j];
            gb_idx++;
        } else {
            for (int j = 0; j < L; j++) Ga[ga_idx][j] = G[i][j];
            ga_idx++;
        }
    }
}

// 配列の一致確認
bool arrays_equal(int a[L], int b[L]) {
    for (int i = 0; i < L; i++) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

// Z3を用いた探索
bool find_b_from_random_a(Z3_context ctx, int out_a[L], int out_b[L]) {
    int seen_avec[MAX_SEEN][L];
    int seen_count = 0;
    int attempt = 0;
    
    Z3_sort int_sort = Z3_mk_int_sort(ctx);
    Z3_ast zero = Z3_mk_int(ctx, 0, int_sort);
    Z3_ast p_ast = Z3_mk_int(ctx, P, int_sort);

    while (1) {
        attempt++;
        printf("試行回数: %d\n", attempt);
        
        int a_vec[L];
        gen_coprime_array(a_vec);
        
        // 生成済み配列のチェック
        bool seen = false;
        for (int i = 0; i < seen_count; i++) {
            if (arrays_equal(seen_avec[i], a_vec)) {
                seen = true;
                break;
            }
        }
        if (seen) continue;
        
        if (seen_count < MAX_SEEN) {
            for (int i = 0; i < L; i++) seen_avec[seen_count][i] = a_vec[i];
            seen_count++;
        }
        
        int Ga[34][L];
        int Gb[2][L];
        gen_g_mat(a_vec, Ga, Gb);
        
        // Z3 ソルバーの設定
        Z3_solver solver = Z3_mk_solver(ctx);
        Z3_solver_inc_ref(ctx, solver);
        
        Z3_params p_params = Z3_mk_params(ctx);
        Z3_params_inc_ref(ctx, p_params);
        Z3_symbol timeout_sym = Z3_mk_string_symbol(ctx, "timeout");
        Z3_params_set_uint(ctx, p_params, timeout_sym, 30000);
        Z3_solver_set_params(ctx, solver, p_params);
        Z3_params_dec_ref(ctx, p_params);

        // 変数 b の定義と範囲制約 (0 <= b_i < P)
        Z3_ast b[L];
        for (int i = 0; i < L; i++) {
            char name[10];
            sprintf(name, "b_%d", i);
            Z3_symbol sym = Z3_mk_string_symbol(ctx, name);
            b[i] = Z3_mk_const(ctx, sym, int_sort);
            
            Z3_ast uge = Z3_mk_ge(ctx, b[i], zero);
            Z3_ast ult = Z3_mk_lt(ctx, b[i], p_ast);
            Z3_solver_assert(ctx, solver, uge);
            Z3_solver_assert(ctx, solver, ult);
        }
        
        // 等式制約: Ga * b == 0 mod P
        for (int i = 0; i < 34; i++) {
            Z3_ast terms[L];
            int num_terms = 0;
            for (int j = 0; j < L; j++) {
                if (Ga[i][j] != 0) {
                    Z3_ast coeff = Z3_mk_int(ctx, Ga[i][j], int_sort);
                    Z3_ast mul_args[2] = {coeff, b[j]};
                    terms[num_terms++] = Z3_mk_mul(ctx, 2, mul_args);
                }
            }
            Z3_ast sum;
            if (num_terms == 0) sum = zero;
            else if (num_terms == 1) sum = terms[0];
            else sum = Z3_mk_add(ctx, num_terms, terms);
            
            Z3_ast mod_expr = Z3_mk_mod(ctx, sum, p_ast);
            Z3_ast eq = Z3_mk_eq(ctx, mod_expr, zero);
            Z3_solver_assert(ctx, solver, eq);
        }
        
        // 不等式制約: Gb * b != 0 mod P
        for (int i = 0; i < 2; i++) {
            Z3_ast terms[L];
            int num_terms = 0;
            for (int j = 0; j < L; j++) {
                if (Gb[i][j] != 0) {
                    Z3_ast coeff = Z3_mk_int(ctx, Gb[i][j], int_sort);
                    Z3_ast mul_args[2] = {coeff, b[j]};
                    terms[num_terms++] = Z3_mk_mul(ctx, 2, mul_args);
                }
            }
            Z3_ast sum;
            if (num_terms == 0) sum = zero;
            else if (num_terms == 1) sum = terms[0];
            else sum = Z3_mk_add(ctx, num_terms, terms);
            
            Z3_ast mod_expr = Z3_mk_mod(ctx, sum, p_ast);
            Z3_ast not_eq = Z3_mk_not(ctx, Z3_mk_eq(ctx, mod_expr, zero));
            Z3_solver_assert(ctx, solver, not_eq);
        }
        
        // モデルのチェック
        Z3_lbool res = Z3_solver_check(ctx, solver);
        if (res == Z3_L_TRUE) {
            Z3_model model = Z3_solver_get_model(ctx, solver);
            Z3_model_inc_ref(ctx, model);
            
            for (int i = 0; i < L; i++) {
                Z3_ast v_ast;
                Z3_model_eval(ctx, model, b[i], Z3_TRUE, &v_ast);
                int val = 0;
                Z3_get_numeral_int(ctx, v_ast, &val);
                out_b[i] = val;
                out_a[i] = a_vec[i];
            }
            
            Z3_model_dec_ref(ctx, model);
            Z3_solver_dec_ref(ctx, solver);
            return true;
        }
        
        // メモリリーク防止のため参照カウントを下げる
        Z3_solver_dec_ref(ctx, solver);
    }
}

int main() {
    // 乱数初期化
    srand((unsigned int)time(NULL));
    init_candidates();
    
    // Z3 コンテキストの初期化
    Z3_config cfg = Z3_mk_config();
    Z3_context ctx = Z3_mk_context(cfg);
    Z3_del_config(cfg);
    
    int a_res[L], b_res[L];
    
    printf("探索を開始...\n");
    if (find_b_from_random_a(ctx, a_res, b_res)) {
        printf("解が見つかりました！\n");
        printf("a_vec = [");
        for (int i=0; i<L; i++) printf("%d%s", a_res[i], i==L-1 ? "" : ", ");
        printf("]\n");
        
        printf("b_list = [");
        for (int i=0; i<L; i++) printf("%d%s", b_res[i], i==L-1 ? "" : ", ");
        printf("]\n");
    }
    
    Z3_del_context(ctx);
    return 0;
}