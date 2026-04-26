#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include "cycles_data012.h"

// 実行時にコマンドライン引数から設定するため、グローバル変数に変更
int P;
#define L 12
#define L_H 6

/*
 * 改良および修正点:
 * 1. Pを実行時引数で受け取るように変更。
 * 2. 線形部aの完全全探索: FIXED_Aを削除し、条件を満たす全てのaを動的に探索。
 * 3. Fail-First順序の維持: 厳しい制約から先に処理し、探索空間の爆発を抑える。
 * 4. 動的制約充足: 選んだaに対して、既存のすべての決定済み変数との合同式を解き、bを絞り込む。
 */

// 探索順序: F0(0), G3(9), F1(1), G2(8)...
const int search_order[L] = {0, 9, 1, 8, 2, 7, 3, 6, 4, 11, 5, 10};
int idx_to_step[L];

const Cycle** active_cycles_x[L];
const Cycle** active_cycles_z[L];
int num_active_cycles_x[L] = {0};
int num_active_cycles_z[L] = {0};

typedef struct {
    int a;
    int b;
} Pair;

typedef struct {
    Pair *cands;
    int num_cands;
    int current_cand_idx;
} SearchState;

inline int mod(int a, int b) {
    int r = a % b;
    return r < 0 ? r + b : r;
}

int gcd(int a, int b) {
    a = abs(a); b = abs(b);
    while (b) {
        a %= b;
        int tmp = a; a = b; b = tmp;
    }
    return a;
}

int mod_inverse(int a, int m) {
    int m0 = m, t, q;
    int x0 = 0, x1 = 1;
    if (m == 1) return 0;
    while (a > 1) {
        q = a / m;
        t = m;
        m = mod(a, m); a = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }
    if (x1 < 0) x1 += m0;
    return x1;
}

int solve_linear_congruence(int A, int B, int M, int *results) {
    A = mod(A, M);
    B = mod(B, M);
    int g = gcd(A, M);
    if (B % g != 0) return 0; 

    int a_prime = A / g;
    int b_prime = B / g;
    int m_prime = M / g;
    int x0 = mod(mod_inverse(a_prime, m_prime) * b_prime, m_prime);

    for (int k = 0; k < g; k++) {
        results[k] = x0 + k * m_prime;
    }
    return g;
}

inline void composite(int a1, int b1, int a2, int b2, int *out_a, int *out_b) {
    *out_a = mod(a1 * a2, P);
    *out_b = mod(a1 * b2 + b1, P);
}

inline void func_inv(int a, int b, int *out_a, int *out_b) {
    int inv = mod_inverse(a, P);
    *out_a = inv;
    *out_b = mod(-inv * b, P);
}

inline bool is_closed(int a, int b) {
    if (a == 1) return b == 0;
    return mod(b, gcd(a - 1, P)) == 0;
}

bool check_cycles(int step, const int *a_vec, const int *b_vec) {
    for(int c = 0; c < num_active_cycles_x[step]; c++) {
        const Cycle* cyc = active_cycles_x[step][c];
        int len = cyc->len;
        const int *nodes = cyc->nodes;
        int c_a = 1, c_b = 0;
        for(int i = 0; i < len; i++) {
            int n = nodes[i];
            int f_a = a_vec[n];
            int f_b = b_vec[n];
            if(i % 2 == 0) composite(c_a, c_b, f_a, f_b, &c_a, &c_b);
            else {
                int inv_a, inv_b;
                func_inv(f_a, f_b, &inv_a, &inv_b);
                composite(c_a, c_b, inv_a, inv_b, &c_a, &c_b);
            }
        }
        if(is_closed(c_a, c_b)) return false;
    }
    for(int c = 0; c < num_active_cycles_z[step]; c++) {
        const Cycle* cyc = active_cycles_z[step][c];
        int len = cyc->len;
        const int *nodes = cyc->nodes;
        int c_a = 1, c_b = 0;
        for(int i = 0; i < len; i++) {
            int n = nodes[i];
            int f_a = a_vec[n];
            int f_b = b_vec[n];
            if(i % 2 == 1) composite(c_a, c_b, f_a, f_b, &c_a, &c_b);
            else {
                int inv_a, inv_b;
                func_inv(f_a, f_b, &inv_a, &inv_b);
                composite(c_a, c_b, inv_a, inv_b, &c_a, &c_b);
            }
        }
        if(is_closed(c_a, c_b)) return false;
    }
    return true;
}

unsigned int xorshift32(unsigned int *state) {
    unsigned int x = *state;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    return *state = x;
}

void shuffle_pairs(Pair *array, int n, unsigned int *state) {
    if (n > 1) {
        for (int i = 0; i < n - 1; i++) {
            int j = i + xorshift32(state) % (n - i);
            Pair t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}

void init_state(int step, const int *a_vec, const int *b_vec, SearchState *state, unsigned int *seed, const int *valid_a_list, int num_valid_a) {
    state->num_cands = 0;
    state->current_cand_idx = 0;
    
    int cur_idx = search_order[step];

    // 全ての有効なaについて探索
    for (int a_i = 0; a_i < num_valid_a; a_i++) {
        int a_val = valid_a_list[a_i];
        bool b_valid[P];
        for(int k=0; k<P; k++) b_valid[k] = true;

        for(int s=0; s<step; s++) {
            int prev_idx = search_order[s];
            
            // F同士、G同士の場合は可換条件の制約を受けない
            if ((cur_idx < L_H && prev_idx < L_H) || (cur_idx >= L_H && prev_idx >= L_H)) {
                continue;
            }

            int i, j;
            int a_F, a_G, b_F_known, b_G_known;
            bool solving_for_F = (cur_idx < L_H);

            if (solving_for_F) {
                i = cur_idx;
                j = prev_idx - L_H;
                a_F = a_val;
                a_G = a_vec[prev_idx];
                b_G_known = b_vec[prev_idx];
            } else {
                i = prev_idx;
                j = cur_idx - L_H;
                a_F = a_vec[prev_idx];
                a_G = a_val;
                b_F_known = b_vec[prev_idx];
            }

            int r = (i + j) % L_H;
            bool commute_required = true;
            if (r == 3 && ((i == 0 && j == 3) || (i == 1 && j == 2))) {
                commute_required = false; 
            }

            // Ax = B (mod P) の構築
            int A, B;
            if (solving_for_F) {
                A = a_G - 1;
                B = (a_F - 1) * b_G_known;
            } else {
                A = a_F - 1;
                B = (a_G - 1) * b_F_known;
            }

            int sol_buffer[P];
            int num_sol = solve_linear_congruence(A, B, P, sol_buffer);

            if (commute_required) {
                bool current_mask[P];
                for(int k=0; k<P; k++) current_mask[k] = false; // ループで明示的に初期化
                for(int k=0; k<num_sol; k++) current_mask[sol_buffer[k]] = true;
                for(int k=0; k<P; k++) b_valid[k] &= current_mask[k];
            } else {
                for(int k=0; k<num_sol; k++) b_valid[sol_buffer[k]] = false;
            }
        }

        for(int d=0; d<P; d++) {
            if(b_valid[d] && state->num_cands < 200000) {
                state->cands[state->num_cands].a = a_val;
                state->cands[state->num_cands].b = d;
                state->num_cands++;
            }
        }
    }
    shuffle_pairs(state->cands, state->num_cands, seed);
}

int main(int argc, char *argv[]) {
    // 実行時引数のチェック
    if (argc != 2) {
        fprintf(stderr, "使用法: %s <Pの値>\n", argv[0]);
        fprintf(stderr, "例: %s 384\n", argv[0]);
        return 1;
    }

    P = atoi(argv[1]);
    if (P <= 0) {
        fprintf(stderr, "エラー: Pは正の整数である必要がある。\n");
        return 1;
    }

    int num_workers = omp_get_max_threads();
    // printf("P = %d, %d 個の並列プロセスで完全全探索を開始 (Fail-First順序、線形部aの固定解除)...\n", P, num_workers);
    fflush(stdout);

    int valid_a_list[P];
    int num_valid_a = 0;
    for(int a = 1; a < P; a++) if(gcd(a, P) == 1) valid_a_list[num_valid_a++] = a;

    for (int step = 0; step < L; step++) {
        idx_to_step[search_order[step]] = step;
    }
    for (int i=0; i<L; i++) {
        active_cycles_x[i] = (const Cycle**)malloc(sizeof(const Cycle*) * 100000);
        active_cycles_z[i] = (const Cycle**)malloc(sizeof(const Cycle*) * 100000);
    }
    
    for (int orig_idx = 0; orig_idx < L; orig_idx++) {
        for (int c = 0; c < num_cycles_x_by_idx[orig_idx]; c++) {
            const Cycle* cyc = &cycles_x_by_idx[orig_idx][c];
            int max_step = -1;
            for (int i = 0; i < cyc->len; i++) {
                int step = idx_to_step[cyc->nodes[i]];
                if (step > max_step) max_step = step;
            }
            active_cycles_x[max_step][num_active_cycles_x[max_step]++] = cyc;
        }
        for (int c = 0; c < num_cycles_z_by_idx[orig_idx]; c++) {
            const Cycle* cyc = &cycles_z_by_idx[orig_idx][c];
            int max_step = -1;
            for (int i = 0; i < cyc->len; i++) {
                int step = idx_to_step[cyc->nodes[i]];
                if (step > max_step) max_step = step;
            }
            active_cycles_z[max_step][num_active_cycles_z[max_step]++] = cyc;
        }
    }
    
    volatile bool found = false;
    int sol_a[L], sol_b[L];
    double start_time = omp_get_wtime();
    
    #pragma omp parallel
    {
        int worker_id = omp_get_thread_num();
        unsigned int seed = 123456789 ^ (worker_id * 9876543) ^ (unsigned int)time(NULL);
        int a_vec[L] = {0}, b_vec[L] = {0};
        SearchState states[L];
        for(int i=0; i<L; i++) states[i].cands = (Pair*)malloc(sizeof(Pair) * 200000);
        
        int step = 0, max_reached = 0;
        const int MAX_EVALS = 500;
        const int RESET_LIMIT = 50;
        int fail_count = 0;
        long loop_count = 0;
        long reset_count = 0;

        init_state(0, a_vec, b_vec, &states[0], &seed, valid_a_list, num_valid_a);
        
        while(step < L && !found) {
            loop_count++;
            
            // if (worker_id == 0 && loop_count % 100000 == 0) {
            //     // printf("\r[W00] 探索中... step:%d | 試行:%ld | Max:%d | リセット:%ld | Time:%.1fs", 
            //         //    step, loop_count, max_reached, reset_count, omp_get_wtime() - start_time);
            //     fflush(stdout);
            // }

            if (states[step].current_cand_idx < states[step].num_cands && states[step].current_cand_idx < MAX_EVALS) {
                Pair cand = states[step].cands[states[step].current_cand_idx++];
                int cur_idx = search_order[step];
                a_vec[cur_idx] = cand.a; 
                b_vec[cur_idx] = cand.b;
                
                if (check_cycles(step, a_vec, b_vec)) {
                    step++;
                    if (step > max_reached) {
                        max_reached = step;
                        // if(worker_id == 0) {
                        //     // printf("\r[W00] ★深層到達更新: step %d | リセット: %ld | Time: %.1fs          \n", 
                        //            max_reached, reset_count, omp_get_wtime()-start_time);
                        //     fflush(stdout);
                        // }
                    }
                    if (step < L) init_state(step, a_vec, b_vec, &states[step], &seed, valid_a_list, num_valid_a);
                }
            } else {
                fail_count++;
                if (fail_count >= RESET_LIMIT || step == 0) {
                    step = 0; 
                    fail_count = 0;
                    reset_count++;
                    init_state(0, a_vec, b_vec, &states[0], &seed, valid_a_list, num_valid_a);
                } else {
                    step--;
                }
            }
        }
        
        if (step == L) {
            #pragma omp critical
            { if(!found) { found = true; for(int i=0; i<L; i++){ sol_a[i]=a_vec[i]; sol_b[i]=b_vec[i]; } } }
        }
        for(int i=0; i<L; i++) free(states[i].cands);
    }
    
    for(int i=0; i<L; i++) {
        free((void*)active_cycles_x[i]);
        free((void*)active_cycles_z[i]);
    }

    if (found) {
        printf("\n\n【完全探索成功】合計時間: %.3f 秒\n", omp_get_wtime() - start_time);
        printf("a_vec: [");
        for(int i=0; i<L; i++) printf("%d%s", sol_a[i], i==L-1?"]\n":", ");
        printf("b_vec: [");
        for(int i=0; i<L; i++) printf("%d%s", sol_b[i], i==L-1?"]\n":", ");
    }
    return 0;
}