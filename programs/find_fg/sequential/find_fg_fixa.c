#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include "cycles_data.h"

#define P 768
#define L 12
#define L_H 6

/**
 * 改良および修正点:
 * 1. 論文Table 1に基づく線形部(a)の固定。
 * 2. 探索空間の削減: 並進部(b)のみを探索。
 * 3. 沼対策: MAX_EVALSとRESET_LIMITを極端に下げ、ランダムリスタートを促進。
 * 4. リセット回数の表示: bを一から作り直した回数(reset_count)をログに追加。
 */

const int FIXED_A[L] = {
    763, 679, 397, 61, 697, 373,
    289, 257, 625, 41, 193, 449
};

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

bool check_cycles(int idx, const int *a_vec, const int *b_vec) {
    for(int c = 0; c < num_cycles_x_by_idx[idx]; c++) {
        int len = cycles_x_by_idx[idx][c].len;
        const int *nodes = cycles_x_by_idx[idx][c].nodes;
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
    for(int c = 0; c < num_cycles_z_by_idx[idx]; c++) {
        int len = cycles_z_by_idx[idx][c].len;
        const int *nodes = cycles_z_by_idx[idx][c].nodes;
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

void init_state(int idx, const int *a_vec, const int *b_vec, SearchState *state, unsigned int *seed) {
    state->num_cands = 0;
    state->current_cand_idx = 0;
    
    int a_val = FIXED_A[idx];

    if (idx < L_H) {
        for(int b = 0; b < P; b++) {
            state->cands[state->num_cands].a = a_val;
            state->cands[state->num_cands].b = b;
            state->num_cands++;
        }
    } else {
        int j = idx - L_H;
        int sol_buffer[P];
        
        bool d_valid[P];
        for(int k=0; k<P; k++) d_valid[k] = true;

        for(int i = 0; i < L_H; i++) {
            int r = (i + j) % L_H;
            int A = a_vec[i] - 1;
            int B = (a_val - 1) * b_vec[i];
            int num_sol = solve_linear_congruence(A, B, P, sol_buffer);

            bool commute_required = true;
            if (r == 3 && ((i == 0 && j == 3) || (i == 1 && j == 2))) {
                commute_required = false; 
            }

            if (commute_required) {
                bool current_mask[P] = {false};
                for(int s=0; s<num_sol; s++) current_mask[sol_buffer[s]] = true;
                for(int k=0; k<P; k++) d_valid[k] &= current_mask[k];
            } else {
                for(int s=0; s<num_sol; s++) d_valid[sol_buffer[s]] = false;
            }
        }

        for(int d=0; d<P; d++) {
            if(d_valid[d]) {
                state->cands[state->num_cands].a = a_val;
                state->cands[state->num_cands].b = d;
                state->num_cands++;
            }
        }
    }
    shuffle_pairs(state->cands, state->num_cands, seed);
}

int main() {
    int num_workers = omp_get_max_threads();
    printf("%d 個の並列プロセスで探索を開始 (線形部固定、超低粘着化)...\n", num_workers);
    fflush(stdout);
    
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
        
        int idx = 0, max_reached = 0;
        const int MAX_EVALS = 100;     // 局所解を避けるため極端に小さく設定
        const int RESET_LIMIT = 50;    // 局所解を避けるため極端に小さく設定
        int fail_count = 0;
        long loop_count = 0;
        long reset_count = 0;         // リセット回数を記録する変数を追加

        init_state(0, a_vec, b_vec, &states[0], &seed);
        
        while(idx < L && !found) {
            loop_count++;
            
            if (worker_id == 0 && loop_count % 100000 == 0) {
                // リセット回数を出力に追加
                printf("\r[W00] 探索中... idx:%d | 試行:%ld | Max:%d | リセット:%ld | Time:%.1fs", 
                       idx, loop_count, max_reached, reset_count, omp_get_wtime() - start_time);
                fflush(stdout);
            }

            if (states[idx].current_cand_idx < states[idx].num_cands && states[idx].current_cand_idx < MAX_EVALS) {
                Pair cand = states[idx].cands[states[idx].current_cand_idx++];
                a_vec[idx] = cand.a; b_vec[idx] = cand.b;
                
                if (check_cycles(idx, a_vec, b_vec)) {
                    idx++;
                    if (idx > max_reached) {
                        max_reached = idx;
                        if(worker_id == 0) {
                            printf("\r[W00] ★深層到達更新: %d | リセット: %ld | Time: %.1fs          \n", 
                                   max_reached, reset_count, omp_get_wtime()-start_time);
                            fflush(stdout);
                        }
                    }
                    if (idx < L) init_state(idx, a_vec, b_vec, &states[idx], &seed);
                }
            } else {
                fail_count++;
                if (fail_count >= RESET_LIMIT || idx == 0) {
                    idx = 0; 
                    fail_count = 0;
                    reset_count++; // ここでリセット回数をカウントアップ
                    init_state(0, a_vec, b_vec, &states[0], &seed);
                } else {
                    idx--;
                }
            }
        }
        
        if (idx == L) {
            #pragma omp critical
            { if(!found) { found = true; for(int i=0; i<L; i++){ sol_a[i]=a_vec[i]; sol_b[i]=b_vec[i]; } } }
        }
        for(int i=0; i<L; i++) free(states[i].cands);
    }
    
    if (found) {
        printf("\n\n【探索成功】合計時間: %.3f 秒\n", omp_get_wtime() - start_time);
        printf("a_vec: [");
        for(int i=0; i<L; i++) printf("%d%s", sol_a[i], i==L-1?"]\n":", ");
        printf("b_vec: [");
        for(int i=0; i<L; i++) printf("%d%s", sol_b[i], i==L-1?"]\n":", ");
    }
    return 0;
}