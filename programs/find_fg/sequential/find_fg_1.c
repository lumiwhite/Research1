// デフォルトの実装

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

// aとbのペアを格納する構造体
typedef struct {
    int a;
    int b;
} Pair;

// 各インデックスの探索状態
typedef struct {
    Pair *cands;
    int num_cands;
    int current_cand_idx;
} SearchState;

// 負の剰余をPythonと同じ挙動にする関数
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

inline bool is_commute(int a1, int b1, int a2, int b2) {
    return mod(a1 * b2 + b1, P) == mod(a2 * b1 + b2, P);
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
            if(i % 2 == 0) {
                composite(c_a, c_b, f_a, f_b, &c_a, &c_b);
            } else {
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
            if(i % 2 == 1) {
                composite(c_a, c_b, f_a, f_b, &c_a, &c_b);
            } else {
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

void init_state(int idx, const int *a_vec, const int *b_vec, SearchState *state, unsigned int *seed, const int *valid_a_list, int num_valid_a) {
    state->num_cands = 0;
    state->current_cand_idx = 0;
    
    if (idx < L_H) {
        for(int i = 0; i < num_valid_a; i++) {
            int a_val = valid_a_list[i];
            for(int b = 0; b < P; b++) {
                state->cands[state->num_cands].a = a_val;
                state->cands[state->num_cands].b = b;
                state->num_cands++;
            }
        }
    } else {
        int j = idx - L_H;
        for(int a_i = 0; a_i < num_valid_a; a_i++) {
            int a_val = valid_a_list[a_i];
            for(int b = 0; b < P; b++) {
                bool valid = true;
                for(int i = 0; i < L_H; i++) {
                    int r = (i + j) % L_H;
                    if (r != 3) {
                        if (!is_commute(a_val, b, a_vec[i], b_vec[i])) { valid = false; break; }
                    } else {
                        if ((i == 0 && j == 3) || (i == 1 && j == 2)) {
                            if (is_commute(a_val, b, a_vec[i], b_vec[i])) { valid = false; break; }
                        }
                    }
                }
                if (valid) {
                    state->cands[state->num_cands].a = a_val;
                    state->cands[state->num_cands].b = b;
                    state->num_cands++;
                }
            }
        }
    }
    shuffle_pairs(state->cands, state->num_cands, seed);
}

int main() {
    int num_workers = omp_get_max_threads();
    printf("%d 個のプロセスでランダムA並列探索を開始する...\n", num_workers);
    for(int i=0; i<num_workers; i++) printf("\n");
    
    // aの有効な候補を事前計算 (Pythonの VALID_A_LIST 相当)
    int valid_a_list[P];
    int num_valid_a = 0;
    for(int a = 1; a < P; a++) {
        if(gcd(a, P) == 1) valid_a_list[num_valid_a++] = a;
    }
    
    volatile bool found = false;
    int sol_a[L], sol_b[L];
    double start_time = omp_get_wtime();
    
    #pragma omp parallel
    {
        int worker_id = omp_get_thread_num();
        unsigned int seed = 123456789 ^ (worker_id * 9876543) ^ time(NULL);
        
        int a_vec[L] = {0}, b_vec[L] = {0};
        SearchState states[L];
        
        // 最大19万個の候補配列をヒープに割り当て（スタックオーバーフロー防止）
        for(int i=0; i<L; i++) {
            states[i].cands = (Pair*)malloc(sizeof(Pair) * 200000);
        }
        
        int idx = 0, max_reached = 0, reset_count = 0;
        int active_fail_count = 0, latent_fail_count = 0;
        
        // ランダム探索用のパラメータ調整
        const int MAX_EVALS_PER_STATE = 2000; 
        const int ACTIVE_STALL_LIMIT = 500;
        const int LATENT_STALL_LIMIT = 1000;
        
        init_state(0, a_vec, b_vec, &states[0], &seed, valid_a_list, num_valid_a);
        int loop_counter = 0;
        
        while(idx < L && !found) {
            loop_counter++;
            if (loop_counter % 200000 == 0) {
                #pragma omp critical
                {
                    printf("\033[%dA\r\033[K[W%02d] idx:%2d | Max:%2d | Rst:%3d | %5.1fs\n\033[%dB",
                           num_workers - worker_id, worker_id, idx, max_reached, reset_count, 
                           omp_get_wtime() - start_time, num_workers - worker_id - 1);
                    fflush(stdout);
                }
            }
            
            // 一定数の候補を試したら、見込みなしと判断してバックトラックする
            if (states[idx].current_cand_idx < states[idx].num_cands && states[idx].current_cand_idx < MAX_EVALS_PER_STATE) {
                Pair cand = states[idx].cands[states[idx].current_cand_idx++];
                a_vec[idx] = cand.a;
                b_vec[idx] = cand.b;
                
                if (check_cycles(idx, a_vec, b_vec)) {
                    idx++;
                    if (idx > max_reached) max_reached = idx;
                    if (idx < L) init_state(idx, a_vec, b_vec, &states[idx], &seed, valid_a_list, num_valid_a);
                }
            } else {
                if (idx >= L_H) {
                    latent_fail_count++;
                    if (latent_fail_count >= LATENT_STALL_LIMIT) {
                        idx = 0; max_reached = 0; reset_count++;
                        latent_fail_count = 0; active_fail_count = 0;
                        init_state(0, a_vec, b_vec, &states[0], &seed, valid_a_list, num_valid_a);
                    } else idx--;
                } else {
                    active_fail_count++;
                    if (active_fail_count >= ACTIVE_STALL_LIMIT) {
                        idx = 0; max_reached = 0; reset_count++;
                        active_fail_count = 0; latent_fail_count = 0;
                        init_state(0, a_vec, b_vec, &states[0], &seed, valid_a_list, num_valid_a);
                    } else idx = (idx > 0) ? idx - 1 : 0;
                }
            }
        }
        
        if (idx == L) {
            #pragma omp critical
            {
                if (!found) {
                    found = true;
                    for(int i=0; i<L; i++) { sol_a[i] = a_vec[i]; sol_b[i] = b_vec[i]; }
                }
            }
        }
        
        for(int i=0; i<L; i++) free(states[i].cands);
    }
    
    if (found) {
        printf("\n\n【探索成功！】計算時間: %.3f 秒\n", omp_get_wtime() - start_time);
        printf("a_vec: ["); for(int i=0; i<L; i++) printf("%d%s", sol_a[i], i==L-1?"":", "); printf("]\n");
        printf("b_vec: ["); for(int i=0; i<L; i++) printf("%d%s", sol_b[i], i==L-1?"":", "); printf("]\n");
    }
    return 0;
}