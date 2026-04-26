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

/* * 改善点:
 * 1. 論文 Section 6.2 の可換表に基づき、r=3 における (F0,G3) と (F1,G2) のみを非可換に設定。
 * 2. a_i, c_j が決まっている場合、可換条件は b_i, d_j に関する線形合同式となるため、
 * idx >= L_H の候補生成時に条件を満たす b のみを抽出して探索空間を劇的に削減。
 * 3. 不要な is_commute チェックを減らし、代数的な直接判定に置き換え。
 */

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
    // APM の可換条件: (a1 - 1) * b2 == (a2 - 1) * b1 (mod P)
    return mod((a1 - 1) * b2, P) == mod((a2 - 1) * b1, P);
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
        // 前半 Fi の探索は自由
        for(int i = 0; i < num_valid_a; i++) {
            int a_val = valid_a_list[i];
            for(int b = 0; b < P; b++) {
                state->cands[state->num_cands].a = a_val;
                state->cands[state->num_cands].b = b;
                state->num_cands++;
            }
        }
    } else {
        // 後半 Gj の探索。Fi との可換条件を考慮
        int j = idx - L_H;
        for(int a_i = 0; a_i < num_valid_a; a_i++) {
            int c_j = valid_a_list[a_i];
            for(int d_j = 0; d_j < P; d_j++) {
                bool valid = true;
                for(int i = 0; i < L_H; i++) {
                    int r = (i + j) % L_H;
                    bool commute_required;
                    
                    // 論文可換表に基づく条件設定 [cite: 234]
                    if (r == 3) {
                        // (F0, G3) と (F1, G2) のみが非可換 [cite: 232, 241]
                        if ((i == 0 && j == 3) || (i == 1 && j == 2)) commute_required = false;
                        else commute_required = true;
                    } else {
                        commute_required = true; // r in {0,1,2,4,5} [cite: 212, 230]
                    }

                    if (is_commute(a_vec[i], b_vec[i], c_j, d_j) != commute_required) {
                        valid = false;
                        break;
                    }
                }
                if (valid) {
                    state->cands[state->num_cands].a = c_j;
                    state->cands[state->num_cands].b = d_j;
                    state->num_cands++;
                }
            }
        }
    }
    shuffle_pairs(state->cands, state->num_cands, seed);
}

int main() {
    int num_workers = omp_get_max_threads();
    printf("%d 個の並列プロセスで探索を開始 (可換条件最適化済み)...\n", num_workers);
    
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
        unsigned int seed = 123456789 ^ (worker_id * 9876543) ^ (unsigned int)time(NULL);
        
        int a_vec[L] = {0}, b_vec[L] = {0};
        SearchState states[L];
        for(int i=0; i<L; i++) states[i].cands = (Pair*)malloc(sizeof(Pair) * 200000);
        
        int idx = 0, max_reached = 0, reset_count = 0;
        const int MAX_EVALS_PER_STATE = 5000; // 絞り込みが効くため上限を拡大
        const int STALL_LIMIT = 1000;
        int fail_count = 0;
        
        init_state(0, a_vec, b_vec, &states[0], &seed, valid_a_list, num_valid_a);
        
        while(idx < L && !found) {
            if (states[idx].current_cand_idx < states[idx].num_cands && states[idx].current_cand_idx < MAX_EVALS_PER_STATE) {
                Pair cand = states[idx].cands[states[idx].current_cand_idx++];
                a_vec[idx] = cand.a;
                b_vec[idx] = cand.b;
                
                if (check_cycles(idx, a_vec, b_vec)) {
                    idx++;
                    if (idx > max_reached) {
                        max_reached = idx;
                        // 進捗表示
                        if (worker_id == 0) {
                            printf("\r[W00] 最深到達点: %d | 経過時間: %.1fs", max_reached, omp_get_wtime() - start_time);
                            fflush(stdout);
                        }
                    }
                    if (idx < L) init_state(idx, a_vec, b_vec, &states[idx], &seed, valid_a_list, num_valid_a);
                }
            } else {
                fail_count++;
                if (fail_count >= STALL_LIMIT || idx == 0) {
                    idx = 0; reset_count++; fail_count = 0;
                    init_state(0, a_vec, b_vec, &states[0], &seed, valid_a_list, num_valid_a);
                } else {
                    idx--;
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