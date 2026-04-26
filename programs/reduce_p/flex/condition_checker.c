#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

#ifndef CYCLES_HEADER
#define CYCLES_HEADER "cycles_data_024.h"
#endif

#include CYCLES_HEADER

// ========================================================
// グローバル変数 (コマンドラインから設定)
// ========================================================
int ACTIVE_ROWS[6];
int NUM_ACTIVE_ROWS = 0;

typedef struct {
    int i;
    int j;
} NCPair;

NCPair NC_PAIRS[36];
int NUM_NC_PAIRS = 0;

int P = 0;
#define L 12
#define L_H 6

// ========================================================
// 補助関数群
// ========================================================
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
    return mod((a1 - 1) * b2, P) == mod((a2 - 1) * b1, P);
}

// ========================================================
// 検証関数群
// ========================================================
bool check_condition_AB(const int *a_vec, const int *b_vec) {
    for (int i = 0; i < L_H; i++) {
        for (int j = 0; j < L_H; j++) {
            int r = (i + j) % L_H;
            bool should_commute = true;
            
            for (int p_idx = 0; p_idx < NUM_NC_PAIRS; p_idx++) {
                if (i == NC_PAIRS[p_idx].i && j == NC_PAIRS[p_idx].j) {
                    should_commute = false;
                    break;
                }
            }

            if (!should_commute) {
                bool is_orthogonal_r = false;
                for (int k = 0; k < NUM_ACTIVE_ROWS; k++) {
                    for (int m = 0; m < NUM_ACTIVE_ROWS; m++) {
                        int diff = (ACTIVE_ROWS[k] - ACTIVE_ROWS[m] + L_H) % L_H;
                        if (r == diff) {
                            is_orthogonal_r = true;
                            break;
                        }
                    }
                    if (is_orthogonal_r) break;
                }
                
                if (is_orthogonal_r) {
                    printf("設定エラー: 直交性に必須な差分 (r=%d) を持つペア (%d, %d) が非可換に設定されています。\n", r, i, j);
                    return false;
                }
            }
            
            int a_F = a_vec[i];
            int b_F = b_vec[i];
            int a_G = a_vec[L_H + j];
            int b_G = b_vec[L_H + j];
            
            if (is_commute(a_F, b_F, a_G, b_G) != should_commute) {
                printf("検証失敗: ペア (%d, %d) [r=%d] の可換性が期待値 (%s) と一致しません。\n", 
                       i, j, r, should_commute ? "可換" : "非可換");
                return false;
            }
        }
    }
    return true;
}

bool check_condition_C(const int *a_vec, const int *b_vec) {
    for (int idx = 0; idx < L; idx++) {
        for(int c = 0; c < num_cycles_x_by_idx[idx]; c++) {
            int len = cycles_x_by_idx[idx][c].len;
            const int *nodes = cycles_x_by_idx[idx][c].nodes;
            int c_a = 1, c_b = 0;
            for(int i = 0; i < len; i++) {
                int n = nodes[i];
                if(i % 2 == 0) {
                    composite(c_a, c_b, a_vec[n], b_vec[n], &c_a, &c_b);
                } else {
                    int inv_a, inv_b;
                    func_inv(a_vec[n], b_vec[n], &inv_a, &inv_b);
                    composite(c_a, c_b, inv_a, inv_b, &c_a, &c_b);
                }
            }
            if(is_closed(c_a, c_b)) {
                printf("検証失敗: X側サイクル (長さ%d) が閉じています。\n", len);
                return false;
            }
        }
        
        for(int c = 0; c < num_cycles_z_by_idx[idx]; c++) {
            int len = cycles_z_by_idx[idx][c].len;
            const int *nodes = cycles_z_by_idx[idx][c].nodes;
            int c_a = 1, c_b = 0;
            for(int i = 0; i < len; i++) {
                int n = nodes[i];
                if(i % 2 == 1) {
                    composite(c_a, c_b, a_vec[n], b_vec[n], &c_a, &c_b);
                } else {
                    int inv_a, inv_b;
                    func_inv(a_vec[n], b_vec[n], &inv_a, &inv_b);
                    composite(c_a, c_b, inv_a, inv_b, &c_a, &c_b);
                }
            }
            if(is_closed(c_a, c_b)) {
                printf("検証失敗: Z側サイクル (長さ%d) が閉じています。\n", len);
                return false;
            }
        }
    }
    return true;
}

bool verify_all_conditions(const int *a_vec, const int *b_vec) {
    if (!check_condition_AB(a_vec, b_vec)) {
        printf("結果: ❌ 条件B (可換・非可換条件) を満たしていません。\n");
        return false;
    }
    if (!check_condition_C(a_vec, b_vec)) {
        printf("結果: ❌ 条件C (サイクル回避条件) を満たしていません。\n");
        return false;
    }
    
    printf("結果: ✅ 全ての条件を満たしています。\n");
    return true;
}

// ========================================================
// メイン部
// ========================================================
int main(int argc, char *argv[]) {
    int test_a_vec[L] = {0};
    int test_b_vec[L] = {0};
    int a_count = 0, b_count = 0;

    int arg_idx = 1;
    while (arg_idx < argc) {
        if (strcmp(argv[arg_idx], "--p") == 0) {
            arg_idx++;
            if (arg_idx < argc) { P = atoi(argv[arg_idx]); arg_idx++; }
        } else if (strcmp(argv[arg_idx], "--active") == 0) {
            arg_idx++;
            while (arg_idx < argc && argv[arg_idx][0] != '-') {
                if (NUM_ACTIVE_ROWS < 6) ACTIVE_ROWS[NUM_ACTIVE_ROWS++] = atoi(argv[arg_idx]);
                arg_idx++;
            }
        } else if (strcmp(argv[arg_idx], "--nc") == 0) {
            arg_idx++;
            while (arg_idx < argc && argv[arg_idx][0] != '-') {
                int r, c;
                if (sscanf(argv[arg_idx], "%d,%d", &r, &c) == 2) {
                    if (NUM_NC_PAIRS < 36) {
                        NC_PAIRS[NUM_NC_PAIRS].i = r;
                        NC_PAIRS[NUM_NC_PAIRS].j = c;
                        NUM_NC_PAIRS++;
                    }
                }
                arg_idx++;
            }
        } else if (strcmp(argv[arg_idx], "--a") == 0) {
            arg_idx++;
            while (arg_idx < argc && argv[arg_idx][0] != '-' && a_count < L) {
                test_a_vec[a_count++] = atoi(argv[arg_idx]);
                arg_idx++;
            }
        } else if (strcmp(argv[arg_idx], "--b") == 0) {
            arg_idx++;
            while (arg_idx < argc && argv[arg_idx][0] != '-' && b_count < L) {
                test_b_vec[b_count++] = atoi(argv[arg_idx]);
                arg_idx++;
            }
        } else {
            arg_idx++;
        }
    }

    if (P <= 0 || a_count != L || b_count != L) {
        fprintf(stderr, "使用法: %s --p <P> --active <rows> [--nc <pairs>] --a <12 vals> --b <12 vals>\n", argv[0]);
        return 1;
    }

    printf("========================================\n");
    printf("パラメータ検証テストを開始 (P = %d)\n", P);
    
    printf("\na_vec: [");
    for(int i=0; i<L; i++) printf("%d%s", test_a_vec[i], i==L-1?"]\n":", ");
    printf("b_vec: [");
    for(int i=0; i<L; i++) printf("%d%s", test_b_vec[i], i==L-1?"]\n":", ");
    printf("----------------------------------------\n");

    verify_all_conditions(test_a_vec, test_b_vec);
    printf("========================================\n");

    return 0;
}