#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "cycles_data012.h"

#define P 768
#define L 12
#define L_H 6

// --- 補助関数群 ---
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

// --- 検証関数群 ---

// 条件B: F_i と G_j の可換・非可換条件を判定 (直交性条件)
bool check_condition_AB(const int *a_vec, const int *b_vec) {
    for (int i = 0; i < L_H; i++) {
        for (int j = 0; j < L_H; j++) {
            int r = (i + j) % L_H;
            bool should_commute = true;
            
            // r=3 の特定のペア(F0とG3、F1とG2)のみ非可換
            if (r == 3 && ((i == 0 && j == 3) || (i == 1 && j == 2))) {
                should_commute = false;
            }
            
            int a_F = a_vec[i];
            int b_F = b_vec[i];
            int a_G = a_vec[L_H + j];
            int b_G = b_vec[L_H + j];
            
            if (is_commute(a_F, b_F, a_G, b_G) != should_commute) {
                return false;
            }
        }
    }
    return true;
}

// 条件C: サイクル回避条件を判定
bool check_condition_C(const int *a_vec, const int *b_vec) {
    // 探索用の動的マッピング(active_cycles)ではなく、
    // 静的な cycles_data.h の構造を用いて全インデックスを独立してチェックする
    for (int idx = 0; idx < L; idx++) {
        
        // X側のサイクルチェック
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
            if(is_closed(c_a, c_b)) return false;
        }
        
        // Z側のサイクルチェック
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
            if(is_closed(c_a, c_b)) return false;
        }
    }
    return true;
}

// 総合判定関数 (外部からはこの関数を呼び出す)
bool verify_all_conditions(const int *a_vec, const int *b_vec) {
    if (!check_condition_AB(a_vec, b_vec)) {
        printf("検証失敗: 条件B (可換・非可換条件) を満たしていません。\n");
        return false;
    }
    if (!check_condition_C(a_vec, b_vec)) {
        printf("検証失敗: 条件C (サイクル回避条件) を満たしていません。\n");
        return false;
    }
    
    printf("検証成功: 全ての条件を満たしています。\n");
    return true;
}

// --- 実行部 ---
int main() {
    // 検証対象のパラメータをここに設定する
    // 今回は先ほど高速探索で発見された「新しい解」をテストデータとして使用
    const int test_a_vec[L] = {569, 233, 1, 449, 721, 481, 517, 421, 331, 295, 541, 697};
    const int test_b_vec[L] = {396, 4, 640, 224, 8, 48, 426, 90, 57, 207, 582, 204};

    printf("========================================\n");
    printf("パラメータ検証テストを開始します。\n");
    printf("a_vec: [");
    for(int i=0; i<L; i++) printf("%d%s", test_a_vec[i], i==L-1?"]\n":", ");
    printf("b_vec: [");
    for(int i=0; i<L; i++) printf("%d%s", test_b_vec[i], i==L-1?"]\n":", ");
    printf("----------------------------------------\n");

    verify_all_conditions(test_a_vec, test_b_vec);
    
    printf("========================================\n");

    return 0;
}