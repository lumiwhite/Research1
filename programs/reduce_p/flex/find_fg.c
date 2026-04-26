#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <string.h>

// コンパイルオプションで指定されなかった場合のデフォルト設定
#ifndef CYCLES_HEADER
#define CYCLES_HEADER "cycles_data_024.h"
#endif

// マクロを使用して動的にインクルード
#include CYCLES_HEADER

// ========================================================
// グローバル変数: アクティブ行と非可換ペア (コマンドラインから設定)
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

    for (int a_i = 0; a_i < num_valid_a; a_i++) {
        int a_val = valid_a_list[a_i];
        bool b_valid[P];
        for(int k=0; k<P; k++) b_valid[k] = true;

        for(int s=0; s<step; s++) {
            int prev_idx = search_order[s];
            
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

            for (int p_idx = 0; p_idx < NUM_NC_PAIRS; p_idx++) {
                if (i == NC_PAIRS[p_idx].i && j == NC_PAIRS[p_idx].j) {
                    commute_required = false; 
                    break;
                }
            }

            if (!commute_required) {
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
                    fprintf(stderr, "\n致命的なエラー: 直交性に必須なインデックス差分 (r=%d) を持つペア (%d, %d) が非可換に設定されている。\n", r, i, j);
                    exit(1);
                }
            }

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
                for(int k=0; k<P; k++) current_mask[k] = false;
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
    // argparse のような引数解析の実装
    int arg_idx = 1;
    bool save_txt = false; // ★追加: ファイル保存フラグ
    char *out_dir = "out"; // ★追加: 出力ディレクトリのデフォルト値

    while (arg_idx < argc) {
        if (strcmp(argv[arg_idx], "--p") == 0) {
            arg_idx++;
            if (arg_idx < argc) {
                P = atoi(argv[arg_idx]);
                arg_idx++;
            }
        } else if (strcmp(argv[arg_idx], "--active") == 0) {
            arg_idx++;
            while (arg_idx < argc && argv[arg_idx][0] != '-') {
                if (NUM_ACTIVE_ROWS < 6) {
                    ACTIVE_ROWS[NUM_ACTIVE_ROWS++] = atoi(argv[arg_idx]);
                }
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
        } else if (strcmp(argv[arg_idx], "--outdir") == 0) { // ★追加: フォルダ指定オプション
            arg_idx++;
            if (arg_idx < argc) {
                out_dir = argv[arg_idx];
                arg_idx++;
            }
        } else if (strcmp(argv[arg_idx], "--save") == 0) { // ★追加: 保存オプション
            save_txt = true;
            arg_idx++;
        } else {
            arg_idx++;
        }
    }
    // エラー処理とデフォルト値の設定
    if (P <= 0) {
        fprintf(stderr, "使用法: %s --p <Pの値> [--active <行1> <行2> ...] [--nc <i,j> <k,l> ...] [--save]\n", argv[0]);
        return 1;
    }
    if (NUM_ACTIVE_ROWS == 0) {
        ACTIVE_ROWS[0] = 0; ACTIVE_ROWS[1] = 2; ACTIVE_ROWS[2] = 4;
        NUM_ACTIVE_ROWS = 3;
    }

    int num_workers = omp_get_max_threads();
    printf("P = %d, %d 個の並列プロセスで完全全探索を開始...\n", P, num_workers);
    
    printf("【設定】アクティブ行: {");
    for(int i=0; i<NUM_ACTIVE_ROWS; i++) printf("%d%s", ACTIVE_ROWS[i], i==NUM_ACTIVE_ROWS-1?"}\n":", ");
    printf("【設定】非可換ペア: ");
    if (NUM_NC_PAIRS == 0) {
        printf("なし\n");
    } else {
        for(int i=0; i<NUM_NC_PAIRS; i++) printf("(%d, %d)%s", NC_PAIRS[i].i, NC_PAIRS[i].j, i==NUM_NC_PAIRS-1?"\n":", ");
    }
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

            if (states[step].current_cand_idx < states[step].num_cands && states[step].current_cand_idx < MAX_EVALS) {
                Pair cand = states[step].cands[states[step].current_cand_idx++];
                int cur_idx = search_order[step];
                a_vec[cur_idx] = cand.a; 
                b_vec[cur_idx] = cand.b;
                
                if (check_cycles(step, a_vec, b_vec)) {
                    step++;
                    if (step > max_reached) {
                        max_reached = step;
                        if(worker_id == 0) {
                            fflush(stdout);
                        }
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

        // ★追加: テキストファイルへの書き出し処理
        if (save_txt) {
            // ---------- ここから差し替え ----------
            // 出力用のフォルダを自動作成 (既に存在する場合は何もしない)
            char mkdir_cmd[1024];
            sprintf(mkdir_cmd, "mkdir -p %s", out_dir);
            int ret = system(mkdir_cmd);
            if (ret != 0) {
                fprintf(stderr, "警告: ディレクトリ '%s' の作成に失敗した可能性があります。\n", out_dir);
            }

            // 現在の日時を取得して文字列化
            time_t t = time(NULL);
            struct tm *tm_info = localtime(&t);
            char time_str[32];
            strftime(time_str, sizeof(time_str), "%Y%m%d_%H%M%S", tm_info);

            // アクティブ行の文字列化 (例: "024")
            char act_str[32] = "";
            for (int k = 0; k < NUM_ACTIVE_ROWS; k++) {
                char tmp[16];
                sprintf(tmp, "%d", ACTIVE_ROWS[k]);
                strcat(act_str, tmp);
            }

            // 非可換ペアの文字列化 (例: "0-1_2-3")
            char nc_str[256] = "";
            if (NUM_NC_PAIRS == 0) {
                strcpy(nc_str, "none");
            } else {
                for (int k = 0; k < NUM_NC_PAIRS; k++) {
                    char tmp[32];
                    // 最初の要素以外はアンダースコアで繋ぐ
                    sprintf(tmp, "%s%d-%d", (k == 0) ? "" : "_", NC_PAIRS[k].i, NC_PAIRS[k].j);
                    strcat(nc_str, tmp);
                }
            }

            char filename[2048]; // パスが長くなる可能性を考慮して少し拡張
            // 指定されたフォルダ内にファイルを生成
            sprintf(filename, "%s/P_%d_act_%s_nc_%s_%s.txt", out_dir, P, act_str, nc_str, time_str); 
            // ---------- ここまで差し替え ----------
            FILE *fp = fopen(filename, "w");
            if (fp) {
                // 1行目: P J L
                fprintf(fp, "%d %d %d\n", P, NUM_ACTIVE_ROWS, L);
                
                // 2行目: Seed等の代用ID (とりあえず a_vec の先頭を使用)
                fprintf(fp, "%d\n", sol_a[0]);
                
                // 3行目: a_vecの前半
                for (int i = 0; i < L_H; i++) fprintf(fp, "%d%s", sol_a[i], (i == L_H - 1) ? "" : " ");
                fprintf(fp, "\n");
                
                // 4行目: b_vecの前半
                for (int i = 0; i < L_H; i++) fprintf(fp, "%d%s", sol_b[i], (i == L_H - 1) ? "" : " ");
                fprintf(fp, "\n");
                
                // 5行目: a_vecの後半
                for (int i = L_H; i < L; i++) fprintf(fp, "%d%s", sol_a[i], (i == L - 1) ? "" : " ");
                fprintf(fp, "\n");
                
                // 6行目: b_vecの後半
                for (int i = L_H; i < L; i++) fprintf(fp, "%d%s", sol_b[i], (i == L - 1) ? "" : " ");
                fprintf(fp, "\n");
                
                // 7行目: active_set
                fprintf(fp, "active_set=");
                for (int i = 0; i < NUM_ACTIVE_ROWS; i++) {
                    fprintf(fp, "%d%s", ACTIVE_ROWS[i], (i == NUM_ACTIVE_ROWS - 1) ? "" : ",");
                }
                fprintf(fp, "\n");

                fclose(fp);
                printf("シミュレータ用パラメータを '%s' に保存しました。\n", filename);
            } else {
                fprintf(stderr, "エラー: '%s' への書き込みに失敗しました。\n", filename);
            }
        }
    }
    return 0;
}