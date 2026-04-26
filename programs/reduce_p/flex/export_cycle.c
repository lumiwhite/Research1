#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_CYCLES 30000
#define L 12

// ==========================================
// データ構造
// ==========================================
typedef struct {
    int len;
    int nodes[6];
} Cycle;

Cycle cx[MAX_CYCLES];
Cycle cz[MAX_CYCLES];
int total_cycles = 0;

typedef struct {
    Cycle *cycles[MAX_CYCLES];
    int count;
} ByIdx;

ByIdx by_idx_x[L];
ByIdx by_idx_z[L];

// ==========================================
// 補助関数群
// ==========================================
int func_X(int c, int v) {
    if (v < 6) {
        int res = (v - c) % 6;
        return res < 0 ? res + 6 : res;
    } else {
        int res = (v - 6 - c) % 6;
        return 6 + (res < 0 ? res + 6 : res);
    }
}

int func_Z(int c, int v) {
    if (v < 6) {
        int res = (c - v) % 6;
        return 6 + (res < 0 ? res + 6 : res);
    } else {
        int res = (c - (v - 6)) % 6;
        return res < 0 ? res + 6 : res;
    }
}

int get_max(Cycle *cyc) {
    int max_val = -1;
    for (int i = 0; i < cyc->len; i++) {
        if (cyc->nodes[i] > max_val) {
            max_val = cyc->nodes[i];
        }
    }
    return max_val;
}

// ==========================================
// サイクル生成ロジック
// ==========================================
void add_cycle(int *c_arr, int *v_arr, int k) {
    Cycle *x_cyc = &cx[total_cycles];
    Cycle *z_cyc = &cz[total_cycles];
    
    x_cyc->len = k * 2;
    z_cyc->len = k * 2;
    
    int idx = 0;
    for (int i = 0; i < k; i++) {
        int c_current = c_arr[i];
        int v_current = v_arr[i];
        int c_next = c_arr[(i + 1) % k];
        
        x_cyc->nodes[idx]     = func_X(c_current, v_current);
        x_cyc->nodes[idx + 1] = func_X(c_next, v_current);
        
        z_cyc->nodes[idx]     = func_Z(c_current, v_current);
        z_cyc->nodes[idx + 1] = func_Z(c_next, v_current);
        
        idx += 2;
    }
    total_cycles++;
}

void get_cycles(int *active, int num_active) {
    // 長さ4のサイクル抽出
    for (int i = 0; i < num_active - 1; i++) {
        for (int j = i + 1; j < num_active; j++) {
            for (int v0 = 0; v0 < 11; v0++) {
                for (int v1 = v0 + 1; v1 < 12; v1++) {
                    int c_arr[2] = {active[i], active[j]};
                    int v_arr[2] = {v0, v1};
                    add_cycle(c_arr, v_arr, 2);
                }
            }
        }
    }
    
    // 長さ6のサイクル抽出
    for (int v0 = 0; v0 < 10; v0++) {
        for (int v1 = v0 + 1; v1 < 11; v1++) {
            for (int v2 = v1 + 1; v2 < 12; v2++) {
                for (int i = 0; i < num_active - 2; i++) {
                    for (int j = i + 1; j < num_active - 1; j++) {
                        for (int k = j + 1; k < num_active; k++) {
                            int c_arr[3] = {active[i], active[j], active[k]};
                            int v_pats[6][3] = {
                                {v0, v1, v2}, {v0, v2, v1},
                                {v1, v0, v2}, {v1, v2, v0},
                                {v2, v0, v1}, {v2, v1, v0}
                            };
                            for (int p = 0; p < 6; p++) {
                                add_cycle(c_arr, v_pats[p], 3);
                            }
                        }
                    }
                }
            }
        }
    }
}

// ==========================================
// ファイル出力
// ==========================================
void write_cycles_header(const char *filename) {
    char guard_name[256];
    int idx = 0;
    for (int i = 0; filename[i] != '\0'; i++) {
        char ch = filename[i];
        if (isalpha(ch)) ch = toupper(ch);
        else if (ch == '.') ch = '_';
        guard_name[idx++] = ch;
    }
    guard_name[idx] = '\0';
    
    FILE *f = fopen(filename, "w");
    if (!f) {
        fprintf(stderr, "エラー: ファイル '%s' を開くことができない。\n", filename);
        exit(1);
    }
    
    fprintf(f, "#ifndef %s\n#define %s\n\n", guard_name, guard_name);
    fprintf(f, "typedef struct {\n    int len;\n    int nodes[6];\n} Cycle;\n\n");
    
    // X, Zの2回ループ
    for (int t = 0; t < 2; t++) {
        char name = t == 0 ? 'x' : 'z';
        ByIdx *by_idx = t == 0 ? by_idx_x : by_idx_z;
        
        int max_c = 1;
        for (int i = 0; i < L; i++) {
            if (by_idx[i].count > max_c) {
                max_c = by_idx[i].count;
            }
        }
        
        fprintf(f, "const int num_cycles_%c_by_idx[%d] = {", name, L);
        for (int i = 0; i < L; i++) {
            fprintf(f, "%d%s", by_idx[i].count, (i == L - 1) ? "" : ", ");
        }
        fprintf(f, "};\n\n");
        
        fprintf(f, "const Cycle cycles_%c_by_idx[%d][%d] = {\n", name, L, max_c);
        for (int i = 0; i < L; i++) {
            fprintf(f, "    {\n");
            if (by_idx[i].count == 0) {
                fprintf(f, "        {0, {0, 0, 0, 0, 0, 0}},\n");
            }
            for (int c = 0; c < by_idx[i].count; c++) {
                Cycle *cyc = by_idx[i].cycles[c];
                int padded[6] = {0};
                for (int k = 0; k < cyc->len; k++) padded[k] = cyc->nodes[k];
                
                fprintf(f, "        {%d, {%d, %d, %d, %d, %d, %d}},\n",
                        cyc->len, padded[0], padded[1], padded[2], padded[3], padded[4], padded[5]);
            }
            fprintf(f, "    },\n");
        }
        fprintf(f, "};\n\n");
    }
    
    fprintf(f, "#endif\n");
    fclose(f);
}

// ==========================================
// メインルーチン
// ==========================================
int main(int argc, char *argv[]) {
    int active[12];
    int num_active = 0;
    char out_file[256] = "";
    
    // 1. コマンドライン引数のパース
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--active") == 0) {
            i++;
            while (i < argc && argv[i][0] != '-') {
                if (num_active < 12) {
                    active[num_active++] = atoi(argv[i]);
                }
                i++;
            }
            i--; // --active の読み取りループ調整
        } else if (strcmp(argv[i], "--out") == 0) {
            i++;
            if (i < argc) {
                strncpy(out_file, argv[i], sizeof(out_file) - 1);
            }
        }
    }
    
    // デフォルト値の設定
    if (num_active == 0) {
        active[0] = 0; active[1] = 2; active[2] = 4;
        num_active = 3;
    }
    
    // ファイル名の自動決定
    if (out_file[0] == '\0') {
        char active_str[128] = "";
        for (int i = 0; i < num_active; i++) {
            char buf[16];
            sprintf(buf, "%d", active[i]);
            strcat(active_str, buf);
        }
        sprintf(out_file, "cycles_data_%s.h", active_str);
    }
    
    printf("アクティブ行 [");
    for (int i = 0; i < num_active; i++) printf("%d%s", active[i], i == num_active - 1 ? "" : ", ");
    printf("] のサイクルデータを計算中...\n");
    
    // 2. サイクルデータの生成
    get_cycles(active, num_active);
    printf("総サイクル数: %d\n", total_cycles);
    
    // 3. 最大インデックスによる分類
    for (int i = 0; i < L; i++) {
        by_idx_x[i].count = 0;
        by_idx_z[i].count = 0;
    }
    
    for (int c = 0; c < total_cycles; c++) {
        int max_x = get_max(&cx[c]);
        by_idx_x[max_x].cycles[by_idx_x[max_x].count++] = &cx[c];
        
        int max_z = get_max(&cz[c]);
        by_idx_z[max_z].cycles[by_idx_z[max_z].count++] = &cz[c];
    }
    
    // 4. ヘッダーファイルのエクスポート
    write_cycles_header(out_file);
    printf("--> %s のエクスポートが完了した。\n", out_file);
    
    return 0;
}