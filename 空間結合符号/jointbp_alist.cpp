#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

struct EdgeRef {
    int check = -1;
    int pos = -1;
};

using Msg = std::array<double, 4>;

// ==============================================================================
// 1. AListファイルの読み込み関数
// ==============================================================================
static bool load_alist(const std::string& filename, int& N, int& M,
                       std::vector<std::vector<int>>& checks,
                       std::vector<std::vector<EdgeRef>>& var_to_checks) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open " << filename << "\n";
        return false;
    }

    std::string line;
    std::getline(file, line);
    std::stringstream ss1(line);
    ss1 >> N >> M;

    std::getline(file, line);
    std::getline(file, line);
    std::getline(file, line);

    for (int i = 0; i < N; ++i) {
        std::getline(file, line);
    }

    checks.assign(M, std::vector<int>());
    for (int i = 0; i < M; ++i) {
        std::getline(file, line);
        std::stringstream ss(line);
        int val;
        while (ss >> val) {
            if (val > 0) {
                checks[i].push_back(val - 1);
            }
        }
    }

    var_to_checks.assign(N, std::vector<EdgeRef>());
    for (int c = 0; c < M; ++c) {
        for (int pos = 0; pos < static_cast<int>(checks[c].size()); ++pos) {
            int v = checks[c][pos];
            if (v >= 0 && v < N) {
                var_to_checks[v].push_back({c, pos});
            }
        }
    }

    return true;
}

// ==============================================================================
// 2. BP復号に関連するユーティリティ
// ==============================================================================
static void normalize_msg(Msg &m) {
    double sum = m[0] + m[1] + m[2] + m[3];
    if (sum <= 0.0) {
        m = {0.25, 0.25, 0.25, 0.25};
        return;
    }
    for (double &v : m) v /= sum;
}

static double llr_from_probs(double p0, double p1, double eps) {
    if (p0 < eps) p0 = eps;
    if (p1 < eps) p1 = eps;
    return std::log(p0 / p1);
}

static double safe_atanh(double x, double eps) {
    if (x >= 1.0) return std::atanh(1.0 - eps);
    if (x <= -1.0) return std::atanh(-1.0 + eps);
    return std::atanh(x);
}

static double logsumexp2(double a, double b) {
    if (a < b) std::swap(a, b);
    return a + std::log1p(std::exp(b - a));
}

static int xbit(int state) { return (state == 1 || state == 3) ? 1 : 0; }
static int zbit(int state) { return (state == 2 || state == 3) ? 1 : 0; }

static std::vector<int> compute_syndrome(
    const std::vector<std::vector<int>> &checks,
    const std::vector<int> &err,
    bool use_xbit) {
    std::vector<int> synd(checks.size(), 0);
    for (int c = 0; c < static_cast<int>(checks.size()); ++c) {
        int parity = 0;
        for (int v : checks[c]) {
            int bit = use_xbit ? xbit(err[v]) : zbit(err[v]);
            parity ^= bit;
        }
        synd[c] = parity;
    }
    return synd;
}

static void update_check_llr(
    const std::vector<std::vector<double>> &v2c,
    const std::vector<int> &synd,
    std::vector<std::vector<double>> &c2v,
    double eps) {
    for (size_t c = 0; c < v2c.size(); ++c) {
        size_t deg = v2c[c].size();
        if (deg == 0) continue;
        std::vector<double> t(deg);
        for (size_t i = 0; i < deg; ++i) {
            t[i] = std::tanh(0.5 * v2c[c][i]);
        }
        std::vector<double> prefix(deg + 1, 1.0);
        std::vector<double> suffix(deg + 1, 1.0);
        for (size_t i = 0; i < deg; ++i) prefix[i + 1] = prefix[i] * t[i];
        for (size_t i = deg; i-- > 0;) suffix[i] = suffix[i + 1] * t[i];
        
        double sign = (synd[c] == 0) ? 1.0 : -1.0;
        for (size_t i = 0; i < deg; ++i) {
            double prod_excl = prefix[i] * suffix[i + 1];
            double arg = sign * prod_excl;
            c2v[c][i] = 2.0 * safe_atanh(arg, eps);
        }
    }
}

struct JointBPResult {
    std::vector<int> est;
    int iterations = 0;
    bool syndrome_match = false;
};

// ==============================================================================
// 3. 純粋な Joint BP デコーダ (詳細ログ出力付き)
// ==============================================================================
static JointBPResult joint_bp_decode(
    const std::vector<std::vector<int>> &x_checks,
    const std::vector<std::vector<int>> &z_checks,
    const std::vector<std::vector<EdgeRef>> &var_to_x,
    const std::vector<std::vector<EdgeRef>> &var_to_z,
    const std::vector<int> &sx,
    const std::vector<int> &sz,
    const Msg &prior,
    int max_iter,
    bool verbose) {
    
    int mX = static_cast<int>(x_checks.size());
    int mZ = static_cast<int>(z_checks.size());
    int nvars = static_cast<int>(var_to_x.size());

    const double prior_eps = 1e-300;
    double log_phi00 = std::log(std::max(prior[0], prior_eps));
    double log_phi10 = std::log(std::max(prior[1], prior_eps));
    double log_phi01 = std::log(std::max(prior[2], prior_eps));
    double log_phi11 = std::log(std::max(prior[3], prior_eps));

    double Lg_x0 = llr_from_probs(prior[0] + prior[2], prior[1] + prior[3], prior_eps);
    double Lg_z0 = llr_from_probs(prior[0] + prior[1], prior[2] + prior[3], prior_eps);

    std::vector<std::vector<double>> x_v2c(mX), x_c2v(mX);
    for (int c = 0; c < mX; ++c) {
        x_v2c[c].assign(x_checks[c].size(), Lg_x0);
        x_c2v[c].assign(x_checks[c].size(), 0.0);
    }

    std::vector<std::vector<double>> z_v2c(mZ), z_c2v(mZ);
    for (int c = 0; c < mZ; ++c) {
        z_v2c[c].assign(z_checks[c].size(), Lg_z0);
        z_c2v[c].assign(z_checks[c].size(), 0.0);
    }

    std::vector<int> est(nvars, 0);

    for (int iter = 0; iter < max_iter; ++iter) {
        update_check_llr(x_v2c, sx, x_c2v, 1e-12);
        update_check_llr(z_v2c, sz, z_c2v, 1e-12);

        for (int v = 0; v < nvars; ++v) {
            double Lx_sum = 0.0;
            for (const auto &e : var_to_x[v]) Lx_sum += x_c2v[e.check][e.pos];
            
            double Lz_sum = 0.0;
            for (const auto &e : var_to_z[v]) Lz_sum += z_c2v[e.check][e.pos];

            double Lx_half = 0.5 * Lx_sum;
            double Lz_half = 0.5 * Lz_sum;

            double log_x0 = Lx_half;
            double log_x1 = -Lx_half;
            double log_z0 = Lz_half;
            double log_z1 = -Lz_half;

            double Lg_to_x = logsumexp2(log_phi00 + log_z0, log_phi01 + log_z1) -
                             logsumexp2(log_phi10 + log_z0, log_phi11 + log_z1);
            double Lg_to_z = logsumexp2(log_phi00 + log_x0, log_phi10 + log_x1) -
                             logsumexp2(log_phi01 + log_x0, log_phi11 + log_x1);

            double logI = log_phi00 + log_x0 + log_z0;
            double logX = log_phi10 + log_x1 + log_z0;
            double logZ = log_phi01 + log_x0 + log_z1;
            double logY = log_phi11 + log_x1 + log_z1;

            int best = 0;
            double best_val = logI;
            if (logX > best_val) { best = 1; best_val = logX; }
            if (logZ > best_val) { best = 2; best_val = logZ; }
            if (logY > best_val) { best = 3; }
            est[v] = best;

            for (const auto &e : var_to_x[v]) {
                x_v2c[e.check][e.pos] = Lg_to_x + (Lx_sum - x_c2v[e.check][e.pos]);
            }
            for (const auto &e : var_to_z[v]) {
                z_v2c[e.check][e.pos] = Lg_to_z + (Lz_sum - z_c2v[e.check][e.pos]);
            }
        }

        auto sx_hat = compute_syndrome(x_checks, est, true);
        auto sz_hat = compute_syndrome(z_checks, est, false);
        
        int us_x = 0, us_z = 0;
        for (size_t c = 0; c < x_checks.size(); ++c) {
            if (sx_hat[c] != sx[c]) us_x++;
        }
        for (size_t c = 0; c < z_checks.size(); ++c) {
            if (sz_hat[c] != sz[c]) us_z++;
        }

        if (verbose) {
            std::cout << "  Iter " << std::setw(3) << iter + 1 
                      << " | Unsatisfied X: " << std::setw(5) << us_x 
                      << " | Unsatisfied Z: " << std::setw(5) << us_z << "\n";
        }
        
        if (us_x == 0 && us_z == 0) {
            return {est, iter + 1, true};
        }
    }
    return {est, max_iter, false};
}

// ==============================================================================
// メイン関数
// ==============================================================================
int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " --hx <H_X.alist> --hz <H_Z.alist> --p <error_rate> --trials <trials> [--max-iter <max_iter>] [--verbose]\n";
        return 1;
    }

    std::string hx_path, hz_path;
    double p_err = 0.05;
    int max_iter = 50;
    int trials = 1000;
    bool verbose = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--hx" && i + 1 < argc) hx_path = argv[++i];
        else if (arg == "--hz" && i + 1 < argc) hz_path = argv[++i];
        else if (arg == "--p" && i + 1 < argc) p_err = std::stod(argv[++i]);
        else if (arg == "--trials" && i + 1 < argc) trials = std::stoi(argv[++i]);
        else if (arg == "--max-iter" && i + 1 < argc) max_iter = std::stoi(argv[++i]);
        else if (arg == "--verbose") verbose = true;
    }

    if (hx_path.empty() || hz_path.empty()) {
        std::cerr << "Error: --hx and --hz are required.\n";
        return 1;
    }

    int N_x, M_x, N_z, M_z;
    std::vector<std::vector<int>> x_checks, z_checks;
    std::vector<std::vector<EdgeRef>> var_to_x, var_to_z;

    if (!load_alist(hx_path, N_x, M_x, x_checks, var_to_x)) return 1;
    if (!load_alist(hz_path, N_z, M_z, z_checks, var_to_z)) return 1;

    int nvars = N_x;
    std::cout << "Loaded H_X: " << M_x << " x " << nvars << "\n";
    std::cout << "Loaded H_Z: " << M_z << " x " << nvars << "\n";
    std::cout << "p = " << p_err << ", trials = " << trials << ", max_iter = " << max_iter << "\n\n";

    Msg prior{1.0 - p_err, p_err / 3.0, p_err / 3.0, p_err / 3.0};
    normalize_msg(prior);

    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    int failures = 0;
    long long total_iters = 0;

    auto start_time = std::chrono::steady_clock::now();

    for (int t = 0; t < trials; ++t) {
        std::vector<int> err(nvars, 0);
        int actual_errors = 0;
        for (int i = 0; i < nvars; ++i) {
            double u = dist(rng);
            if (u >= (1.0 - p_err)) {
                double v = (u - (1.0 - p_err)) / p_err;
                if (v < 1.0 / 3.0) err[i] = 1;
                else if (v < 2.0 / 3.0) err[i] = 2;
                else err[i] = 3;
                actual_errors++;
            }
        }

        if (verbose) {
            std::cout << "=== Trial " << t + 1 << " / " << trials << " ===\n";
            std::cout << "Initial Error Weight: " << actual_errors << "\n";
        }

        auto sx = compute_syndrome(x_checks, err, true);
        auto sz = compute_syndrome(z_checks, err, false);

        // --verbose フラグが立っていれば、内部の進捗を表示する
        auto res = joint_bp_decode(x_checks, z_checks, var_to_x, var_to_z, sx, sz, prior, max_iter, verbose);

        total_iters += res.iterations;

        if (!res.syndrome_match) {
            failures++;
        }

        if (!verbose && (t + 1) % 100 == 0) {
            std::cout << "Progress: " << (t + 1) << " / " << trials 
                      << " | Failures: " << failures 
                      << " | FER: " << (static_cast<double>(failures) / (t + 1)) << "\n";
        }
    }

    auto end_time = std::chrono::steady_clock::now();
    double elapsed_sec = std::chrono::duration<double>(end_time - start_time).count();
    double fer = static_cast<double>(failures) / trials;

    std::cout << "\n=== Simulation Finished ===\n";
    std::cout << "Trials     : " << trials << "\n";
    std::cout << "Failures   : " << failures << "\n";
    std::cout << "FER        : " << fer << "\n";
    std::cout << "Elapsed    : " << elapsed_sec << " s\n";

    return 0;
}