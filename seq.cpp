#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <chrono>
#include "common.h"

int main(int argc, char** argv) {
    int M = argc > 1 ? atoi(argv[1]) : 80;
    int N = argc > 2 ? atoi(argv[2]) : 80;
    double delta = argc > 3 ? atof(argv[3]) : 1e-8;
    const char* out = (argc > 4 && argv[4][0] && argv[4][0] != '-') ? argv[4] : nullptr;
    const char* trace = (argc > 5 && argv[5][0] && argv[5][0] != '-') ? argv[5] : nullptr;
    FILE* ft = trace ? fopen(trace, "w") : nullptr;
    const char* ref = argc > 6 ? argv[6] : nullptr;

    auto t0 = std::chrono::steady_clock::now();

    const double h1 = (prob::B1 - prob::A1) / M, h2 = (prob::B2 - prob::A2) / N;
    const double h = std::max(h1, h2), eps = h * h;
    const double ih1 = 1.0 / (h1 * h1), ih2 = 1.0 / (h2 * h2), hh = h1 * h2;
    const int ld = N + 1;
    const size_t sz = (size_t)(M + 1) * (N + 1);
    auto id = [ld](int i, int j) { return (size_t)i * ld + j; };

    std::vector<double> Q(sz, 0), B(sz, 0), Dg(sz, 1), w(sz, 0), r(sz, 0),
                        z(sz, 0), p(sz, 0), Ap(sz, 0);

    for (int i = 1; i < M; ++i)
        for (int j = 1; j < N; ++j) {
            double x = prob::A1 + i * h1, y = prob::A2 + j * h2;
            prob::Coef c = prob::coef(x, y, h1, h2, eps);
            double b = c.rhs;
            if (i == 1)     b += ih1 * prob::v0(prob::A1, y);
            if (i == M - 1) b += ih1 * prob::v0(prob::B1, y);
            if (j == 1)     b += ih2 * prob::v0(x, prob::A2);
            if (j == N - 1) b += ih2 * prob::v0(x, prob::B2);
            Q[id(i, j)] = c.Q;
            B[id(i, j)] = b;
            Dg[id(i, j)] = 2 * ih1 + 2 * ih2 + c.Q;
        }

    std::vector<double> R;
    std::vector<char> mask(sz, 0);
    for (int i = 1; i < M; ++i)
        for (int j = 1; j < N; ++j)
            mask[id(i, j)] = prob::inD(prob::A1 + i * h1, prob::A2 + j * h2);
    if (ref) {
        R.assign(sz, 0.0);
        if (std::string(ref) == "exact") {
            for (int i = 0; i <= M; ++i)
                for (int j = 0; j <= N; ++j)
                    R[id(i, j)] = prob::exact(prob::A1 + i * h1, prob::A2 + j * h2);
        } else {
            FILE* fr = fopen(ref, "r");
            int rm = 0, rn = 0;
            if (!fr || fscanf(fr, "%d %d", &rm, &rn) != 2 || rm != M || rn != N) {
                fprintf(stderr, "bad reference file %s\n", ref);
                return 1;
            }
            double xx, yy;
            for (size_t k = 0; k < sz; ++k)
                if (fscanf(fr, "%lf %lf %lf", &xx, &yy, &R[k]) != 3) { fprintf(stderr, "short reference\n"); return 1; }
            fclose(fr);
        }
    }
    auto errD = [&](const std::vector<double>& u, double& eC, double& eE) {
        eC = 0; eE = 0;
        for (int i = 1; i < M; ++i)
            for (int j = 1; j < N; ++j) {
                size_t k = id(i, j);
                if (!mask[k]) continue;
                double d = std::fabs(u[k] - R[k]);
                eC = std::max(eC, d);
                eE += d * d;
            }
        eE = std::sqrt(eE * hh);
    };

    auto applyA = [&](const std::vector<double>& u, std::vector<double>& res) {
        for (int i = 1; i < M; ++i)
            for (int j = 1; j < N; ++j) {
                size_t k = id(i, j);
                res[k] = -(u[k + ld] - 2 * u[k] + u[k - ld]) * ih1
                         -(u[k + 1]  - 2 * u[k] + u[k - 1])  * ih2 + Q[k] * u[k];
            }
    };
    auto dot = [&](const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0;
        for (int i = 1; i < M; ++i)
            for (int j = 1; j < N; ++j) s += a[id(i, j)] * b[id(i, j)];
        return s * hh;
    };

    for (size_t k = 0; k < sz; ++k) { r[k] = B[k]; z[k] = r[k] / Dg[k]; p[k] = z[k]; }
    double zr = dot(z, r);
    long maxit = (long)(M - 1) * (N - 1);
    int it = 0;
    double diff = 0;
    while (it < maxit) {
        applyA(p, Ap);
        double alpha = zr / dot(Ap, p);
        for (int i = 1; i < M; ++i)
            for (int j = 1; j < N; ++j) {
                size_t k = id(i, j);
                w[k] += alpha * p[k];
                r[k] -= alpha * Ap[k];
                z[k] = r[k] / Dg[k];
            }
        ++it;
        diff = std::fabs(alpha) * std::sqrt(dot(p, p));
        if (ft) {
            double J = 0;
            for (int i = 1; i < M; ++i)
                for (int j = 1; j < N; ++j) J += (B[id(i, j)] + r[id(i, j)]) * w[id(i, j)];
            double eC = -1, eE = -1;
            if (ref) errD(w, eC, eE);
            fprintf(ft, "%d %.6e %.15e %.6e %.6e %.6e\n", it, diff, -0.5 * J * hh,
                    std::sqrt(dot(r, r)), eE, eC);
        }
        if (diff < delta) break;
        double zr_new = dot(z, r);
        double beta = zr_new / zr;
        zr = zr_new;
        for (int i = 1; i < M; ++i)
            for (int j = 1; j < N; ++j) {
                size_t k = id(i, j);
                p[k] = z[k] + beta * p[k];
            }
    }
    for (int i = 0; i <= M; ++i) {
        w[id(i, 0)] = prob::v0(prob::A1 + i * h1, prob::A2);
        w[id(i, N)] = prob::v0(prob::A1 + i * h1, prob::B2);
    }
    for (int j = 0; j <= N; ++j) {
        w[id(0, j)] = prob::v0(prob::A1, prob::A2 + j * h2);
        w[id(M, j)] = prob::v0(prob::B1, prob::A2 + j * h2);
    }
    if (ft) fclose(ft);
    double time = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    printf("SEQ M=%d N=%d iters=%d time=%.4f diff=%.3e\n", M, N, it, time, diff);
    if (ref) {
        double eC, eE;
        errD(w, eC, eE);
        printf("ERR_D C=%.6e E=%.6e\n", eC, eE);
    }
    if (out) {
        FILE* fo = fopen(out, "w");
        fprintf(fo, "%d %d\n", M, N);
        for (int i = 0; i <= M; ++i)
            for (int j = 0; j <= N; ++j)
                fprintf(fo, "%.6f %.6f %.12e\n", prob::A1 + i * h1, prob::A2 + j * h2, w[id(i, j)]);
        fclose(fo);
    }
    return 0;
}
