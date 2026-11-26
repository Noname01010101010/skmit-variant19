#include <cstdio>
#include <cstdlib>
#include <vector>
#include <omp.h>
#include "common.h"

int main(int argc, char** argv) {
    int M = argc > 1 ? atoi(argv[1]) : 80;
    int N = argc > 2 ? atoi(argv[2]) : 80;
    double delta = argc > 3 ? atof(argv[3]) : 1e-8;
    const char* out = (argc > 4 && argv[4][0] && argv[4][0] != '-') ? argv[4] : nullptr;

    double t0 = omp_get_wtime();

    const double h1 = (prob::B1 - prob::A1) / M, h2 = (prob::B2 - prob::A2) / N;
    const double h = std::max(h1, h2), eps = h * h;
    const double ih1 = 1.0 / (h1 * h1), ih2 = 1.0 / (h2 * h2), hh = h1 * h2;
    const int ld = N + 1;
    const size_t sz = (size_t)(M + 1) * (N + 1);
    auto id = [ld](int i, int j) { return (size_t)i * ld + j; };

    std::vector<double> Q(sz, 0), B(sz, 0), Dg(sz, 1), w(sz, 0), r(sz, 0),
                        z(sz, 0), p(sz, 0), Ap(sz, 0);

#pragma omp parallel for schedule(static)
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

    auto applyA = [&](const std::vector<double>& u, std::vector<double>& res) {
#pragma omp parallel for schedule(static)
        for (int i = 1; i < M; ++i)
            for (int j = 1; j < N; ++j) {
                size_t k = id(i, j);
                res[k] = -(u[k + ld] - 2 * u[k] + u[k - ld]) * ih1
                         -(u[k + 1]  - 2 * u[k] + u[k - 1])  * ih2 + Q[k] * u[k];
            }
    };
    auto dot = [&](const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0;
#pragma omp parallel for schedule(static) reduction(+:s)
        for (int i = 1; i < M; ++i)
            for (int j = 1; j < N; ++j) s += a[id(i, j)] * b[id(i, j)];
        return s * hh;
    };

#pragma omp parallel for schedule(static)
    for (size_t k = 0; k < sz; ++k) { r[k] = B[k]; z[k] = r[k] / Dg[k]; p[k] = z[k]; }
    double zr = dot(z, r);
    long maxit = (long)(M - 1) * (N - 1);
    int it = 0;
    double diff = 0;
    while (it < maxit) {
        applyA(p, Ap);
        double alpha = zr / dot(Ap, p);
#pragma omp parallel for schedule(static)
        for (int i = 1; i < M; ++i)
            for (int j = 1; j < N; ++j) {
                size_t k = id(i, j);
                w[k] += alpha * p[k];
                r[k] -= alpha * Ap[k];
                z[k] = r[k] / Dg[k];
            }
        ++it;
        diff = std::fabs(alpha) * std::sqrt(dot(p, p));
        if (diff < delta) break;
        double zr_new = dot(z, r);
        double beta = zr_new / zr;
        zr = zr_new;
#pragma omp parallel for schedule(static)
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
    double time = omp_get_wtime() - t0;

    printf("OMP threads=%d M=%d N=%d iters=%d time=%.4f diff=%.3e\n", omp_get_max_threads(), M, N, it, time, diff);
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
