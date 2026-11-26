#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "common.h"
#include "partition.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int M = argc > 1 ? atoi(argv[1]) : 80;
    int N = argc > 2 ? atoi(argv[2]) : 80;
    double delta = argc > 3 ? atof(argv[3]) : 1e-8;
    const char* out = (argc > 4 && argv[4][0] && argv[4][0] != '-') ? argv[4] : nullptr;

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    int px, py;
    choose_grid(size, M - 1, N - 1, px, py);
    int dims[2] = {px, py}, periods[2] = {0, 0}, coords[2];
    MPI_Comm cart;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart);
    MPI_Comm_rank(cart, &rank);
    MPI_Cart_coords(cart, rank, 2, coords);
    int left, right, down, up;
    MPI_Cart_shift(cart, 0, 1, &left, &right);
    MPI_Cart_shift(cart, 1, 1, &down, &up);

    int i0, i1, j0, j1;
    split(M - 1, px, coords[0], 1, i0, i1);
    split(N - 1, py, coords[1], 1, j0, j1);
    const int nx = i1 - i0 + 1, ny = j1 - j0 + 1;
    const int ld = ny + 2;
    const size_t sz = (size_t)(nx + 2) * (ny + 2);
    auto id = [ld](int li, int lj) { return (size_t)li * ld + lj; };

    const double h1 = (prob::B1 - prob::A1) / M, h2 = (prob::B2 - prob::A2) / N;
    const double h = std::max(h1, h2), eps = h * h;
    const double ih1 = 1.0 / (h1 * h1), ih2 = 1.0 / (h2 * h2), hh = h1 * h2;

    std::vector<double> Q(sz, 0), B(sz, 0), Dg(sz, 1), w(sz, 0), r(sz, 0),
                        z(sz, 0), p(sz, 0), Ap(sz, 0);

    for (int li = 1; li <= nx; ++li)
        for (int lj = 1; lj <= ny; ++lj) {
            int i = i0 + li - 1, j = j0 + lj - 1;
            double x = prob::A1 + i * h1, y = prob::A2 + j * h2;
            prob::Coef c = prob::coef(x, y, h1, h2, eps);
            double b = c.rhs;
            if (i == 1)     b += ih1 * prob::v0(prob::A1, y);
            if (i == M - 1) b += ih1 * prob::v0(prob::B1, y);
            if (j == 1)     b += ih2 * prob::v0(x, prob::A2);
            if (j == N - 1) b += ih2 * prob::v0(x, prob::B2);
            size_t k = id(li, lj);
            Q[k] = c.Q; B[k] = b; Dg[k] = 2 * ih1 + 2 * ih2 + c.Q;
        }

    MPI_Datatype colType;
    MPI_Type_vector(nx, 1, ld, MPI_DOUBLE, &colType);
    MPI_Type_commit(&colType);

    auto exchange = [&](std::vector<double>& u) {
        double* a = u.data();
        MPI_Sendrecv(a + id(nx, 1), ny, MPI_DOUBLE, right, 0,
                     a + id(0, 1),  ny, MPI_DOUBLE, left,  0, cart, MPI_STATUS_IGNORE);
        MPI_Sendrecv(a + id(1, 1),      ny, MPI_DOUBLE, left,  1,
                     a + id(nx + 1, 1), ny, MPI_DOUBLE, right, 1, cart, MPI_STATUS_IGNORE);
        MPI_Sendrecv(a + id(1, ny), 1, colType, up,   2,
                     a + id(1, 0),  1, colType, down, 2, cart, MPI_STATUS_IGNORE);
        MPI_Sendrecv(a + id(1, 1),      1, colType, down, 3,
                     a + id(1, ny + 1), 1, colType, up,   3, cart, MPI_STATUS_IGNORE);
    };
    auto applyA = [&](std::vector<double>& u, std::vector<double>& res) {
        exchange(u);
        for (int li = 1; li <= nx; ++li)
            for (int lj = 1; lj <= ny; ++lj) {
                size_t k = id(li, lj);
                res[k] = -(u[k + ld] - 2 * u[k] + u[k - ld]) * ih1
                         -(u[k + 1]  - 2 * u[k] + u[k - 1])  * ih2 + Q[k] * u[k];
            }
    };
    auto dot = [&](const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0, g = 0;
        for (int li = 1; li <= nx; ++li)
            for (int lj = 1; lj <= ny; ++lj) s += a[id(li, lj)] * b[id(li, lj)];
        MPI_Allreduce(&s, &g, 1, MPI_DOUBLE, MPI_SUM, cart);
        return g * hh;
    };

    for (size_t k = 0; k < sz; ++k) { r[k] = B[k]; z[k] = r[k] / Dg[k]; p[k] = z[k]; }
    double zr = dot(z, r);
    long maxit = (long)(M - 1) * (N - 1);
    int it = 0;
    double diff = 0;
    while (it < maxit) {
        applyA(p, Ap);
        double alpha = zr / dot(Ap, p);
        for (int li = 1; li <= nx; ++li)
            for (int lj = 1; lj <= ny; ++lj) {
                size_t k = id(li, lj);
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
        for (int li = 1; li <= nx; ++li)
            for (int lj = 1; lj <= ny; ++lj) {
                size_t k = id(li, lj);
                p[k] = z[k] + beta * p[k];
            }
    }
    double tloc = MPI_Wtime() - t0, time;
    MPI_Reduce(&tloc, &time, 1, MPI_DOUBLE, MPI_MAX, 0, cart);

    if (rank == 0) {
        printf("MPI procs=%d (%dx%d) M=%d N=%d iters=%d time=%.4f diff=%.3e\n",
               size, px, py, M, N, it, time, diff);
    }

    if (out) {
        std::vector<double> loc((size_t)nx * ny);
        for (int li = 1; li <= nx; ++li)
            for (int lj = 1; lj <= ny; ++lj) loc[(size_t)(li - 1) * ny + lj - 1] = w[id(li, lj)];
        if (rank == 0) {
            std::vector<double> W((size_t)(M + 1) * (N + 1), 0.0);
            auto put = [&](int a0, int a1, int b0, int b1, const double* d) {
                int cy = b1 - b0 + 1;
                for (int i = a0; i <= a1; ++i)
                    for (int j = b0; j <= b1; ++j)
                        W[(size_t)i * (N + 1) + j] = d[(size_t)(i - a0) * cy + (j - b0)];
            };
            put(i0, i1, j0, j1, loc.data());
            for (int src = 1; src < size; ++src) {
                int b[4];
                MPI_Recv(b, 4, MPI_INT, src, 10, cart, MPI_STATUS_IGNORE);
                std::vector<double> buf((size_t)(b[1] - b[0] + 1) * (b[3] - b[2] + 1));
                MPI_Recv(buf.data(), (int)buf.size(), MPI_DOUBLE, src, 11, cart, MPI_STATUS_IGNORE);
                put(b[0], b[1], b[2], b[3], buf.data());
            }
            for (int i = 0; i <= M; ++i) {
                W[(size_t)i * (N + 1)]     = prob::v0(prob::A1 + i * h1, prob::A2);
                W[(size_t)i * (N + 1) + N] = prob::v0(prob::A1 + i * h1, prob::B2);
            }
            for (int j = 0; j <= N; ++j) {
                W[j]                       = prob::v0(prob::A1, prob::A2 + j * h2);
                W[(size_t)M * (N + 1) + j] = prob::v0(prob::B1, prob::A2 + j * h2);
            }
            FILE* fo = fopen(out, "w");
            fprintf(fo, "%d %d\n", M, N);
            for (int i = 0; i <= M; ++i)
                for (int j = 0; j <= N; ++j)
                    fprintf(fo, "%.6f %.6f %.12e\n", prob::A1 + i * h1, prob::A2 + j * h2,
                            W[(size_t)i * (N + 1) + j]);
            fclose(fo);
        } else {
            int b[4] = {i0, i1, j0, j1};
            MPI_Send(b, 4, MPI_INT, 0, 10, cart);
            MPI_Send(loc.data(), (int)loc.size(), MPI_DOUBLE, 0, 11, cart);
        }
    }
    MPI_Type_free(&colType);
    MPI_Comm_free(&cart);
    MPI_Finalize();
    return 0;
}

