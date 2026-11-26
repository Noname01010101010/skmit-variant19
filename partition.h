#pragma once
#include <cmath>
#include <algorithm>

// Двумерное разбиение P = px * py для сетки nx x ny внутренних узлов.
inline void choose_grid(int P, int nx, int ny, int& px, int& py) {
    double best = 1e300; bool bestOk = false;
    px = P; py = 1;
    for (int a = 1; a <= P; ++a) {
        if (P % a) continue;
        int b = P / a;
        if (a > nx || b > ny) continue;
        int xmin = nx / a, xmax = xmin + (nx % a ? 1 : 0);
        int ymin = ny / b, ymax = ymin + (ny % b ? 1 : 0);
        bool ok = (double)xmin / ymax >= 0.5 && (double)xmax / ymin <= 2.0;
        double score = std::fabs(std::log(((double)nx / a) / ((double)ny / b)));
        if ((ok && !bestOk) || (ok == bestOk && score < best)) {
            best = score; bestOk = ok; px = a; py = b;
        }
    }
}

// Узлы first..first+n-1 делятся на p частей; часть k получает [lo, hi].
inline void split(int n, int p, int k, int first, int& lo, int& hi) {
    int base = n / p, rem = n % p;
    lo = first + k * base + std::min(k, rem);
    hi = lo + base + (k < rem ? 1 : 0) - 1;
}
