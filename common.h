#pragma once
#include <cmath>
#include <algorithm>

namespace prob {

const double A1 = -2.0, B1 = 2.0;
const double A2 = -2.0, B2 = 2.0;

#ifdef TEST_EXACT
// тестовая задача с известным решением u = sin(x+y): -Lap u = 2 sin(x+y)
inline double f(double x, double y)   { return 2.0 * std::sin(x + y); }
#else
inline double f(double, double)       { return 1.0; }
#endif
inline double exact(double x, double y) { return std::sin(x + y); }
inline double phi(double x, double y) { return std::sin(x + y); }
inline double v0(double x, double y)  { return std::sin(x + y); }

inline bool inD(double x, double y) {
    return x > A1 && x < B1 && y > A2 && y < B2 &&
           x * x + (y - 2.0) * (y - 2.0) > 1.0;
}

inline double P(double x) {
    x = std::max(-1.0, std::min(1.0, x));
    return 0.5 * (x * std::sqrt(1.0 - x * x) + std::asin(x));
}

inline double I_t(double t, double a, double b) {
    if (t <= 1.0) return 0.0;
    double d = 2.0 - t;
    if (d < 0.0) d = 0.0;
    double s = std::sqrt(1.0 - d * d);
    double lo = std::max(a, -s), hi = std::min(b, s);
    if (hi <= lo) return 0.0;
    return (t - 2.0) * (hi - lo) + P(hi) - P(lo);
}

// Площадь пересечения ячейки [xa,xb]x[ya,yb] с полукругом x^2+(y-2)^2<1, y<2
inline double areaH(double xa, double xb, double ya, double yb) {
    yb = std::min(yb, 2.0);
    if (yb <= ya) return 0.0;
    return I_t(yb, xa, xb) - I_t(ya, xa, xb);
}

struct Coef { double Q, rhs; };

inline Coef coef(double x, double y, double h1, double h2, double eps) {
    double xa = x - 0.5 * h1, xb = x + 0.5 * h1;
    double ya = y - 0.5 * h2, yb = y + 0.5 * h2;
    double cell = h1 * h2;
    double S = cell - areaH(xa, xb, ya, yb);
    if (S < 0) S = 0;
    double frac = S / cell;
    double Q = (1.0 - frac) / (eps * eps);
    double F;
    if (frac >= 1.0 - 1e-14) {
        F = f(x, y);
    } else if (frac <= 1e-14) {
        F = 0.0;
    } else {
        F = frac * f(x, y);
    }
    return {Q, F + Q * v0(x, y)};
}

}
