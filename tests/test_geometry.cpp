// Проверка вычисления площадей и коэффициентов схемы (common.h).
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "../common.h"

static int fails = 0;
static void check(bool ok, const char* what, double got, double want) {
    printf("%s %-58s got=%.12g want=%.12g\n", ok ? "[ OK ]" : "[FAIL]", what, got, want);
    if (!ok) ++fails;
}

// площадь пересечения прямоугольника с полукругом по средним точкам мелкой сетки
static double areaBrute(double xa, double xb, double ya, double yb, int n) {
    double hx = (xb - xa) / n, hy = (yb - ya) / n, s = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double x = xa + (i + 0.5) * hx, y = ya + (j + 0.5) * hy;
            if (y < 2.0 && x * x + (y - 2) * (y - 2) < 1.0) s += hx * hy;
        }
    return s;
}

int main() {
    const double PI = std::acos(-1.0);
    double a = prob::areaH(-1, 1, 1, 2);
    check(std::fabs(a - PI / 2) < 1e-14, "площадь всего полукруга = pi/2", a, PI / 2);
    a = prob::areaH(-2, 2, -2, 2);
    check(std::fabs(a - PI / 2) < 1e-14, "площадь полукруга внутри всего квадрата = pi/2", a, PI / 2);
    a = prob::areaH(-2, 2, -2, 1);
    check(a == 0.0, "ниже y=1 полукруга нет", a, 0);
    a = prob::areaH(1, 2, 1, 2);
    check(a == 0.0, "правее x=1 полукруга нет", a, 0);
    a = prob::areaH(-1, 1, 2, 3);
    check(a == 0.0, "выше y=2 (верхний полукруг вне квадрата) не считается", a, 0);
    a = prob::areaH(0, 1, 1, 2);
    check(std::fabs(a - PI / 4) < 1e-14, "четверть круга = pi/4", a, PI / 4);

    for (int n : {20, 37, 80, 400}) {
        double h1 = 4.0 / n, h2 = 4.0 / (n + 3), s = 0;
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n + 3; ++j) {
                double x = -2 + i * h1, y = -2 + j * h2;
                s += prob::areaH(x - h1 / 2, x + h1 / 2, y - h2 / 2, y + h2 / 2);
            }
        char buf[96];
        snprintf(buf, sizeof buf, "сумма площадей по ячейкам сетки %dx%d = pi/2", n, n + 3);
        check(std::fabs(s - PI / 2) < 1e-12, buf, s, PI / 2);
    }

    srand(12345);
    double worst = 0;   // отношение ошибки к допустимой погрешности численного подсчета
    for (int t = 0; t < 200; ++t) {
        double xa = -1.3 + 2.4 * rand() / RAND_MAX, ya = 0.8 + 1.1 * rand() / RAND_MAX;
        double w = 0.05 + 0.4 * rand() / RAND_MAX, h = 0.05 + 0.4 * rand() / RAND_MAX;
        double e = prob::areaH(xa, xa + w, ya, ya + h);
        const int n = 1500;
        double b = areaBrute(xa, xa + w, ya, ya + h, n);
        double tol = 2.0 * (w + h) * std::max(w, h) / n;   // погрешность подсчета по средним точкам
        worst = std::max(worst, std::fabs(e - b) / tol);
    }
    check(worst < 1.0, "200 случайных ячеек: |точная - численная| / допуск < 1", worst, 1);

    double h = 0.1, eps = h * h;
    prob::Coef c = prob::coef(-1.5, -1.5, h, h, eps);
    check(c.Q == 0.0 && std::fabs(c.rhs - prob::f(-1.5, -1.5)) < 1e-15,
          "ячейка внутри D: Q=0, правая часть = f", c.Q, 0);
    c = prob::coef(0.0, 1.5, h, h, eps);
    check(std::fabs(c.Q - 1 / (eps * eps)) < 1e-6 / (eps * eps) &&
          std::fabs(c.rhs - c.Q * prob::v0(0.0, 1.5)) < 1e-9 * c.Q,
          "ячейка в фиктивной области: Q=1/eps^2, F=0", c.Q * eps * eps, 1);
    c = prob::coef(0.0, 1.0, h, h, eps);
    check(c.Q > 0 && c.Q < 1 / (eps * eps), "ячейка на дуге: 0 < Q < 1/eps^2", c.Q * eps * eps, 0.5);

    printf(fails ? "geometry: %d FAILED\n" : "geometry: all passed\n", fails);
    return fails ? 1 : 0;
}
