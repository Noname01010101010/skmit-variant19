// Проверка алгоритма разбиения на домены (partition.h).
#include <cstdio>
#include <vector>
#include "../partition.h"

int main() {
    int fails = 0, cases = 0, notOk = 0;
    const int grids[][2] = {{79, 79}, {399, 599}, {799, 1199}, {36, 52}, {100, 30}, {7, 300}};
    for (auto& g : grids)
        for (int P = 1; P <= 64; ++P) {
            int nx = g[0], ny = g[1];
            if (P > nx * ny) continue;
            int px, py;
            choose_grid(P, nx, ny, px, py);
            ++cases;
            bool ok = px * py == P;
            // покрытие: домены по x и y идут подряд без пропусков и пересечений
            int prev = 0, mnx = 1 << 30, mxx = 0, mny = 1 << 30, mxy = 0;
            for (int k = 0; k < px; ++k) {
                int lo, hi; split(nx, px, k, 1, lo, hi);
                ok &= lo == prev + 1; prev = hi;
                mnx = std::min(mnx, hi - lo + 1); mxx = std::max(mxx, hi - lo + 1);
            }
            ok &= prev == nx;
            prev = 0;
            for (int k = 0; k < py; ++k) {
                int lo, hi; split(ny, py, k, 1, lo, hi);
                ok &= lo == prev + 1; prev = hi;
                mny = std::min(mny, hi - lo + 1); mxy = std::max(mxy, hi - lo + 1);
            }
            ok &= prev == ny;
            ok &= mxx - mnx <= 1 && mxy - mny <= 1;            // размеры отличаются не более чем на 1
            bool ratio = (double)mnx / mxy >= 0.5 && (double)mxx / mny <= 2.0;
            if (!ratio) {
                ++notOk;
                // алгоритм не должен пропускать допустимое разбиение, если оно существует
                for (int a = 1; a <= P; ++a) {
                    if (P % a || a > nx || P / a > ny) continue;
                    int b = P / a, x0 = nx / a, x1 = x0 + (nx % a != 0), y0 = ny / b, y1 = y0 + (ny % b != 0);
                    if ((double)x0 / y1 >= 0.5 && (double)x1 / y0 <= 2.0) {
                        ok = false;
                        printf("[FAIL] grid %dx%d P=%d: missed valid split %dx%d\n", nx, ny, P, a, b);
                    }
                }
            }
            if (!ok) {
                ++fails;
                printf("[FAIL] grid %dx%d P=%d -> %dx%d\n", nx, ny, P, px, py);
            }
        }
    // конфигурации из таблиц отчета обязаны удовлетворять условию [1/2, 2]
    const int tab[][3] = {{79,79,4},{399,599,2},{399,599,4},{399,599,8},{399,599,16},
                          {799,1199,4},{799,1199,8},{799,1199,16},{799,1199,32}};
    for (auto& t : tab) {
        int px, py; choose_grid(t[2], t[0], t[1], px, py);
        int mnx = t[0] / px, mxx = mnx + (t[0] % px != 0), mny = t[1] / py, mxy = mny + (t[1] % py != 0);
        bool r = (double)mnx / mxy >= 0.5 && (double)mxx / mny <= 2.0;
        printf("%s grid %4dx%-4d P=%2d -> %dx%d, x/y in [%.2f, %.2f]\n", r ? "[ OK ]" : "[FAIL]",
               t[0], t[1], t[2], px, py, (double)mnx / mxy, (double)mxx / mny);
        if (!r) ++fails;
    }
    printf("partition: %d cases checked, %d without exact [1/2,2] (checked: no such split exists), %d failed\n",
           cases, notOk, fails);
    return fails ? 1 : 0;
}
