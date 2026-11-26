# Графики сходимости и ошибки (раздел отчета "Сходимость и погрешность").
# Запуск из корня репозитория:  python3 scripts/plot_convergence.py
# Собирает seq (обычный и с -DTEST_EXACT) в build/, считает решения и строит
# report/fig/{conv_cg,conv_iter_err,conv_grid,error_map,exact_grid}.png.
# Эталон 1280x1280 считается около минуты.
import os, subprocess, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

B, F = 'build', 'report/fig'
os.makedirs(B, exist_ok=True); os.makedirs(F, exist_ok=True)
subprocess.run(['g++', '-O1', '-std=c++11', 'seq.cpp', '-o', f'{B}/seq'], check=True)
subprocess.run(['g++', '-O1', '-std=c++11', '-DTEST_EXACT', 'seq.cpp', '-o', f'{B}/seq_exact'], check=True)

def seq(*args, exe='seq'):
    return subprocess.run([f'{B}/{exe}'] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout

def load(fn):
    with open(fn) as f: M, N = map(int, f.readline().split())
    d = np.loadtxt(fn, skiprows=1)
    return [d[:, c].reshape(M + 1, N + 1) for c in range(3)]

# 1. Журналы итераций: ||w(k+1)-w(k)||, J, невязка, ошибка относительно точного решения схемы
for n in (20, 40, 80):
    seq(n, n, 1e-12, f'{B}/ref{n}.txt')                           # решение схемы w*
    seq(n, n, 1e-8, '-', f'{B}/it{n}.txt', f'{B}/ref{n}.txt')
for n in (20, 40, 80, 160):
    seq(n, n, 1e-8, '-', f'{B}/itex{n}.txt', 'exact', exe='seq_exact')

fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.8))
for n in (20, 40, 80):
    k, d, J, r, eE, eC = np.loadtxt(f'{B}/it{n}.txt').T
    axs[0].semilogy(k, d, label=f'{n}×{n}')
    g = J - J[-1]; m = g > 1e-13 * abs(J[-1]); axs[1].semilogy(k[m], g[m], label=f'{n}×{n}')
axs[0].axhline(1e-8, color='k', ls='--', lw=0.8); axs[0].text(2, 1.5e-8, r'$\delta=10^{-8}$', fontsize=8)
axs[0].set(xlabel='итерация k', ylabel=r'$\|w^{(k+1)}-w^{(k)}\|_E$', title='Критерий остановки')
axs[1].set(xlabel='итерация k', ylabel=r'$J(w^{(k)})-J(w^{*})$', title='Монотонное убывание функционала J')
for a in axs: a.grid(alpha=.3); a.legend()
fig.tight_layout(); fig.savefig(f'{F}/conv_cg.png', dpi=160); plt.close(fig)

fig, axs = plt.subplots(1, 3, figsize=(13, 3.9))
for n in (20, 40, 80):
    k, d, J, r, eE, eC = np.loadtxt(f'{B}/it{n}.txt').T
    l, = axs[0].semilogy(k, eE, label=f'{n}×{n}, E')
    axs[0].semilogy(k, eC, '--', color=l.get_color(), label=f'{n}×{n}, C')
    axs[1].semilogy(k, r, label=f'{n}×{n}')
for n in (20, 40, 80, 160):
    k, d, J, r, eE, eC = np.loadtxt(f'{B}/itex{n}.txt').T
    axs[2].semilogy(k, eC, label=f'{n}×{n}')
axs[0].set_title(r'Ошибка итерации $\|w^{(k)}-w^{*}\|$ в $D$')
axs[1].set_title(r'Невязка $\|r^{(k)}\|_E=\|B-Aw^{(k)}\|_E$')
axs[2].set_title(r'Тест $u=\sin(x+y)$: $\|w^{(k)}-u\|_C$ в $D$')
for a in axs: a.set_xlabel('итерация k'); a.grid(alpha=.3); a.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f'{F}/conv_iter_err.png', dpi=160); plt.close(fig)

# 2. Сходимость по сетке в D относительно эталона 1280x1280
grids = (20, 40, 80, 160, 320, 640)
for n in grids + (1280,):
    seq(n, n, 1e-10, f'{B}/sol{n}.txt')
Xr, Yr, Wr = load(f'{B}/sol1280.txt')
res = []
for n in grids:
    X, Y, W = load(f'{B}/sol{n}.txt'); R = Wr[::1280 // n, ::1280 // n]
    inD = (X**2 + (Y - 2)**2 > 1) & (abs(X) < 2) & (abs(Y) < 2)
    e = (W - R)[inD]; h = 4 / n
    res.append((n, h, abs(e).max(), np.sqrt((e**2).sum() * h * h)))
res = np.array(res)
print('сетка  h       ||e||_C     ||e||_E')
for n, h, C, E in res: print(f'{int(n):5d} {h:.4f} {C:.3e} {E:.3e}')
fig, ax = plt.subplots(figsize=(5.6, 4.2))
ax.loglog(res[:, 1], res[:, 2], 'o-', label=r'$\|w^h-w^{ref}\|_C$ (в $D$)')
ax.loglog(res[:, 1], res[:, 3], 's-', label=r'$\|w^h-w^{ref}\|_E$ (в $D$)')
ax.loglog(res[:, 1], 0.8 * res[:, 1], 'k--', lw=0.9, label=r'$\sim h$ (1-й порядок)')
for n, h, C, E in res: ax.annotate(f'{int(n)}', (h, C), textcoords='offset points', xytext=(3, 5), fontsize=7)
ax.invert_xaxis(); ax.set(xlabel='шаг сетки h', ylabel='погрешность', title='Сходимость по сетке (эталон 1280×1280)')
ax.grid(alpha=.3, which='both'); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f'{F}/conv_grid.png', dpi=160); plt.close(fig)

X, Y, W = load(f'{B}/sol80.txt'); E = np.ma.masked_where(X**2 + (Y - 2)**2 < 1, abs(W - Wr[::16, ::16]))
t = np.linspace(np.pi, 2 * np.pi, 200)
fig, ax = plt.subplots(figsize=(5.6, 4.7))
cs = ax.pcolormesh(X, Y, E, shading='auto', cmap='magma_r'); ax.plot(np.cos(t), 2 + np.sin(t), 'c-', lw=1)
ax.set_aspect('equal'); ax.set(xlabel='x', ylabel='y', title='Погрешность на сетке 80×80')
fig.colorbar(cs, ax=ax, label='|w - w_ref|'); fig.tight_layout(); fig.savefig(f'{F}/error_map.png', dpi=160); plt.close(fig)

# 3. Тест с точным решением: порядок сходимости по сетке
H, C, E = [], [], []
for n in grids:
    out = seq(n, n, 1e-10, '-', '-', 'exact', exe='seq_exact')
    H.append(4 / n); C.append(float(out.split('C=')[1].split()[0])); E.append(float(out.split('E=')[1].split()[0]))
H, C, E = map(np.array, (H, C, E))
print('точное решение: порядок C', np.log2(C[:-1] / C[1:]).round(2), ' E', np.log2(E[:-1] / E[1:]).round(2))
fig, ax = plt.subplots(figsize=(5.6, 4.2))
ax.loglog(H, C, 'o-', label=r'$\|w^h-u\|_C$'); ax.loglog(H, E, 's-', label=r'$\|w^h-u\|_E$')
ax.loglog(H, 0.08 * H**2, 'k--', lw=.9, label=r'$\sim h^2$')
ax.invert_xaxis(); ax.set(xlabel='шаг сетки h', ylabel='ошибка', title='Тест с точным решением u = sin(x+y)')
ax.grid(alpha=.3, which='both'); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f'{F}/exact_grid.png', dpi=160)
