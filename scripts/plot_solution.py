# Рисунок приближенного решения. Запуск из корня репозитория: python3 scripts/plot_solution.py sol.txt
import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
fn = sys.argv[1] if len(sys.argv) > 1 else 'sol_800x1200.txt'
M, N = map(int, open(fn).readline().split())
d = np.loadtxt(fn, skiprows=1); X, Y, W = (d[:, k].reshape(M+1, N+1) for k in range(3))
Wm = np.ma.masked_where(X**2 + (Y-2)**2 < 1, W)          # только область D
t = np.linspace(np.pi, 2*np.pi, 200)
fig, ax = plt.subplots(figsize=(6.4, 5.4))
cs = ax.contourf(X, Y, Wm, 40, cmap='viridis'); ax.contour(X, Y, Wm, 15, colors='k', linewidths=0.4)
ax.plot(np.cos(t), 2+np.sin(t), 'r-', lw=1.2); ax.set_aspect('equal')
ax.set_xlabel('x'); ax.set_ylabel('y'); fig.colorbar(cs, ax=ax, label='u(x,y)')
ax.set_title(f'Приближенное решение, сетка {M}×{N}')
fig.tight_layout(); fig.savefig('report/fig/solution_map.png', dpi=160)
fig = plt.figure(figsize=(7, 5.4)); ax = fig.add_subplot(projection='3d')
s = max(1, M//100)
ax.plot_surface(X[::s, ::s], Y[::s, ::s], Wm[::s, ::s], cmap='viridis', linewidth=0, antialiased=True)
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('u'); ax.view_init(30, -120)
ax.set_title(f'Приближенное решение, сетка {M}×{N}')
fig.tight_layout(); fig.savefig('report/fig/solution_3d.png', dpi=160)
