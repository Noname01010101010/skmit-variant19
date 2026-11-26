# Сбор результатов: results/polus/*.out -> report/tables/table{1,2,3}.tex и report/fig/speedup_*.png
# Запуск: python3 collect.py [папка_с_результатами]
import re, glob, os, sys
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
src = sys.argv[1] if len(sys.argv) > 1 else 'results/polus'
ITERS = {'400x600': 1499, '800x1200': 2892}     # число итераций (не зависит от числа процессов)
seq, omp, mpi, hyb = {}, {}, {}, {}
for fn in glob.glob(os.path.join(src, '*.out')):
    for l in open(fn, errors='ignore'):
        m = re.search(r'M=(\d+) N=(\d+) iters=(\d+) time=([\d.]+)', l)
        if not m: continue
        g, it, t = f'{m[1]}x{m[2]}', int(m[3]), float(m[4])
        if l.startswith('SEQ'): d, k = seq, g
        elif l.startswith('OMP'): d, k = omp, (int(re.search(r'threads=(\d+)', l)[1]), g)
        elif 'threads=' in l: d, k = hyb, (int(re.search(r'procs=(\d+)', l)[1]), int(re.search(r'threads=(\d+)', l)[1]), g)
        else: d, k = mpi, (int(re.search(r'procs=(\d+)', l)[1]), g)
        if k not in d or t < d[k][1]: d[k] = (it, t)

def cell(d, k, g):
    if k in d and g in seq:
        it, t = d[k]; return f'{it} & {t:.3f} & {seq[g][1]/t:.2f}'
    if k in d:
        it, t = d[k]; return f'{it} & {t:.3f} & --'
    return f'{ITERS[g]} & & '
os.makedirs('report/tables', exist_ok=True)
def gname(g): return g.replace('x', r' $\times$ ')
rows = [(2,'400x600'),(4,'400x600'),(8,'400x600'),(16,'400x600'),(4,'800x1200'),(8,'800x1200'),(16,'800x1200'),(32,'800x1200')]
open('report/tables/table1.tex','w').write('\n'.join(f'{n} & {gname(g)} & {cell(omp,(n,g),g)} \\\\' for n,g in rows)+'\n')
open('report/tables/table2.tex','w').write('\n'.join(f'{n} & {gname(g)} & {cell(mpi,(n,g),g)} \\\\' for n,g in rows)+'\n')
rows3 = [(2,t,'400x600') for t in (1,2,4,8)] + [(4,t,'800x1200') for t in (1,2,4,8)]
open('report/tables/table3.tex','w').write('\n'.join(f'{p} & {t} & {gname(g)} & {cell(hyb,(p,t,g),g)} \\\\' for p,t,g in rows3)+'\n')

os.makedirs('report/fig', exist_ok=True)
for name, d, rr in [('omp', omp, rows), ('mpi', mpi, rows)]:
    fig, ax = plt.subplots(figsize=(5.5, 4)); ok = False
    for g in ITERS:
        pts = [(n, seq[g][1]/d[(n, g)][1]) for n, gg in rr if gg == g and (n, g) in d and g in seq]
        if pts: ok = True; ax.plot(*zip(*pts), 'o-', label=gname(g).replace('$\\times$', '×'))
    if ok:
        ax.set_xscale('log', base=2); ax.set_xlabel('нити' if name == 'omp' else 'процессы')
        ax.set_ylabel('ускорение'); ax.grid(alpha=.3); ax.legend()
        fig.tight_layout(); fig.savefig(f'report/fig/speedup_{name}.png', dpi=160)
    plt.close(fig)
fig, ax = plt.subplots(figsize=(5.5, 4)); ok = False
for p, g in [(2, '400x600'), (4, '800x1200')]:
    pts = [(t, seq[g][1]/hyb[(p, t, g)][1]) for t in (1, 2, 4, 8) if (p, t, g) in hyb and g in seq]
    if pts: ok = True; ax.plot(*zip(*pts), 'o-', label=f'{p} MPI, ' + gname(g).replace('$\\times$', '×'))
if ok:
    ax.set_xscale('log', base=2); ax.set_xlabel('нитей на процесс'); ax.set_ylabel('ускорение')
    ax.grid(alpha=.3); ax.legend(); fig.tight_layout(); fig.savefig('report/fig/speedup_hybrid.png', dpi=160)
print('seq:', seq); print('omp:', omp); print('mpi:', mpi); print('hyb:', hyb)
