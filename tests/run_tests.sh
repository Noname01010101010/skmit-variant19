#!/bin/bash
# Набор тестов. Запуск из корня репозитория:  make test   (или bash tests/run_tests.sh)
# Переменные окружения:
#   MPIRUN   - команда запуска MPI (по умолчанию "mpirun"; на своей машине может понадобиться
#              MPIRUN="mpirun --oversubscribe")
#   CXX, MPICXX - компиляторы
set -u
cd "$(dirname "$0")/.."
CXX=${CXX:-g++}
MPICXX=${MPICXX:-mpicxx}
MPIRUN=${MPIRUN:-mpirun}
FLAGS="-O1 -std=c++11"
B=tests/build
mkdir -p $B
PASS=0; FAIL=0
ok()   { echo "[ OK ] $*"; PASS=$((PASS+1)); }
bad()  { echo "[FAIL] $*"; FAIL=$((FAIL+1)); }
# максимальная разность третьих столбцов двух файлов решения
maxdiff() { paste "$1" "$2" | awk 'NR>1{d=$3-$6; if(d<0)d=-d; if(d>m)m=d} END{printf "%.3e", m}'; }
lt() { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a<b)}'; }
iters() { sed -n 's/.*iters=\([0-9]*\).*/\1/p' | head -1; }

echo "=== Сборка"
$CXX $FLAGS seq.cpp -o $B/seq || exit 1
$CXX $FLAGS -DTEST_EXACT seq.cpp -o $B/seq_exact || exit 1
$CXX $FLAGS -fopenmp omp.cpp -o $B/omp || exit 1
$MPICXX $FLAGS mpi.cpp -o $B/mpi || exit 1
$MPICXX $FLAGS -fopenmp hybrid.cpp -o $B/hybrid || exit 1
$CXX $FLAGS tests/test_geometry.cpp -o $B/test_geometry || exit 1
$CXX $FLAGS tests/test_partition.cpp -o $B/test_partition || exit 1

echo; echo "=== 1. Площади и коэффициенты схемы"
if $B/test_geometry > $B/geometry.log; then ok "test_geometry ($(grep -c OK $B/geometry.log) проверок)"
else bad "test_geometry"; grep FAIL $B/geometry.log; fi

echo; echo "=== 2. Разбиение на домены"
if $B/test_partition > $B/partition.log; then ok "test_partition: $(tail -1 $B/partition.log)"
else bad "test_partition"; grep FAIL $B/partition.log; fi

echo; echo "=== 3. Тест с точным решением u = sin(x+y), f = 2 sin(x+y): второй порядок"
prev=""
for n in 40 80 160; do
  e=$($B/seq_exact $n $n 1e-10 - - exact | sed -n 's/.*C=\([^ ]*\).*/\1/p')
  echo "    сетка ${n}x${n}: ||w-u||_C в D = $e"
  if [ -n "$prev" ]; then
    p=$(awk -v a="$prev" -v b="$e" 'BEGIN{printf "%.2f", log(a/b)/log(2)}')
    if lt 1.8 "$p" && lt "$p" 2.3; then ok "порядок сходимости $p (ожидается 2)"
    else bad "порядок сходимости $p (ожидается 2)"; fi
  fi
  prev=$e
done

echo; echo "=== 4. Метод сопряженных градиентов"
for n in 20 40 80; do
  $B/seq $n $n 1e-8 - $B/tr_$n.txt > $B/seq_$n.log
  it=$(iters < $B/seq_$n.log)
  # J(w(k)) не возрастает (допуск - ошибки округления 1e-12 относительно |J|)
  bad_j=$(awk 'NR>1 && $3 > p + 1e-12*(p<0?-p:p) {c++} {p=$3} END{print c+0}' $B/tr_$n.txt)
  last=$(tail -1 $B/tr_$n.txt | awk '{print $2}')
  if [ "$bad_j" = 0 ]; then ok "${n}x${n}: функционал J монотонно убывает ($it итераций)"
  else bad "${n}x${n}: J возрастает на $bad_j итерациях"; fi
  if lt "$last" 1e-8; then ok "${n}x${n}: критерий остановки ||w(k+1)-w(k)|| = $last < 1e-8"
  else bad "${n}x${n}: критерий остановки не выполнен ($last)"; fi
  if [ "$it" -le $(( (n-1)*(n-1) )) ]; then ok "${n}x${n}: итераций $it <= числа неизвестных"
  else bad "${n}x${n}: итераций больше числа неизвестных"; fi
done
# регрессия: число итераций на сетках из задания
exp="66 130 256"; got="$(for n in 20 40 80; do iters < $B/seq_$n.log; done | xargs)"
[ "$got" = "$exp" ] && ok "число итераций на сетках 20/40/80: $got" || bad "число итераций: $got (ожидалось $exp)"

echo; echo "=== 5. Совпадение параллельных версий с последовательной"
for g in "80 80" "37 53"; do
  set -- $g; M=$1; N=$2
  $B/seq $M $N 1e-8 $B/ref.txt > $B/ref.log; itref=$(iters < $B/ref.log)
  run() {  # run <метка> <команда...>
    local tag=$1; shift
    "$@" > $B/run.log 2>&1
    local it=$(iters < $B/run.log)
    local d=$(maxdiff $B/ref.txt $B/out.txt)
    if [ "$it" = "$itref" ] && lt "$d" 1e-10; then ok "${M}x${N} $tag: итераций $it, max|w-w_seq| = $d"
    else bad "${M}x${N} $tag: итераций $it (seq $itref), max|w-w_seq| = $d"; tail -3 $B/run.log; fi
  }
  for t in 1 4 16; do run "OpenMP $t нитей" env OMP_NUM_THREADS=$t $B/omp $M $N 1e-8 $B/out.txt; done
  for p in 1 2 4 6; do run "MPI $p проц." env OMP_NUM_THREADS=1 $MPIRUN -np $p $B/mpi $M $N 1e-8 $B/out.txt; done
  for p in 1 2; do run "MPI+OpenMP ${p}x4" env OMP_NUM_THREADS=4 $MPIRUN -np $p $B/hybrid $M $N 1e-8 $B/out.txt; done
done

echo; echo "=== Итог: пройдено $PASS, не пройдено $FAIL"
[ $FAIL = 0 ]
