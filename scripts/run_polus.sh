#!/bin/bash
# Постановка всех расчетов в очередь LSF на IBM Polus (запуск на Polus в каталоге с собранными
# программами). Ключ -n: только создать файлы заданий jobs/*.lsf. REPS - число повторов.
D=1e-8
QUEUE=${QUEUE:-normal}
W=00:15
REPS=${REPS:-3}
SUBMIT=1; [ "$1" = "-n" ] && SUBMIT=0
rm -rf jobs results; mkdir -p jobs results

job() {
  local name=$1 np=$2 nt=$3; shift 3
  local reps=$REPS; case $name in seqsol_*) reps=1;; esac
  local f=jobs/$name.lsf
  local core=$nt;  [ $nt -gt 20 ] && core=20
  local ptile=$((20 / core)); [ $ptile -gt $np ] && ptile=$np
  {
    echo "#BSUB -J $name"
    echo "#BSUB -q $QUEUE"
    echo "#BSUB -n $np"
    echo "#BSUB -x"
    echo "#BSUB -W $W"
    echo "#BSUB -o results/$name.out"
    echo "#BSUB -e results/$name.err"
    echo "#BSUB -R \"span[ptile=$ptile]\""
    echo "#BSUB -R \"affinity[core($core)]\""
    echo "export OMP_NUM_THREADS=$nt"
    echo "for r in \$(seq $reps); do $*; done"
  } > $f
  if [ $SUBMIT = 1 ]; then bsub < $f; fi
}

for g in "80 80" "400 600" "800 1200"; do job seq_${g/ /x} 1 1 ./seq $g $D; done
for n in 20 40 80; do job seqsol_$n 1 1 ./seq $n $n $D sol_$n.txt; done
job seqsol_800x1200 1 1 ./seq 800 1200 $D sol_800x1200.txt

for t in 1 4 16;    do job omp_${t}_80x80    1 $t ./omp 80 80 $D;    done
for t in 2 4 8 16;  do job omp_${t}_400x600  1 $t ./omp 400 600 $D;  done
for t in 4 8 16 32; do job omp_${t}_800x1200 1 $t ./omp 800 1200 $D; done

for p in 1 2 4;     do job mpi_${p}_80x80    $p 1 mpiexec ./mpi 80 80 $D;    done
for p in 2 4 8 16;  do job mpi_${p}_400x600  $p 1 mpiexec ./mpi 400 600 $D;  done
for p in 4 8 16 32; do job mpi_${p}_800x1200 $p 1 mpiexec ./mpi 800 1200 $D; done

for p in 1 2;   do job hyb_${p}x4_80x80     $p 4  mpiexec ./hybrid 80 80 $D;     done
for t in 1 2 4 8; do job hyb_2x${t}_400x600  2 $t mpiexec ./hybrid 400 600 $D;  done
for t in 1 2 4 8; do job hyb_4x${t}_800x1200 4 $t mpiexec ./hybrid 800 1200 $D; done
