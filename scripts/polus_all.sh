#!/bin/bash
# Полный цикл на IBM Polus с локальной машины (запуск из корня репозитория):
#   bash scripts/polus_all.sh
# Заливает код, компилирует, ставит все задания в очередь, ждет их завершения,
# скачивает результаты в results/polus/ и пересобирает таблицы и графики отчета.
# Логин и адрес можно переопределить: POLUS_LOGIN=... POLUS_HOST=... bash scripts/polus_all.sh
set -euo pipefail
cd "$(dirname "$0")/.."
LOGIN=${POLUS_LOGIN:-edu-cmc-skmodel26-608-02}
HOST=${POLUS_HOST:-polus.hpc.cs.msu.ru}
KEYID=TfmRzD0RLt          # хвост публичного ключа из таблицы доступа

# ---------- 1. Поиск приватного ключа, зарегистрированного на Polus ----------
KEY=""
for pub in ~/.ssh/*.pub ~/Downloads/*.pub ~/Desktop/*.pub ~/Documents/*.pub; do
  [ -f "$pub" ] || continue
  grep -q "$KEYID" "$pub" || continue
  priv="${pub%.pub}"
  if [ -f "$priv" ]; then KEY="$priv"; echo ">>> Ключ Polus: $KEY"; fi
  break
done
if [ -z "$KEY" ]; then
  echo "!!! Не нашел пару к ключу ...$KEYID. Публичные ключи на этом Маке:"
  for pub in ~/.ssh/*.pub; do
    [ -f "$pub" ] && echo "    $pub -> $(awk '{print $1, substr($2,length($2)-9)}' "$pub")"
  done
  exit 1
fi

# ---------- 2. Проверка: ключ с паролем нужно один раз разблокировать ----------
# Если ключ зашифрован, а в ssh-agent его нет и терминала нет - ssh молча
# получит "Permission denied". Проверяем это заранее и объясняем, что делать.
ENCRYPTED=0
ssh-keygen -y -P "" -f "$KEY" >/dev/null 2>&1 || ENCRYPTED=1
FP=$(ssh-keygen -lf "$KEY" 2>/dev/null | awk '{print $2}')
IN_AGENT=0
if [ "$ENCRYPTED" = 1 ] && ssh-add -l 2>/dev/null | grep -qF "$FP"; then IN_AGENT=1; fi
if [ "$ENCRYPTED" = 1 ] && [ "$IN_AGENT" = 0 ] && [ ! -t 0 ]; then
  cat <<MSG
!!! Ключ $KEY защищен паролем, в ssh-agent его нет, а ввести пароль негде
    (скрипт запущен без терминала). Сделайте один раз:

        ssh-add --apple-use-keychain $KEY

    и запустите скрипт снова. Либо запустите его прямо в Терминале,
    тогда пароль спросят здесь же.
MSG
  exit 2
fi

SOCK=/tmp/polus_ctl_$$
SSHOPT="-i $KEY -o IdentitiesOnly=yes -o AddKeysToAgent=yes
        -o ControlMaster=auto -o ControlPath=$SOCK -o ControlPersist=2h
        -o ServerAliveInterval=30 -o ConnectTimeout=30"
# shellcheck disable=SC2086
SSH="ssh $SSHOPT $LOGIN@$HOST"
SCP="scp $SSHOPT"
cleanup() { ssh -o ControlPath=$SOCK -O exit "$LOGIN@$HOST" 2>/dev/null || true; }
trap cleanup EXIT

echo ">>> Подключаюсь к $LOGIN@$HOST (пароль от ключа спросят один раз)"
$SSH true

# ---------- 3. Заливка исходников ----------
echo ">>> Копирую код"
$SSH "rm -rf ~/variant19 && mkdir -p ~/variant19"
$SCP common.h partition.h seq.cpp omp.cpp mpi.cpp hybrid.cpp Makefile scripts/run_polus.sh \
     "$LOGIN@$HOST:variant19/"

# ---------- 4. Сборка и постановка заданий в очередь ----------
echo ">>> Компилирую и ставлю задания в очередь"
$SSH 'bash -lc "set -e; cd ~/variant19
  source /polusfs/setenv/setup.SMPI
  chmod +x run_polus.sh
  make clean >/dev/null 2>&1 || true
  make
  ls -l seq omp mpi hybrid
  ./run_polus.sh"'

# ---------- 5. Ожидание завершения ----------
echo ">>> Жду завершения заданий (проверка раз в минуту)"
STUCK=0
while true; do
  # без bash -lc: иначе MOTD из профиля попадает в stdout и портит число
  n=$($SSH 'bjobs -noheader 2>/dev/null | wc -l' 2>/dev/null | tail -1 | tr -dc 0-9) || n=""
  if [ -z "$n" ]; then
    STUCK=$((STUCK+1))
    echo "    $(date +%H:%M) не смог опросить очередь (попытка $STUCK/5)"
    [ $STUCK -ge 5 ] && { echo "!!! Очередь не опрашивается, прерываю ожидание"; break; }
    sleep 60; continue
  fi
  STUCK=0
  echo "    $(date +%H:%M) в очереди/выполняется: $n"
  [ "$n" = "0" ] && break
  sleep 60
done

# ---------- 6. Скачивание результатов ----------
echo ">>> Скачиваю результаты"
OUT=results/polus
[ -d $OUT ] && mv $OUT "${OUT}_$(date +%m%d_%H%M)"; mkdir -p $OUT
$SCP "$LOGIN@$HOST:variant19/results/*" $OUT/ || echo "!!! results пусты"
$SCP "$LOGIN@$HOST:variant19/sol_*.txt" $OUT/ 2>/dev/null \
  || echo "    (файлы решения sol_*.txt не найдены)"
echo ">>> Скачано $(ls $OUT | wc -l | tr -d ' ') файлов в $OUT"

echo ">>> Непустые .err (если есть):"
find $OUT -name '*.err' -size +0 -exec echo "    --- {}" \; -exec head -5 {} \; || true

echo ">>> Результаты:"
grep -h "iters=" $OUT/*.out 2>/dev/null | sort || echo "    (нет строк с iters=)"

# ---------- 7. Таблицы и графики отчета ----------
echo ">>> Пересобираю таблицы и графики"
python3 scripts/collect.py $OUT
SOL=$(ls $OUT/sol_800x1200.txt $OUT/sol_80.txt 2>/dev/null | head -1)
if [ -n "$SOL" ]; then python3 scripts/plot_solution.py "$SOL"; fi
echo ">>> Готово. Отчет: cd report && xelatex report.tex && xelatex report.tex"
