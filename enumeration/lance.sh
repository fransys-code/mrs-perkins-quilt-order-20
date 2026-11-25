#!/bin/bash
# Usage : ./lance.sh ORDRE PARTS COEURS [DOSSIER_SORTIE] [COTE_MIN_ECRIT]
# Traite l'ordre ORDRE (nombre de carrés) en PARTS parts (res/mod de plantri),
# COEURS parts à la fois, en priorité basse. Reprise : une part déjà finie
# (fichier .log contenant la ligne PERKINS) n'est pas relancée.
# Variables : PERKINS_EXACT=1 (tout en entiers), PERKINS_NICE (priorité, 10 par défaut).
set -e
ici="$(dirname "$(readlink -f "$0")")"
N="${1:?ordre}"; M="${2:?parts}"; J="${3:?coeurs}"
out="${4:-$ici/run/N$N}"; dmin="${5:-0}"
mkdir -p "$out"
export ici N M out dmin
une_part() {
    r="$1"
    log="$out/part_$r.log"
    grep -q "^PERKINS ordre" "$log" 2>/dev/null && return 0
    PERKINS_DUMP="$out/part_$r.dump" PERKINS_DUMP_MIN="$dmin" \
        nice -n "${PERKINS_NICE:-10}" "$ici/perkins" -qc2m2 "$((N + 3))" "$r/$M" > /dev/null 2> "$log"
}
export -f une_part
seq 0 "$((M - 1))" | xargs -P "$J" -I{} bash -c 'une_part {}'
echo "ordre $N : $M parts terminees dans $out"
