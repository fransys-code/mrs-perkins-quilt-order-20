#!/bin/bash
# Usage : ./fusion.sh DOSSIER_RUN ORDRE...
# Pour chaque ordre : additionne les journaux des parts (graphes, solutions, temps CPU)
# et fusionne les découpages en supprimant les doublons entre parts.
set -e
run="${1:?dossier run}"; shift
for N in "$@"; do
    d="$run/N$N"
    parts=$(ls "$d"/part_*.log | wc -l)
    finies=$(grep -l "^PERKINS ordre" "$d"/part_*.log | wc -l)
    cat "$d"/part_*.log | awk -v N="$N" -v parts="$parts" -v finies="$finies" '
        /^PERKINS ordre/ { for (i = 2; i <= NF; i++) { split($i, a, "="); s[a[1]] += a[2] } }
        / cpu=/ { for (i = 1; i <= NF; i++) if ($i ~ /^cpu=/) { split($i, a, "="); cpu += a[2] } }
        END { printf "ordre %d : parts %d/%d graphes %d candidats %d courant_nul %d solutions %d faux_positifs_flottant %d cpu_s %.1f microsecondes_par_graphe %.3f\n",
              N, finies, parts, s["graphes"], s["candidats"], s["courant_nul"], s["solutions"], s["faux_positifs_flottant"], cpu, (s["graphes"] ? 1e6 * cpu / s["graphes"] : 0) }'
    LC_ALL=C sort -u "$d"/part_*.dump > "$run/decoupages_N$N.txt"
    echo "         decoupages distincts : $(wc -l < "$run/decoupages_N$N.txt")"
done
