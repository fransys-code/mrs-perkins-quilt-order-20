#!/bin/bash
# Usage : ./compile.sh /chemin/vers/plantri55   (dossier contenant plantri.c)
# Produit ./perkins : plantri avec le greffon perkins.c.
set -e
ici="$(dirname "$(readlink -f "$0")")"
src="${1:?dossier des sources de plantri}"
gcc -O3 -march=native -w "-DPLUGIN=\"$ici/perkins.c\"" -o "$ici/perkins" "$src/plantri.c" -lm
