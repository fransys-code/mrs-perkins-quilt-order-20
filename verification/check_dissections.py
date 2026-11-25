#!/usr/bin/env python3
"""Contre-vérification indépendante des découpages annoncés (aucun code commun avec
perkins.c ni bilan.py).

Usage : python3 contre_verif.py FICHIER ORDRE COTE_MIN
Lit les lignes « N n x,y,c x,y,c ... » de FICHIER et, pour chaque ligne de côté
n >= COTE_MIN, vérifie exactement, en entiers :
  - le nombre de carrés vaut ORDRE ;
  - chaque carré a un côté >= 1 et des coordonnées entières dans [0, n] ;
  - balayage par colonnes : dans chaque colonne unité [x, x+1], les carrés qui la
    traversent, triés par ordonnée, s'enchaînent exactement de 0 à n (ni trou ni
    chevauchement) ;
  - la somme des aires vaut n² ;
  - le pgcd des côtés vaut 1.
Affiche, par côté, le nombre de découpages vérifiés, et écrit un exemple par côté.
"""
import sys
from functools import reduce
from math import gcd


def controle(ordre, ligne):
    champs = ligne.split()
    nb, cote = int(champs[0]), int(champs[1])
    carres = []
    for c in champs[2:]:
        x, y, s = c.split(",")
        carres.append((int(x), int(y), int(s)))
    if nb != ordre or len(carres) != ordre:
        raise ValueError("nombre de carres")
    for x, y, s in carres:
        if s < 1 or x < 0 or y < 0 or x + s > cote or y + s > cote:
            raise ValueError("carre hors cadre")
    if sum(s * s for _, _, s in carres) != cote * cote:
        raise ValueError("aire")
    for col in range(cote):
        haut = 0
        for y, s in sorted((y, s) for x, y, s in carres if x <= col < x + s):
            if y != haut:
                raise ValueError(f"trou ou chevauchement, colonne {col}, ordonnee {haut}")
            haut = y + s
        if haut != cote:
            raise ValueError(f"colonne {col} incomplete")
    if reduce(gcd, (s for _, _, s in carres)) != 1:
        raise ValueError("pgcd different de 1")
    return cote


def main():
    fichier, ordre, cote_min = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    compte, exemple, vus = {}, {}, set()
    for ligne in open(fichier):
        champs = ligne.split(None, 2)
        if len(champs) < 3 or int(champs[1]) < cote_min:
            continue
        if ligne in vus:
            raise ValueError("ligne en double")
        vus.add(ligne)
        n = controle(ordre, ligne)
        compte[n] = compte.get(n, 0) + 1
        exemple.setdefault(n, ligne.strip())
    print(f"ordre {ordre}, cotes >= {cote_min} : {sum(compte.values())} decoupages verifies, tous valides")
    for n in sorted(compte):
        print(f"cote {n} : {compte[n]} decoupage(s) ; exemple : {exemple[n]}")
    print("cotes atteints :", " ".join(map(str, sorted(compte))))


if __name__ == "__main__":
    main()
