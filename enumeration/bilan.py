#!/usr/bin/env python3
"""Bilan et validation des découpages produits par perkins (greffon plantri).

Usage : python3 bilan.py decoupages_N*.txt
Chaque fichier est la fusion dédoublonnée (sort -u) des fichiers .dump d'un ordre :
une ligne « N n x,y,s x,y,s ... » par découpage, sous forme canonique.

Le script revérifie les découpages indépendamment du programme C (pavage exact du
carré, pgcd 1, forme canonique), puis compare aux valeurs publiées :
A221841, A221842, A232484, A089047, A005670 et les tables 1 et 3 de Wynn
(arXiv:1308.5420).
"""
import sys
from collections import defaultdict
from math import gcd

A221841 = [1, 0, 0, 1, 0, 1, 2, 6, 16, 56, 183, 657, 2277, 8813, 34178, 137578, 558734, 2285694]
A221842 = [1, 0, 0, 1, 0, 4, 8, 36, 105, 384, 1340, 4975, 17676, 69052, 270716, 1093218, 4455047, 18246018]
A232484 = [1, 0, 0, 1, 0, 1, 1, 2, 4, 7, 18, 40, 119, 323, 1100, 3594, 13068, 47444]
A089047 = [1, 1, 1, 2, 2, 3, 4, 5, 7, 9, 13, 17, 23, 29, 41, 53, 70, 91]
A005670 = [1, 4, 6, 7, 8, 9, 9, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13,
           13, 13, 13, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16,
           15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17,
           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
           18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18] + [19] * 29

# Table 1 de Wynn : nombre de découpages premiers d'ordre N (7..18) du carré de côté n (4..18).
WYNN_T1 = {
    4: [2, 1, 0, 4, 0, 0, 3, 0, 0, 1, 0, 0],
    5: [0, 5, 0, 1, 10, 0, 14, 6, 0, 20, 3, 0],
    6: [0, 0, 11, 0, 9, 28, 12, 47, 58, 85, 27, 151],
    7: [0, 0, 4, 14, 0, 42, 40, 83, 175, 103, 509, 620],
    8: [0, 0, 0, 28, 17, 13, 99, 156, 325, 492, 1291, 1135],
    9: [0, 0, 0, 9, 28, 25, 129, 142, 336, 1357, 1254, 3256],
    10: [0, 0, 0, 0, 75, 13, 86, 474, 371, 1910, 2701, 7150],
    11: [0, 0, 0, 0, 34, 67, 30, 342, 709, 1341, 4377, 9712],
    12: [0, 0, 0, 0, 9, 176, 97, 327, 1617, 1850, 6622, 15402],
    13: [0, 0, 0, 0, 1, 162, 24, 219, 1134, 1759, 7600, 13377],
    14: [0, 0, 0, 0, 0, 93, 244, 260, 1634, 2926, 9271, 21419],
    15: [0, 0, 0, 0, 0, 27, 432, 244, 781, 4853, 7558, 21220],
    16: [0, 0, 0, 0, 0, 9, 493, 393, 852, 5608, 9943, 32117],
    17: [0, 0, 0, 0, 0, 2, 216, 309, 594, 4101, 6076, 27442],
    18: [0, 0, 0, 0, 0, 0, 187, 1135, 1045, 5757, 14362, 36856],
}
# Table 3 de Wynn : ordre N -> (g(N), nombres de découpages de côté g-9 .. g).
WYNN_T3 = {
    12: (17, [13, 25, 13, 67, 176, 162, 93, 27, 9, 2]),
    13: (23, [244, 432, 493, 216, 187, 113, 44, 11, 1, 2]),
    14: (29, [1231, 798, 696, 561, 288, 152, 53, 35, 17, 4]),
    15: (41, [386, 210, 104, 65, 35, 17, 2, 5, 0, 1]),
    16: (53, [231, 127, 81, 62, 27, 16, 6, 3, 1, 1]),
    17: (70, [131, 57, 42, 26, 11, 10, 4, 4, 1, 3]),
    18: (91, [71, 35, 33, 16, 14, 12, 1, 3, 3, 1]),
}

VERIFIE_TOUT_JUSQUA = 400000


def images(n, sq):
    out = []
    for t in range(8):
        im = []
        for x, y, s in sq:
            if t & 1:
                x = n - x - s
            if t & 2:
                y = n - y - s
            if t & 4:
                x, y = y, x
            im.append((x, y, s))
        im.sort()
        out.append(tuple(im))
    return out


def verifie(N, n, sq):
    assert len(sq) == N
    g = 0
    grille = bytearray(n * n)
    for x, y, s in sq:
        assert s >= 1 and 0 <= x and 0 <= y and x + s <= n and y + s <= n
        g = gcd(g, s)
        for yy in range(y, y + s):
            base = yy * n + x
            assert not any(grille[base:base + s]), "chevauchement"
            grille[base:base + s] = b"\x01" * s
    assert all(grille), "trou"
    assert g == 1, "non premier"
    assert min(images(n, sq)) == tuple(sq), "forme non canonique"


def main(fichiers):
    compte = defaultdict(int)
    orbites = defaultdict(int)
    par_cote = defaultdict(int)
    tailles = defaultdict(set)
    verifies = defaultdict(int)
    for f in fichiers:
        lignes = sum(1 for _ in open(f))
        pas = 1 if lignes <= VERIFIE_TOUT_JUSQUA else lignes // VERIFIE_TOUT_JUSQUA + 1
        precedent = None
        for i, ligne in enumerate(open(f)):
            assert ligne != precedent, "doublon : le fichier doit etre passe par sort -u"
            precedent = ligne
            ch = ligne.split()
            N, n = int(ch[0]), int(ch[1])
            sq = [tuple(map(int, c.split(","))) for c in ch[2:]]
            if i % pas == 0:
                verifie(N, n, sq)
                verifies[N] += 1
            compte[N] += 1
            orbites[N] += len(set(images(n, sq)))
            par_cote[(N, n)] += 1
            tailles[N].add(tuple(sorted(s for _, _, s in sq)))

    ordres = sorted(compte)
    Nmax = max(ordres)
    complet = all(N in compte or A221841[N - 1] == 0 for N in range(1, Nmax + 1))
    ok = True

    def cmp(nom, a, b):
        nonlocal ok
        if a != b:
            ok = False
        return f"{a}" + ("" if a == b else f"  *** ECART : publie {b} ({nom})")

    print("ordre  decoupages(A221841)  toutes orientations(A221842)  collections de tailles(A232484)  plus grand cote(A089047)  reverifies")
    gmax = 0
    for N in range(1, Nmax + 1):
        if N not in compte and A221841[N - 1] != 0:
            print(f"{N:5d}  (ordre non fourni)")
            continue
        g = max((n for (NN, n) in par_cote if NN == N), default=0)
        gmax = max(gmax, g)
        pub = lambda t: t[N - 1] if N <= len(t) else None
        print(f"{N:5d}  {cmp('A221841', compte[N], pub(A221841))}  {cmp('A221842', orbites[N], pub(A221842))}"
              f"  {cmp('A232484', len(tailles[N]), pub(A232484))}  {cmp('A089047', gmax, pub(A089047))}  {verifies[N]}")

    print("\nTable 1 de Wynn (cotes 4 a 18, ordres 7 a 18) :")
    ecarts = cases = 0
    for n, rang in WYNN_T1.items():
        for N, v in zip(range(7, 19), rang):
            if N in compte or (N <= Nmax and A221841[N - 1] == 0):
                cases += 1
                if par_cote.get((N, n), 0) != v:
                    ecarts += 1
                    print(f"  *** ECART cote {n} ordre {N} : trouve {par_cote.get((N, n), 0)}, publie {v}")
    print(f"  {cases} cases comparees, {ecarts} ecart(s)")
    ok = ok and ecarts == 0

    print("\nTable 3 de Wynn (cotes g(N)-9 a g(N)) :")
    for N, (g, rang) in WYNN_T3.items():
        if N not in compte:
            continue
        trouve = [par_cote.get((N, g - 9 + i), 0) for i in range(10)]
        audela = sum(v for (NN, n), v in par_cote.items() if NN == N and n > g)
        bon = trouve == rang and audela == 0
        ok = ok and bon
        print(f"  ordre {N} g={g} : {trouve} " + ("identique" if bon else f"*** ECART : publie {rang}, au-dela de g : {audela}"))

    print("\nf(n) = plus petit ordre (A005670) :")
    f = {}
    for (N, n) in par_cote:
        f[n] = min(f.get(n, 99), N)
    if complet:
        ecarts = 0
        dernier = 0
        for n in range(1, len(A005670) + 1):
            if n in f:
                if f[n] != A005670[n - 1]:
                    ecarts += 1
                    print(f"  *** ECART f({n}) : trouve {f[n]}, publie {A005670[n - 1]}")
            elif A005670[n - 1] <= Nmax:
                ecarts += 1
                print(f"  *** ECART f({n}) : aucun decoupage d'ordre <= {Nmax}, publie {A005670[n - 1]}")
        certifies = sorted(f)
        while dernier + 1 in f:
            dernier += 1
        print(f"  ordres 1 a {Nmax} complets : f(n) certifie pour {len(certifies)} cotes ; "
              f"sans trou jusqu'a n = {dernier} ; au-dela : {[n for n in certifies if n > dernier]}")
        print(f"  tout autre cote demande au moins {Nmax + 1} carres")
        print(f"  comparaison a A005670 (n <= {len(A005670)}) : {ecarts} ecart(s)")
        ok = ok and ecarts == 0
        print("  f(n) trouve :", " ".join(f"{n}:{f[n]}" for n in certifies))
    else:
        print("  ordres incomplets : pas de conclusion sur f(n)")
    print("\nVALIDATION :", "tout concorde" if ok else "ECARTS (voir ci-dessus)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
