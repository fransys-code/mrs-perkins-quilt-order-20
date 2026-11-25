# Mrs. Perkins's quilt: all dissections of a square into 19 and into 20 squares

*Mrs. Perkins's quilt* asks for the smallest number f(n) of integer-sided squares, with no common
divisor of their sides, into which an n × n square can be cut
([OEIS A005670](https://oeis.org/A005670)). The *order* of a dissection is its number of squares.
E. Wynn ([arXiv:1308.5420](https://arxiv.org/abs/1308.5420), Discrete Math. 334, 2014) enumerated
all dissections of order at most 18, which proves f(n) for n ≤ 120.

This repository contains an exhaustive enumeration of **orders 19 and 20**, by a different and
much faster method.

| Order | Graphs examined | Dissections (A221841) | Largest side (A089047) |
|---|---|---|---|
| 19 | 10,213,424,233 | **9,392,946** | **126** (the conjectured value) |
| 20 | 57,974,895,671 | **39,823,705** | **158** |

Consequences for f(n):

- **f(n) = 20 for n = 121, 123, 124, 125, for every n from 127 to 154, and for n = 157, 158**
  (no dissection of order ≤ 19 has these sides; order-20 dissections exist). Together with the
  known values, A005670 is determined without a gap up to n = 154.
- **f(155) ≥ 21, f(156) ≥ 21, and f(n) ≥ 21 for every n ≥ 159.**
- f(122) = f(126) = 19 (already given by Wynn's table).

Order 19 also gives A221842(19) = 75,063,024 and A232484(19) = 178,649.

> **Status.** Not checked by a third party; novelty not confirmed. Both orders come from a
> **single implementation** of the enumeration. Exhaustiveness rests on the theorem of Brooks,
> Smith, Stone and Tutte, on the exhaustiveness of plantri, and on full agreement with the
> published results up to order 18 — not on a formal proof of the program or on a second,
> independent enumeration. Corrections and pointers to prior work are welcome: please open an issue.

## Method

1. **Graphs.** `plantri -qc2m2 N+3` (Brinkmann and McKay) generates each simple quadrangulation of
   the sphere with N+3 vertices once. It encodes a non-separable planar map H with N+1 edges and
   its dual.
2. **Electrical network** (Brooks, Smith, Stone, Tutte). Vertices of H are the maximal horizontal
   segments of the dissection, edges are the squares, plus one "battery" edge e = (a, b). The
   squared rectangle obtained by removing e is a square if and only if the effective resistance
   of H between a and b is 1/2.
3. **One matrix inversion per map** gives the N+1 effective resistances. Orders 19 and 20 were run
   entirely in exact integer arithmetic; a floating-point filter run alongside gave no
   disagreement on any graph.
4. **Dissection.** Integer potentials give the coordinates; maps with a zero current are
   discarded; sides are divided by their gcd; the tiling is checked explicitly; duplicates are
   removed with a canonical form under the 8 symmetries of the square.

Validation: for orders 1 to 18 the program reproduces A221841, A221842, A232484 and A089047, the
180 cells of Wynn's Table 1, his Table 3, and A005670 for n ≤ 120, with no discrepancy
(`results/summary_orders_1_to_19.txt`). The numbers of graphs are the terms of A113201. Cost:
about 12 core-hours for order 19 and 74 for order 20.

## Files

- `enumeration/` — `perkins.c` (plantri plug-in), build and run scripts, `bilan.py` (summary and
  comparison with published tables). Comments are in French.
- `verification/check_dissections.py` — independent checker (no code shared with the enumeration).
- `results/sides_order_19.txt`, `results/sides_order_20.txt` — number of dissections per side.
- `results/dissections_order_19_side_92_and_more.txt` (21,944 lines) and
  `results/dissections_order_20_side_121_and_more.txt` (39,771 lines) — the dissections that
  matter for f(n), one per line: `order side x,y,size x,y,size …`.
- `results/sha256_full_lists.txt` — SHA-256 of the complete lists of dissections for orders 1
  to 20 (about 9 GB of text, not stored here; available on request through an issue).

## Verify

```bash
python3 verification/check_dissections.py results/dissections_order_20_side_121_and_more.txt 20 121
python3 verification/check_dissections.py results/dissections_order_19_side_92_and_more.txt 19 92
```

For each line the script checks, in integers, that the squares tile the n × n square exactly
(column sweep, area), that the gcd of the sides is 1 and that the number of squares is the stated
order. This proves the upper bounds f(n) ≤ 20; the lower bounds rest on the exhaustiveness of the
enumeration.

To recompute: get plantri 5.5, then `./compile.sh PLANTRI_DIR`, `PERKINS_EXACT=1 ./lance.sh 19 256 8`
and `./fusion.sh run 19` in `enumeration/`.

*Résumé en français : énumération exhaustive des découpages d'un carré en 19 puis 20 carrés à côtés
entiers premiers entre eux (9 392 946 et 39 823 705). Le nombre minimal de pièces f(n) est connu
sans trou jusqu'à n = 154 : f(n) = 20 pour 121, 123, 124, 125, 127 à 154, 157 et 158 ; au moins 21
pièces pour 155, 156 et tout côté de 159 ou plus. Une seule implémentation de l'énumération ; pas
encore vérifié par un tiers ; nouveauté non confirmée.*

## Licence

Code: MIT (`LICENSE`). Data: CC BY 4.0 (`LICENSE-DATA.md`). Written with the help of Claude.
