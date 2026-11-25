/* Greffon plantri : courtepointe de Mrs Perkins.
 *
 * Compilation :  cc -O3 -march=native '-DPLUGIN="perkins.c"' -o perkins plantri.c -lm
 * Usage :        ./perkins -qc2m2 <N+3> [res/mod] > /dev/null
 *
 * plantri -qc2m2 engendre chaque quadrangulation simple de la sphère une fois
 * (à isomorphisme et miroir près). Une quadrangulation simple à N+3 sommets est
 * la carte d'incidence sommets-faces d'une carte planaire non séparable H à N+1
 * arêtes (et de sa duale) : une classe de couleur = sommets de H, l'autre = faces
 * de H, une face de la quadrangulation = une arête de H.
 *
 * Réseau électrique (Brooks, Smith, Stone, Tutte) : on retire une arête e = (a,b)
 * de H, on branche la pile entre a et b, chaque autre arête est une résistance
 * unité. Les courants donnent un rectangle découpé en N carrés. Ce rectangle est
 * un carré si et seulement si la résistance effective de H entre a et b vaut 1/2
 * (e, résistance unité, en parallèle avec le reste, qui doit valoir 1). On inverse donc le laplacien
 * réduit de H une seule fois par carte et on lit les N+1 résistances.
 *
 * Variables d'environnement :
 *   PERKINS_DUMP=fichier   écrit chaque découpage (forme canonique sous les 8
 *                          symétries du carré), une fois par processus
 *   PERKINS_DUMP_MIN=n     n'écrit que les côtés >= n
 *   PERKINS_EXACT=1        tout en entiers (pas de filtre flottant), et compte
 *                          les désaccords avec le filtre flottant
 */

#include <math.h>
#include <string.h>

#define FILTER perkins_filter
#define SUMMARY perkins_summary
#define PLUGIN_INIT perkins_init();

#define PK_MAXK 15
#define PK_W 16
#define PK_MAXM 64
#define PK_MAXSIDE 8192
#define PK_TOL 1e-6

typedef long long pk_i64;
typedef unsigned long long pk_u64;

static pk_u64 pk_graphs, pk_cand, pk_hits, pk_zero, pk_new, pk_mismatch, pk_falsepos;
static pk_u64 pk_hist_hits[PK_MAXSIDE], pk_hist_new[PK_MAXSIDE];
static int pk_exact, pk_dump_min, pk_maxT_seen;
static FILE *pk_dump;

static pk_u64 *pk_tab;
static size_t pk_cap, pk_fill;

static void pk_die(const char *msg)
{
    fprintf(stderr, ">E perkins: %s\n", msg);
    exit(2);
}

static void perkins_init(void)
{
    const char *s;
    pk_exact = (s = getenv("PERKINS_EXACT")) != NULL && s[0] == '1';
    pk_dump_min = (s = getenv("PERKINS_DUMP_MIN")) != NULL ? atoi(s) : 0;
    if ((s = getenv("PERKINS_DUMP")) != NULL && s[0] != '\0') {
        pk_dump = fopen(s, "w");
        if (pk_dump == NULL) pk_die("impossible d'ouvrir PERKINS_DUMP");
    }
    pk_cap = (size_t)1 << 16;
    pk_tab = (pk_u64 *)calloc(2 * pk_cap, sizeof(pk_u64));
    if (pk_tab == NULL) pk_die("memoire");
}

static pk_u64 pk_mix(pk_u64 h, pk_u64 v)
{
    h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    return h;
}

/* Renvoie 1 si la paire (h1,h2) est nouvelle. */
static int pk_insert(pk_u64 h1, pk_u64 h2)
{
    size_t i, mask;
    if (h1 == 0 && h2 == 0) h2 = 1;
    if (2 * (pk_fill + 1) > pk_cap) {
        size_t ocap = pk_cap, j;
        pk_u64 *old = pk_tab;
        pk_cap *= 2;
        pk_tab = (pk_u64 *)calloc(2 * pk_cap, sizeof(pk_u64));
        if (pk_tab == NULL) pk_die("memoire");
        mask = pk_cap - 1;
        for (j = 0; j < ocap; ++j) {
            if (old[2 * j] == 0 && old[2 * j + 1] == 0) continue;
            i = (size_t)old[2 * j] & mask;
            while (pk_tab[2 * i] != 0 || pk_tab[2 * i + 1] != 0) i = (i + 1) & mask;
            pk_tab[2 * i] = old[2 * j];
            pk_tab[2 * i + 1] = old[2 * j + 1];
        }
        free(old);
    }
    mask = pk_cap - 1;
    i = (size_t)h1 & mask;
    while (pk_tab[2 * i] != 0 || pk_tab[2 * i + 1] != 0) {
        if (pk_tab[2 * i] == h1 && pk_tab[2 * i + 1] == h2) return 0;
        i = (i + 1) & mask;
    }
    pk_tab[2 * i] = h1;
    pk_tab[2 * i + 1] = h2;
    ++pk_fill;
    return 1;
}

static pk_i64 pk_gcd(pk_i64 a, pk_i64 b)
{
    while (b != 0) { pk_i64 t = a % b; a = b; b = t; }
    return a < 0 ? -a : a;
}

/* Adjointe et déterminant du laplacien réduit, Gauss-Jordan sans fraction
   (Bareiss) : à la fin la partie gauche vaut det.I et la droite l'adjointe. */
static pk_i64 pk_adjugate(int k, int m, const int *ea, const int *eb,
                          pk_i64 adj[PK_MAXK + 1][PK_MAXK + 1])
{
    static pk_i64 W[PK_MAXK][2 * PK_MAXK];
    pk_i64 prev = 1, piv, f;
    int i, j, p;

    for (i = 0; i < k; ++i)
        for (j = 0; j < 2 * k; ++j) W[i][j] = (j == k + i);
    for (i = 0; i < m; ++i) {
        int a = ea[i], b = eb[i];
        if (a < k) ++W[a][a];
        if (b < k) ++W[b][b];
        if (a < k && b < k) { --W[a][b]; --W[b][a]; }
    }
    for (p = 0; p < k; ++p) {
        piv = W[p][p];
        if (piv == 0) pk_die("pivot nul");
        for (i = 0; i < k; ++i) {
            if (i == p) continue;
            f = W[i][p];
            for (j = 0; j < 2 * k; ++j)
                W[i][j] = (piv * W[i][j] - f * W[p][j]) / prev;
        }
        prev = piv;
    }
    for (i = 0; i <= k; ++i)
        for (j = 0; j <= k; ++j)
            adj[i][j] = (i < k && j < k) ? W[i][k + j] : 0;
    return prev;
}

/* Construit le découpage de l'arête pile e, le vérifie, l'enregistre.
   Renvoie 1 si c'est un découpage valide, 0 si un courant est nul. */
static int pk_record(int k, int nW, int m, const int *ea, const int *eb,
                     const int *el, const int *er, int e, pk_i64 T,
                     pk_i64 adj[PK_MAXK + 1][PK_MAXK + 1])
{
    pk_i64 phi[PK_MAXK + 1], psi[PK_MAXM], cur[PK_MAXM];
    pk_i64 X[PK_MAXM], Y[PK_MAXM], S[PK_MAXM];
    pk_i64 side, g, x0, y0, xlo, xhi, area;
    pk_u64 key[PK_MAXM], best[PK_MAXM];
    int known[PK_MAXM];
    int i, j, t, n, N, changed, left;

    for (i = 0; i <= k; ++i) phi[i] = adj[i][ea[e]] - adj[i][eb[e]];

    /* contrôle de Kirchhoff : flux net nul hors des pôles */
    {
        pk_i64 net[PK_MAXK + 1];
        for (i = 0; i <= k; ++i) net[i] = 0;
        for (i = 0; i < m; ++i) {
            pk_i64 c = phi[ea[i]] - phi[eb[i]];
            net[ea[i]] += c;
            net[eb[i]] -= c;
        }
        for (i = 0; i <= k; ++i) {
            pk_i64 want = (i == ea[e]) ? T : (i == eb[e]) ? -T : 0;
            if (net[i] != want) pk_die("Kirchhoff viole");
        }
    }

    for (i = 0; i < m; ++i) {
        cur[i] = phi[ea[i]] - phi[eb[i]];
        if (cur[i] == 0) return 0;
    }

    /* la pile referme le circuit : l'arete e porte le courant de retour */
    cur[e] -= T;

    /* potentiels duaux : psi(droite) - psi(gauche) = courant */
    for (i = 0; i < nW; ++i) known[i] = 0;
    known[el[0]] = 1;
    psi[el[0]] = 0;
    left = nW - 1;
    do {
        changed = 0;
        for (i = 0; i < m; ++i) {
            if (known[el[i]] && !known[er[i]]) {
                psi[er[i]] = psi[el[i]] + cur[i]; known[er[i]] = 1; changed = 1; --left;
            } else if (known[er[i]] && !known[el[i]]) {
                psi[el[i]] = psi[er[i]] - cur[i]; known[el[i]] = 1; changed = 1; --left;
            }
        }
    } while (changed && left > 0);
    if (left != 0) pk_die("dual non connexe");
    for (i = 0; i < m; ++i)
        if (psi[er[i]] - psi[el[i]] != cur[i]) pk_die("potentiels duaux incoherents");

    side = phi[ea[e]] - phi[eb[e]];
    if (2 * side != T) pk_die("resistance differente de 1/2");
    y0 = phi[eb[e]];
    xlo = xhi = psi[0];
    for (i = 1; i < nW; ++i) {
        if (psi[i] < xlo) xlo = psi[i];
        if (psi[i] > xhi) xhi = psi[i];
    }
    if (xhi - xlo != side) pk_die("largeur differente de la hauteur");
    x0 = xlo;

    N = 0;
    g = 0;
    area = 0;
    for (i = 0; i < m; ++i) {
        if (i == e) continue;
        S[N] = cur[i] < 0 ? -cur[i] : cur[i];
        X[N] = (psi[el[i]] < psi[er[i]] ? psi[el[i]] : psi[er[i]]) - x0;
        Y[N] = (phi[ea[i]] < phi[eb[i]] ? phi[ea[i]] : phi[eb[i]]) - y0;
        if (X[N] < 0 || Y[N] < 0 || X[N] + S[N] > side || Y[N] + S[N] > side)
            pk_die("carre hors du cadre");
        g = pk_gcd(g, S[N]);
        area += S[N] * S[N];
        ++N;
    }
    if (area != side * side) pk_die("aire incorrecte");
    for (i = 0; i < N; ++i)
        for (j = i + 1; j < N; ++j)
            if (X[i] < X[j] + S[j] && X[j] < X[i] + S[i] &&
                Y[i] < Y[j] + S[j] && Y[j] < Y[i] + S[i])
                pk_die("chevauchement");

    if (side % g != 0) pk_die("pgcd");
    side /= g;
    if (side >= PK_MAXSIDE) pk_die("cote trop grand pour l'histogramme");
    n = (int)side;
    for (i = 0; i < N; ++i) {
        if (X[i] % g != 0 || Y[i] % g != 0) pk_die("coordonnees non multiples du pgcd");
        X[i] /= g; Y[i] /= g; S[i] /= g;
    }

    for (t = 0; t < 8; ++t) {
        int smaller = 0;
        for (i = 0; i < N; ++i) {
            pk_i64 x = X[i], y = Y[i], s = S[i], u;
            if (t & 1) x = n - x - s;
            if (t & 2) y = n - y - s;
            if (t & 4) { u = x; x = y; y = u; }
            key[i] = ((pk_u64)x << 40) | ((pk_u64)y << 20) | (pk_u64)s;
        }
        for (i = 1; i < N; ++i) {
            pk_u64 v = key[i];
            for (j = i; j > 0 && key[j - 1] > v; --j) key[j] = key[j - 1];
            key[j] = v;
        }
        if (t == 0) smaller = 1;
        else
            for (i = 0; i < N; ++i)
                if (key[i] != best[i]) { smaller = key[i] < best[i]; break; }
        if (smaller)
            for (i = 0; i < N; ++i) best[i] = key[i];
    }

    ++pk_hits;
    ++pk_hist_hits[n];
    {
        pk_u64 h1 = 0x243f6a8885a308d3ULL ^ (pk_u64)n, h2 = 0x13198a2e03707344ULL + (pk_u64)N;
        for (i = 0; i < N; ++i) {
            h1 = pk_mix(h1, best[i]);
            h2 = pk_mix(h2, best[i] * 0xc2b2ae3d27d4eb4fULL + (pk_u64)i);
        }
        if (pk_insert(h1, h2)) {
            ++pk_new;
            ++pk_hist_new[n];
            if (pk_dump != NULL && n >= pk_dump_min) {
                fprintf(pk_dump, "%d %d", N, n);
                for (i = 0; i < N; ++i)
                    fprintf(pk_dump, " %d,%d,%d", (int)(best[i] >> 40),
                            (int)((best[i] >> 20) & 0xfffff), (int)(best[i] & 0xfffff));
                fputc('\n', pk_dump);
            }
        }
    }
    return 1;
}

static int perkins_filter(int nbtot, int nbop, int doflip)
{
    int idx[MAXN];
    int ea[PK_MAXM], eb[PK_MAXM], el[PK_MAXM], er[PK_MAXM];
    int nV, nW, k, m, i, j, p, have_adj;
    pk_u64 nbm[MAXN], B, Wm, S, rest;
    double M[PK_W][PK_W] __attribute__((aligned(64)));
    static pk_i64 adj[PK_MAXK + 1][PK_MAXK + 1];
    pk_i64 T = 0;
    EDGE *e;

    (void)nbtot; (void)nbop; (void)doflip;
    ++pk_graphs;

    /* Bicoloration par masques de bits (peu de branchements imprevisibles). */
    for (i = 0; i < nv; ++i) {
        pk_u64 a = 0;
        e = firstedge[i];
        for (j = 0; j < degree[i]; ++j, e = e->next) a |= 1ULL << e->end;
        nbm[i] = a;
    }
    B = 1;
    for (;;) {
        pk_u64 w = 0, b = B;
        for (i = 0; i < nv; ++i) w |= nbm[i] & -((B >> i) & 1);
        for (i = 0; i < nv; ++i) b |= nbm[i] & -((w >> i) & 1);
        Wm = w;
        if (b == B) break;
        B = b;
    }
    nV = __builtin_popcountll(B);
    nW = __builtin_popcountll(Wm);
    if (nV + nW != nv || (B & Wm) != 0) pk_die("bicoloration impossible");
    S = B;
    if (nW < nV) { S = Wm; i = nV; nV = nW; nW = i; }
    for (i = 0; i < nv; ++i) {
        pk_u64 below = (1ULL << i) - 1;
        idx[i] = __builtin_popcountll((((S >> i) & 1) ? S : ~S) & below);
    }
    k = nV - 1;
    if (k > PK_MAXK || nv - 2 >= PK_MAXM) pk_die("graphe trop grand");

    /* Une face (i, w1, u, w2) de la quadrangulation = une arete i-u de H,
       bordee par les faces w1 et w2 de H ; gardee une fois (i < u). */
    m = 0;
    for (rest = S; rest != 0; rest &= rest - 1) {
        i = __builtin_ctzll(rest);
        e = firstedge[i];
        for (j = 0; j < degree[i]; ++j, e = e->next) {
            int u = e->invers->next->end;
            ea[m] = idx[i];
            eb[m] = idx[u];
            el[m] = idx[e->end];
            er[m] = idx[e->prev->end];
            m += (i < u);
        }
    }
    if (m != nv - 2) pk_die("nombre de faces inattendu");

    /* Largeur fixe PK_W : les boucles internes sont de longueur constante et le
       compilateur les vectorise. Les lignes et colonnes au-dela de k restent nulles. */
    memset(M, 0, sizeof(M));
    for (i = 0; i < m; ++i) {
        int a = ea[i], b = eb[i];
        M[a][a] += 1.0;
        M[b][b] += 1.0;
        M[a][b] -= 1.0;
        M[b][a] -= 1.0;
    }
    for (i = 0; i < PK_W; ++i) M[i][k] = M[k][i] = 0.0;
    /* Inversion sur place (matrice symetrique definie positive, pas de pivotage).
       Aucune ecriture scalaire dans une ligne juste avant sa lecture vectorielle :
       la colonne du pivot est traitee par le terme 1 + d de rowq. */
    for (p = 0; p < k; ++p) {
        double d = 1.0 / M[p][p];
        double rowp[PK_W] __attribute__((aligned(64)));
        double rowq[PK_W] __attribute__((aligned(64)));
        double colp[PK_W];
        for (j = 0; j < PK_W; ++j) rowp[j] = M[p][j] * d;
        for (j = 0; j < PK_W; ++j) rowq[j] = rowp[j];
        rowp[p] = d;
        rowq[p] = d + 1.0;
        for (i = 0; i < k; ++i) colp[i] = M[i][p];
        for (i = 0; i < k; ++i) {
            double f = colp[i];
            for (j = 0; j < PK_W; ++j) M[i][j] -= f * rowq[j];
        }
        for (j = 0; j < PK_W; ++j) M[p][j] = rowp[j];
    }

    have_adj = 0;
    if (pk_exact) {
        T = pk_adjugate(k, m, ea, eb, adj);
        have_adj = 1;
    }
    for (i = 0; i < m; ++i) {
        int a = ea[i], b = eb[i];
        double R = M[a][a] + M[b][b] - 2.0 * M[a][b];
        int fl = fabs(R - 0.5) < PK_TOL;
        int ex;
        if (!pk_exact && !fl) continue;
        if (!have_adj) {
            T = pk_adjugate(k, m, ea, eb, adj);
            have_adj = 1;
        }
        ex = 2 * (adj[a][a] + adj[b][b] - 2 * adj[a][b]) == T;
        if (ex != fl) {
            if (pk_exact) ++pk_mismatch;
            else ++pk_falsepos;
        }
        if (!ex) continue;
        ++pk_cand;
        if (T > pk_maxT_seen) pk_maxT_seen = (int)T;
        if (!pk_record(k, nW, m, ea, eb, el, er, i, T, adj)) ++pk_zero;
    }
    return 0;
}

static void perkins_summary(void)
{
    int n;
    fprintf(msgfile, "PERKINS ordre=%d graphes=%llu candidats=%llu courant_nul=%llu "
            "solutions=%llu decoupages_distincts_processus=%llu mode=%s "
            "desaccords_flottant=%llu faux_positifs_flottant=%llu maxT=%d\n",
            maxnv - 3, pk_graphs, pk_cand, pk_zero, pk_hits, pk_new,
            pk_exact ? "exact" : "flottant+exact", pk_mismatch, pk_falsepos, pk_maxT_seen);
    for (n = 0; n < PK_MAXSIDE; ++n)
        if (pk_hist_hits[n] != 0)
            fprintf(msgfile, "PERKINS_N ordre=%d cote=%d solutions=%llu distincts=%llu\n",
                    maxnv - 3, n, pk_hist_hits[n], pk_hist_new[n]);
    if (pk_dump != NULL) fclose(pk_dump);
}
