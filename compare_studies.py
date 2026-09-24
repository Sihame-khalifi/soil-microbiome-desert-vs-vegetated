"""Compare soil bacterial communities from two published MGnify studies.

Study A  MGYS00002018  Alaskan tundra and Oklahoma prairie (Johnston et al., 2016)
Study B  MGYS00000770  Atacama Desert, six sites (Crits-Christoph et al., 2013)

Input:  data/study_A_alaska_oklahoma.tsv and data/study_B_atacama.tsv
        (the taxonomy_abundances_SSU tables from MGnify, pipeline v4.0 and v5.0)
Output: figures in results/figures/, tables in results/tables/

Run with:  python compare_studies.py   (on Windows: py -3 compare_studies.py)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu

MIN_READS = 1000        # samples below this are dropped
RAREFY_DEPTH = 1000     # every sample is subsampled to this many reads
N_PERMUTATIONS = 999    # for the PERMANOVA p-value
SEED = 42               # fixed so that every run gives the same numbers

LABEL_A = "Alaska / Oklahoma"
LABEL_B = "Atacama"
GROUP_ORDER = [LABEL_A, LABEL_B]
GROUP_COLORS = {LABEL_A: "#2E6F95", LABEL_B: "#C8663B"}
PHYLUM_PALETTE = ["#264653", "#2A9D8F", "#8AB17D", "#E9C46A", "#F4A261",
                  "#E76F51", "#6D597A", "#B5838D", "#9C6644", "#577590"]

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
FIG_DIR = HERE / "results" / "figures"
TABLE_DIR = HERE / "results" / "tables"

plt.rcParams.update({"font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13,
                     "axes.spines.top": False, "axes.spines.right": False})


# --- Loading and cleaning ---------------------------------------------------

def load_table(path, label):
    df = pd.read_csv(path, sep="\t").rename(columns={"#SampleID": "taxonomy"})
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")]  # trailing empty column
    df = df.set_index("taxonomy")
    print(f"{label}: {df.shape[0]} taxa, {df.shape[1]} columns")
    return df


def to_family_level(lineage):
    """Cut a lineage after the family rank.

    Pipeline v4.0 often goes down to species while v5.0 mostly stops at family,
    so both tables are brought to the same depth before comparing them.
    """
    kept = [t.strip() for t in str(lineage).split(";")
            if t.strip().split("__", 1)[0] in ("sk", "k", "p", "c", "o", "f")]
    while kept and kept[-1].endswith("__"):  # drop empty trailing ranks
        kept.pop()
    return ";".join(kept)


def harmonize(df, label):
    # Study A also holds some eukaryotic 18S reads; keep prokaryotes only.
    prokaryote = df.index.str.startswith("sk__Bacteria") | df.index.str.startswith("sk__Archaea")
    df = df.loc[prokaryote]
    df = df.groupby(df.index.map(to_family_level)).sum()
    print(f"{label}: {df.shape[0]} taxa after harmonization")
    return df


def quality_filter(df, label):
    """Drop assembly entries (ERZ, only 1-2 reads each) and shallow samples."""
    n_columns = df.shape[1]
    erz = [c for c in df.columns if c.startswith("ERZ")]
    df = df.drop(columns=erz)
    depth = df.sum(axis=0)
    shallow = depth[depth < MIN_READS].index.tolist()
    df = df.drop(columns=shallow)
    print(f"{label}: removed {len(erz)} ERZ entries and {len(shallow)} samples "
          f"under {MIN_READS} reads, kept {df.shape[1]}")
    summary = {"study": label, "columns": n_columns, "erz_removed": len(erz),
               "shallow_removed": len(shallow), "samples_kept": df.shape[1]}
    return df, summary


def rarefy(df, depth, rng):
    """Subsample each sample to `depth` reads without replacement.

    Study A has ~72,000 reads per sample and study B ~2,600. Without this step
    the deeper samples would look more diverse for purely technical reasons.
    """
    out = {}
    n_taxa = df.shape[0]
    for sample in df.columns:
        pool = np.repeat(np.arange(n_taxa), df[sample].to_numpy().astype(int))
        picked = rng.choice(pool, size=depth, replace=False)
        out[sample] = np.bincount(picked, minlength=n_taxa)
    return pd.DataFrame(out, index=df.index)


def get_rank(lineage, prefix):
    """Return the name at a given rank, e.g. prefix 'p' for phylum."""
    for token in str(lineage).split(";"):
        token = token.strip()
        if token.startswith(prefix + "__"):
            return token.split("__", 1)[1] or None
    return None


# --- Statistics ---------------------------------------------------------------

def shannon(p):
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())


def pcoa(dist):
    """Classical principal coordinates analysis (first two axes)."""
    n = dist.shape[0]
    centering = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * centering @ (dist ** 2) @ centering
    eigval, eigvec = np.linalg.eigh(b)
    order = np.argsort(eigval)[::-1]
    eigval, eigvec = eigval[order], eigvec[:, order]
    positive = np.clip(eigval, 0, None)
    coords = eigvec[:, :2] * np.sqrt(positive[:2])
    return coords, positive / positive.sum() * 100


def permanova(dist, labels, n_perm, rng):
    """PERMANOVA (Anderson 2001): pseudo-F, R2 and permutation p-value."""
    labels = np.asarray(labels)
    n = len(labels)
    d2 = dist ** 2
    ss_total = d2[np.triu_indices(n, 1)].sum() / n

    def pseudo_f(lab):
        groups = np.unique(lab)
        ss_within = 0.0
        for g in groups:
            idx = np.where(lab == g)[0]
            ss_within += d2[np.ix_(idx, idx)][np.triu_indices(len(idx), 1)].sum() / len(idx)
        ss_between = ss_total - ss_within
        k = len(groups)
        return (ss_between / (k - 1)) / (ss_within / (n - k)), ss_between / ss_total

    f_obs, r2 = pseudo_f(labels)
    hits = sum(pseudo_f(rng.permutation(labels))[0] >= f_obs for _ in range(n_perm))
    return f_obs, r2, (hits + 1) / (n_perm + 1)


def benjamini_hochberg(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    out = np.empty_like(adjusted)
    out[order] = adjusted
    return out


# --- Figures ------------------------------------------------------------------

def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG_DIR / name, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {name}")


def boxplot_by_group(values, groups, ylabel, title, filename, rng, log=False, hline=None):
    fig, ax = plt.subplots(figsize=(6, 5.5))
    data = [values[groups == g] for g in GROUP_ORDER]
    boxes = ax.boxplot(data, patch_artist=True, widths=0.5, medianprops={"color": "black"},
                       showfliers=log)
    for patch, g in zip(boxes["boxes"], GROUP_ORDER):
        patch.set_facecolor(GROUP_COLORS[g])
        patch.set_alpha(0.6)
    if not log:  # show every sample as a dot
        for i, (d, g) in enumerate(zip(data, GROUP_ORDER), start=1):
            jitter = rng.uniform(-0.12, 0.12, len(d))
            ax.scatter(np.full(len(d), i) + jitter, d, s=16, color=GROUP_COLORS[g],
                       edgecolor="white", zorder=3)
    if hline:
        ax.axhline(hline, color="black", linestyle="--", linewidth=1)
        ax.text(0.55, hline * 1.2, "rarefaction depth (1,000 reads)", fontsize=9)
    if log:
        ax.set_yscale("log")
    ax.set_xticks([1, 2])
    ax.set_xticklabels(GROUP_ORDER)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    save(fig, filename)


def phylum_bars(phyla, groups, top, colors):
    comp = phyla.groupby(groups).mean()[top].copy()
    comp["Other"] = (1 - comp.sum(axis=1)).clip(lower=0)
    comp = comp.loc[GROUP_ORDER] * 100
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    bottom = np.zeros(len(comp))
    for col in comp.columns:
        ax.bar(comp.index, comp[col], bottom=bottom, color=colors.get(col, "#D3D3D3"),
               label=col, width=0.55, edgecolor="white", linewidth=0.8)
        bottom += comp[col].to_numpy()
    ax.set_ylabel("Mean relative abundance (%)")
    ax.set_ylim(0, 100)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], frameon=False, bbox_to_anchor=(1.02, 1),
              loc="upper left", fontsize=12)
    ax.set_title("Phylum composition")
    save(fig, "01_phylum_composition.png")


def phylum_ring(phyla, groups, group, colors, filename):
    mean = phyla.loc[groups == group].drop(columns=["Unassigned"], errors="ignore").mean()
    shown = mean.sort_values(ascending=False).head(6).copy()
    shown["Other"] = 1 - shown.sum()
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    wedges, _, _ = ax.pie(
        shown.values, colors=[colors.get(p, "#BFBFBF") for p in shown.index],
        startangle=90, counterclock=False, pctdistance=0.78,
        autopct=lambda v: f"{v:.0f}%" if v >= 4 else "",
        wedgeprops={"width": 0.45, "edgecolor": "white", "linewidth": 1.5},
        textprops={"fontsize": 13, "color": "white", "weight": "bold"})
    ax.legend(wedges, [f"{p} ({v * 100:.1f}%)" for p, v in shown.items()],
              frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=14)
    ax.set_title(f"{group}: phylum composition", fontsize=15)
    ax.set_aspect("equal")
    save(fig, filename)


# --- Main ---------------------------------------------------------------------

def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    # Load, harmonize, filter
    raw_a = harmonize(load_table(DATA_DIR / "study_A_alaska_oklahoma.tsv", LABEL_A), LABEL_A)
    raw_b = harmonize(load_table(DATA_DIR / "study_B_atacama.tsv", LABEL_B), LABEL_B)
    qc_a, summary_a = quality_filter(raw_a, LABEL_A)
    qc_b, summary_b = quality_filter(raw_b, LABEL_B)
    depth = pd.concat([qc_a.sum().rename("reads").to_frame().assign(group=LABEL_A),
                       qc_b.sum().rename("reads").to_frame().assign(group=LABEL_B)])

    # Rarefy and merge into one sample-by-taxon table
    rar_a = rarefy(qc_a, RAREFY_DEPTH, rng)
    rar_b = rarefy(qc_b, RAREFY_DEPTH, rng)
    counts = pd.concat([rar_a, rar_b], axis=1).fillna(0).T
    groups = pd.Series([LABEL_A] * rar_a.shape[1] + [LABEL_B] * rar_b.shape[1], index=counts.index)
    rel = counts.div(counts.sum(axis=1), axis=0)
    in_a = (groups == LABEL_A).to_numpy()
    print(f"\n{counts.shape[0]} samples x {counts.shape[1]} taxa after rarefaction")

    # Alpha diversity
    alpha = pd.DataFrame({
        "observed_taxa": (counts > 0).sum(axis=1),
        "shannon": rel.apply(lambda row: shannon(row.to_numpy()), axis=1),
        "simpson": 1 - (rel ** 2).sum(axis=1),
        "group": groups,
    })
    alpha_tests = []
    for metric in ["observed_taxa", "shannon", "simpson"]:
        a, b = alpha.loc[in_a, metric], alpha.loc[~in_a, metric]
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        alpha_tests.append({"metric": metric, "median_A": a.median(), "median_B": b.median(),
                            "mean_A": a.mean(), "sd_A": a.std(), "mean_B": b.mean(),
                            "sd_B": b.std(), "U": u, "p_value": p})
    alpha_tests = pd.DataFrame(alpha_tests)
    print("\nAlpha diversity (Mann-Whitney)")
    print(alpha_tests.round(4).to_string(index=False))

    # Beta diversity
    bc = squareform(pdist(rel.to_numpy(), metric="braycurtis"))
    coords, var_exp = pcoa(bc)
    f_stat, r2, p_perm = permanova(bc, groups.to_numpy(), N_PERMUTATIONS, rng)
    bc_within_a = bc[np.ix_(in_a, in_a)][np.triu_indices(in_a.sum(), 1)].mean()
    bc_within_b = bc[np.ix_(~in_a, ~in_a)][np.triu_indices((~in_a).sum(), 1)].mean()
    bc_between = bc[np.ix_(in_a, ~in_a)].mean()
    taxa_a = set(rar_a.index[rar_a.sum(axis=1) > 0])
    taxa_b = set(rar_b.index[rar_b.sum(axis=1) > 0])
    print(f"\nPERMANOVA: pseudo-F = {f_stat:.2f}, R2 = {r2:.3f}, p = {p_perm:.3f}")
    print(f"Mean Bray-Curtis: within A {bc_within_a:.3f}, within B {bc_within_b:.3f}, "
          f"between {bc_between:.3f}")
    print(f"Taxa: {len(taxa_a)} in A, {len(taxa_b)} in B, {len(taxa_a & taxa_b)} shared")

    # Phylum composition
    phylum_of = {t: get_rank(t, "p") or "Unassigned" for t in counts.columns}
    phyla = rel.T.groupby(phylum_of).sum().T
    ranked = [p for p in phyla.mean().sort_values(ascending=False).index if p != "Unassigned"]
    top = ranked[:8]
    phylum_tests = []
    for ph in top:
        a, b = phyla.loc[in_a, ph], phyla.loc[~in_a, ph]
        phylum_tests.append({"phylum": ph, "mean_A_pct": a.mean() * 100,
                             "mean_B_pct": b.mean() * 100,
                             "p_value": mannwhitneyu(a, b, alternative="two-sided")[1]})
    phylum_tests = pd.DataFrame(phylum_tests)
    phylum_tests["p_adjusted_BH"] = benjamini_hochberg(phylum_tests["p_value"])
    phylum_tests = phylum_tests.sort_values("mean_A_pct", ascending=False)
    print("\nMain phyla")
    print(phylum_tests.round(5).to_string(index=False))

    # Figures (the same colour is used for a phylum in every figure)
    colors = {p: PHYLUM_PALETTE[i] for i, p in enumerate(ranked[:10])}
    colors["Other"] = "#D3D3D3"
    phylum_bars(phyla, groups, top, colors)
    boxplot_by_group(alpha["shannon"], groups, "Shannon index (rarefied to 1,000 reads)",
                     f"Shannon diversity (Mann-Whitney p = {alpha_tests.p_value[1]:.1e})",
                     "02_alpha_diversity_shannon.png", rng)
    boxplot_by_group(alpha["observed_taxa"], groups, "Observed taxa (per 1,000 reads)",
                     f"Observed richness (Mann-Whitney p = {alpha_tests.p_value[0]:.1e})",
                     "03_observed_richness.png", rng)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    for g in GROUP_ORDER:
        mask = (groups == g).to_numpy()
        ax.scatter(coords[mask, 0], coords[mask, 1], s=34, alpha=0.85,
                   color=GROUP_COLORS[g], edgecolor="white", label=g)
    ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)")
    ax.set_title(f"PCoA Bray-Curtis (PERMANOVA R² = {r2:.2f}, p = {p_perm:.3f})")
    ax.legend(frameon=False)
    save(fig, "04_pcoa_braycurtis.png")

    phylum_ring(phyla, groups, LABEL_A, colors, "05_pie_phylum_alaska_oklahoma.png")
    phylum_ring(phyla, groups, LABEL_B, colors, "06_pie_phylum_atacama.png")
    boxplot_by_group(depth["reads"], depth["group"], "Assigned reads per sample (log scale)",
                     "Sequencing depth after quality control", "07_sequencing_depth.png",
                     rng, log=True, hline=RAREFY_DEPTH)

    # Tables
    pd.DataFrame([summary_a, summary_b]).to_csv(TABLE_DIR / "quality_control.csv", index=False)
    alpha.to_csv(TABLE_DIR / "alpha_diversity_per_sample.csv")
    alpha_tests.to_csv(TABLE_DIR / "alpha_diversity_tests.csv", index=False)
    phylum_tests.to_csv(TABLE_DIR / "phylum_tests.csv", index=False)
    pd.DataFrame([{"pseudo_F": f_stat, "R2": r2, "p_value": p_perm,
                   "permutations": N_PERMUTATIONS, "bc_within_A": bc_within_a,
                   "bc_within_B": bc_within_b, "bc_between": bc_between,
                   "PC1_pct": var_exp[0], "PC2_pct": var_exp[1], "taxa_A": len(taxa_a),
                   "taxa_B": len(taxa_b), "taxa_shared": len(taxa_a & taxa_b)}]
                 ).to_csv(TABLE_DIR / "beta_diversity.csv", index=False)
    depth.to_csv(TABLE_DIR / "sequencing_depth.csv")
    print("\nDone: 7 figures in results/figures, 6 tables in results/tables")


if __name__ == "__main__":
    main()
