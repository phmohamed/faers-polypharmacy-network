"""Which drugs are reported together in FAERS adverse event reports? Drug co-report network, communities and
a look at whether reports with more drugs involve more deaths. Network code written from scratch (numpy)."""
import itertools, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(5)
r = pd.read_csv("data/reports.csv").fillna("")
r["list"] = r.drugs.map(lambda s: s.split("|") if s else [])
r["n"] = r.list.map(len)
print(len(r), "reports;", (r.n >= 2).sum(), "with 2+ named drugs")

# ---- polypharmacy vs death, within age bands ----
r["bucket"] = pd.cut(r.n, [0, 1, 2, 4, 9, 100], labels=["1", "2", "3-4", "5-9", "10+"])
r["ageband"] = pd.cut(r.age, [17, 44, 64, 79, 120], labels=["18-44", "45-64", "65-79", "80+"])
t = r.groupby("bucket", observed=True).agg(reports=("died", "size"), death_pct=("died", lambda s: round(s.mean() * 100, 1)))
ta = r.pivot_table(index="bucket", columns="ageband", values="died", aggfunc="mean", observed=True).mul(100).round(1)
t.to_csv("data/death_by_drug_count.csv"); ta.to_csv("data/death_by_drug_count_and_age.csv"); print(t); print(ta)

# ---- co-report network ----
multi = r[r.n >= 2]
cnt = pd.Series([d for l in multi.list for d in l]).value_counts()
keep = set(cnt[cnt >= 25].index); N = len(multi)
pair = {}
for l in multi.list:
    l = sorted(set(d for d in l if d in keep))
    for a, b in itertools.combinations(l, 2): pair[(a, b)] = pair.get((a, b), 0) + 1
E = pd.DataFrame([(a, b, c) for (a, b), c in pair.items() if c >= 10], columns=["a", "b", "together"])
E["lift"] = (E.together * N / (E.a.map(cnt) * E.b.map(cnt))).round(1)
E = E.sort_values("together", ascending=False); E.to_csv("data/drug_pairs.csv", index=False)
print("drugs", len(keep), "edges", len(E))
print(E.head(12).to_string()); print(E[E.together >= 20].sort_values("lift", ascending=False).head(12).to_string())

# graph: edges with lift >= 3 and at least 10 joint reports
G = E[E.lift >= 3]; nodes = sorted(set(G.a) | set(G.b)); ix = {n: i for i, n in enumerate(nodes)}
A = np.zeros((len(nodes), len(nodes)))
for a, b, c in zip(G.a, G.b, G.together): A[ix[a], ix[b]] = A[ix[b], ix[a]] = c
deg = (A > 0).sum(1); strength = A.sum(1)
# eigenvector centrality by power iteration
v = np.ones(len(nodes))
for _ in range(200): v = A @ v; v /= np.linalg.norm(v)
cent = pd.DataFrame({"drug": nodes, "partners": deg, "joint_reports": strength.astype(int), "eigenvector": v.round(3), "reports": [cnt[n] for n in nodes]}).sort_values("eigenvector", ascending=False)
cent.to_csv("data/centrality.csv", index=False); print(cent.head(12).to_string())

# communities: weighted label propagation, best of 20 runs by modularity
m2 = A.sum()
def modularity(lab):
    q = 0.0
    for c in set(lab):
        idx = np.where(lab == c)[0]; q += A[np.ix_(idx, idx)].sum() / m2 - (strength[idx].sum() / m2) ** 2
    return q
def lp():
    lab = np.arange(len(nodes))
    for _ in range(50):
        changed = False
        for i in rng.permutation(len(nodes)):
            if strength[i] == 0: continue
            w = {}
            for j in np.where(A[i] > 0)[0]: w[lab[j]] = w.get(lab[j], 0) + A[i, j]
            best = max(w.items(), key=lambda x: (x[1], rng.random()))[0]
            if best != lab[i]: lab[i] = best; changed = True
        if not changed: break
    return lab
labs = [lp() for _ in range(20)]; best = max(labs, key=modularity); Q = modularity(best)
print("modularity", round(Q, 3))
comm = pd.DataFrame({"drug": nodes, "community": best}); sizes = comm.community.value_counts()
big = sizes[sizes >= 4].index.tolist(); comm["community"] = comm.community.map({c: i for i, c in enumerate(big)}).fillna(-1).astype(int)
comm.to_csv("data/communities.csv", index=False)
for c in range(len(big)):
    mem = comm[comm.community == c].drug.tolist(); mem.sort(key=lambda d: -cnt[d]); print(c, len(mem), mem[:10])

# figure: network drawn with a simple spring layout (from scratch)
pos = rng.random((len(nodes), 2)) * 2 - 1
for it in range(300):
    D = pos[:, None, :] - pos[None, :, :]; dist = np.linalg.norm(D, axis=2) + 1e-3
    rep = (D / dist[:, :, None] ** 3).sum(1) * 0.002; att = -(((A > 0) * 1.0)[:, :, None] * D).sum(1) * 0.02
    pos += np.clip(rep + att, -0.05, 0.05)
fig, ax = plt.subplots(figsize=(9, 8)); cols = plt.cm.tab10(np.where(comm.community.values < 0, 9, comm.community.values % 9))
for a, b, c in zip(G.a, G.b, G.together): ax.plot(*zip(pos[ix[a]], pos[ix[b]]), color="grey", lw=0.3 + c / 40, alpha=0.4)
ax.scatter(pos[:, 0], pos[:, 1], s=20 + 3 * np.sqrt(np.array([cnt[n] for n in nodes])), c=cols, zorder=3)
for n in cent.head(18).drug: ax.annotate(n, pos[ix[n]], fontsize=7, zorder=4)
ax.axis("off"); ax.set_title("Drugs reported together in FAERS (2023 sample), colour = community"); fig.tight_layout(); fig.savefig("figures/network.png", dpi=150)
fig, ax = plt.subplots(figsize=(7, 4)); ta.plot(ax=ax, marker="o"); ax.set_ylabel("% of reports with a death"); ax.set_xlabel("Drugs named in the report"); ax.legend(title="Age", frameon=False)
fig.tight_layout(); fig.savefig("figures/death_by_drug_count.png", dpi=150)
