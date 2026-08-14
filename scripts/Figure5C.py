"""
SPHK1_network_revision2.py
Revision 2: same layout/colours as Revision_GeneticNetwork, but use a DIFFERENT
marker shape for GWAS disease-trait genes (squares) vs sphingolipid-measurement
genes / hubs (circles).
Output: Revision2/Figure5C_geneticNetwork_v2.{png,svg}
"""

import os, zipfile, csv, math
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.family"]  = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Nimbus Sans", "Helvetica", "DejaVu Sans"]
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

# ── Paths (relative to this script) ──────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Inputs live in the data/ folder parallel to scripts/
WORK = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "data"))

# INPUT: Cytoscape session file. The whole pipeline starts from this.
CYS  = os.path.join(WORK, "SPHK1_ANDIE.cys")

OUT_DIR = SCRIPT_DIR   # write outputs next to this script (./)
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PNG = os.path.join(OUT_DIR, "Figure5C_geneticNetwork_v2.png")
OUT_SVG = os.path.join(OUT_DIR, "Figure5C_geneticNetwork_v2.svg")

# ── 1. Start from the .cys file ──────────────────────────────────────────────
# Validate the source file exists, then extract its contents into a script-local
# folder. Set FORCE_REEXTRACT=True (or delete .cys_extracted) for a clean re-unpack.
FORCE_REEXTRACT = False
EXTRACT = os.path.join(SCRIPT_DIR, ".cys_extracted")

if not os.path.isfile(CYS):
    raise FileNotFoundError(
        f"Cytoscape session file not found: {CYS}\n"
        f"This script must start from SPHK1_ANDIE.cys located in: {WORK}"
    )

print(f"Source .cys     : {CYS}")
print(f"Extracting into : {EXTRACT}")

if FORCE_REEXTRACT and os.path.isdir(EXTRACT):
    import shutil
    shutil.rmtree(EXTRACT)

needs_extract = (not os.path.isdir(EXTRACT)) or (
    not any(d.startswith("CytoscapeSession") for d in os.listdir(EXTRACT))
)
if needs_extract:
    os.makedirs(EXTRACT, exist_ok=True)
    with zipfile.ZipFile(CYS, "r") as zf:
        zf.extractall(EXTRACT)
    print("Extracted .cys archive.")
else:
    print("Reusing previously extracted .cys archive.")

session_dir = next(
    os.path.join(EXTRACT, d) for d in os.listdir(EXTRACT)
    if d.startswith("CytoscapeSession")
)

# ── 2. Parse view XGMML → exact node positions ────────────────────────────────
NS = "http://www.cs.rpi.edu/XGMML"
view_file = next(
    os.path.join(session_dir, "views", f)
    for f in os.listdir(os.path.join(session_dir, "views"))
    if f.endswith(".xgmml")
)
vtree = ET.parse(view_file)
pos_raw = {}
for n in vtree.getroot().iter(f"{{{NS}}}node"):
    label = n.get("label", "")
    g     = n.find(f"{{{NS}}}graphics")
    if g is not None and label and not label.endswith("network"):
        x, y = float(g.get("x", 0)), float(g.get("y", 0))
        pos_raw[label] = (x, -y)

# ── 3. Parse network XGMML → edges ───────────────────────────────────────────
net_file = next(
    os.path.join(session_dir, "networks", f)
    for f in os.listdir(os.path.join(session_dir, "networks"))
    if f.endswith(".xgmml")
)
ntree = ET.parse(net_file)
net_nodes = {n.get("id"): n.get("label","")
             for n in ntree.getroot().iter(f"{{{NS}}}node")
             if not n.get("label","").endswith("network")}
raw_edges = [
    (e.get("source",""), e.get("target",""))
    for e in ntree.getroot().iter(f"{{{NS}}}edge")
]

# ── 4. Parse shared node cytable → degree + Disease_Trait + group ─────────────
shared_table = None
for dp, _, files in os.walk(os.path.join(session_dir, "tables")):
    for f in files:
        if "SHARED_ATTRS" in f and "CyNode" in f and "root+shared" in f:
            shared_table = os.path.join(dp, f); break
    if shared_table: break

node_attrs = {}
with open(shared_table, encoding="utf-8") as fh:
    lines = fh.readlines()
for row in csv.reader(lines[5:]):
    if len(row) < 7: continue
    suid, degree, disease, group = row[0], row[1], row[2], row[3]
    name = row[6].strip('"')
    disease = disease.strip('"'); group = group.strip('"')
    try:    deg = int(degree.strip('"'))
    except: deg = 1
    if name:
        node_attrs[name] = {"degree": deg, "disease": disease, "group": group}

# ── 5. Build plottable nodes and edges ────────────────────────────────────────
_raw_nodes = {
    lbl: xy for lbl, xy in pos_raw.items()
    if not lbl.lstrip("-").isdigit() and not lbl.endswith("network")
}
_cx = sum(v[0] for v in _raw_nodes.values()) / len(_raw_nodes)
nodes_plot = {lbl: ((_cx + (x - _cx) * 1.20), y)
              for lbl, (x, y) in _raw_nodes.items()}
edges_plot = [
    (net_nodes[s], net_nodes[t])
    for s, t in raw_edges
    if net_nodes.get(s,"") in nodes_plot and net_nodes.get(t,"") in nodes_plot
]
print(f"Nodes: {len(nodes_plot)}  |  Edges: {len(edges_plot)}")

# ── 6. Visual style ──────────────────────────────────────────────────────────
HUB_NODES = {"SPHK1", "S1PR1", "SPHK2"}

HUB_COLOUR = {
    "SPHK1": "#C94060",
    "SPHK2": "#7B40AA",
    "S1PR1": "#E8902A",
}

def node_category(label):
    if label in HUB_NODES:
        return f"hub_{label}"
    attr    = node_attrs.get(label, {})
    disease = attr.get("disease", "")
    grp     = attr.get("group",   "")
    if disease == "sphingolipid measurement":
        return "sphingolipid"
    if grp == "SPHK2" or disease == "SPHK2":
        return "sphk2_specific"
    if grp == "S1PR1" or disease == "S1PR1":
        return "s1pr1_specific"
    return "ibd"

CAT_COLOUR = {
    "hub_SPHK1":       HUB_COLOUR["SPHK1"],
    "hub_S1PR1":       HUB_COLOUR["S1PR1"],
    "hub_SPHK2":       HUB_COLOUR["SPHK2"],
    "sphingolipid":    "#F0C040",
    "ibd":             "#22C8E0",
    "s1pr1_specific":  "#A07BD0",   # purple (lighter than SPHK2 hub purple)
    "sphk2_specific":  "#A2A2A2",   # E8E8E8 darkened ~30%
}
CAT_EDGECOLOUR = {
    "hub_SPHK1":       "#FFFFFF",
    "hub_S1PR1":       "#FFFFFF",
    "hub_SPHK2":       "#FFFFFF",
    "sphingolipid":    "#C8981A",
    "ibd":             "#18A0B8",
    "s1pr1_specific":  "#6E4FA0",
    "sphk2_specific":  "#6E6E6E",
}
CAT_LINEWIDTH = {k: (1.8 if k.startswith("hub") else 0.9) for k in CAT_COLOUR}

# Shape per category — KEY CHANGE FOR REVISION 2
# GWAS disease-trait genes (IBD, S1PR1-specific, SPHK2-specific) → squares
# Sphingolipid-measurement genes → circles
# Hubs → circles
GWAS_CATS = {"ibd", "s1pr1_specific", "sphk2_specific"}
# Shape per category (Revision 2 — per-category shapes):
#   sphingolipid measurement → circle
#   SPHK1 neighbors in IBD   → square
#   S1PR1-specific in IBD    → star
#   SPHK2-specific in IBD    → triangle (up)
#   Hubs                     → circle (kept for visual prominence)
CAT_MARKER = {
    "hub_SPHK1":       "o",
    "hub_S1PR1":       "o",
    "hub_SPHK2":       "o",
    "sphingolipid":    "o",
    "ibd":             "s",
    "s1pr1_specific":  "*",
    "sphk2_specific":  "^",
}

def node_size(label):
    attr = node_attrs.get(label, {})
    deg  = attr.get("degree", 1)
    if label == "SPHK1": return deg * 65 * 1.5
    if label == "S1PR1": return deg * 68 * 1.5
    if label == "SPHK2": return deg * 105 * 1.5
    cat = node_category(label)
    # Match perceptual size across markers: squares are largest at same `s`,
    # triangles smaller, stars smallest. Compensate so all satellites read
    # roughly the same on screen.
    if cat == "ibd":            return 200   # square
    if cat == "s1pr1_specific": return 360   # star
    if cat == "sphk2_specific": return 280   # triangle
    return 230                                # circle (sphingolipid)

# ── 7. Render ─────────────────────────────────────────────────────────────────
print("Rendering …")
fig, ax = plt.subplots(figsize=(21.6, 10))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")
ax.set_aspect("equal")
ax.axis("off")

segs = [
    [nodes_plot[s], nodes_plot[t]]
    for s, t in edges_plot
    if s in nodes_plot and t in nodes_plot
]
ax.add_collection(LineCollection(segs, colors="#999999", linewidths=1.4,
                                 alpha=0.85, zorder=1))

# Group nodes by marker shape so we can call scatter once per shape
node_list = list(nodes_plot.keys())
by_marker = {}
for n in node_list:
    cat = node_category(n)
    by_marker.setdefault(CAT_MARKER[cat], []).append(n)

for marker, names in by_marker.items():
    xs   = [nodes_plot[n][0] for n in names]
    ys   = [nodes_plot[n][1] for n in names]
    cats = [node_category(n) for n in names]
    facecolours = [CAT_COLOUR[c]     for c in cats]
    edgecolours = [CAT_EDGECOLOUR[c] for c in cats]
    linewidths  = [CAT_LINEWIDTH[c]  for c in cats]
    sizes       = [node_size(n)      for n in names]
    ax.scatter(xs, ys,
               s          = sizes,
               c          = facecolours,
               edgecolors = edgecolours,
               linewidths = linewidths,
               marker     = marker,
               zorder     = 2,
               alpha      = 0.92)

# ── 7c. Labels ────────────────────────────────────────────────────────────────
hub_pos = {h: nodes_plot[h] for h in HUB_NODES if h in nodes_plot}

for label in node_list:
    x, y = nodes_plot[label]

    if label in HUB_NODES:
        ax.text(x, y, label,
                fontsize=20, fontweight="bold", color="white",
                ha="center", va="center", zorder=4)
    else:
        nearest_hub = min(hub_pos.keys(),
                          key=lambda h: math.hypot(x - hub_pos[h][0],
                                                   y - hub_pos[h][1]))
        hx, hy = hub_pos[nearest_hub]
        dx, dy = x - hx, y - hy
        dist   = math.hypot(dx, dy) or 1.0
        offset_data = 10
        lx = x + dx / dist * offset_data
        ly = y + dy / dist * offset_data

        ax.text(lx, ly, label,
                fontsize   = 17,
                fontweight = "normal",
                color      = "#222222",
                ha         = "center",
                va         = "center",
                zorder     = 4)

# ── 7d. Title ─────────────────────────────────────────────────────────────────
all_x = [v[0] for v in nodes_plot.values()]
all_y = [v[1] for v in nodes_plot.values()]
cx = (min(all_x) + max(all_x)) / 2

ax.text(cx, max(all_y) + 30,
        "Genetic network enrichment analysis",
        fontsize=13, color="#888888",
        ha="center", va="bottom", zorder=5,
        fontfamily="Arial")

# ── 7e. Legend ────────────────────────────────────────────────────────────────
# Two visual cues now: colour AND shape. Legend reflects both.
legend_elements = [
    Line2D([0], [0], marker="o", color="w",
           label="Sphingolipid measurement",
           markerfacecolor="#F0C040", markersize=10,
           markeredgecolor="#C8981A", markeredgewidth=0.9, linewidth=0),
    Line2D([0], [0], marker="s", color="w",
           label="SPHK1 neighbors in IBD",
           markerfacecolor="#22C8E0", markersize=10,
           markeredgecolor="#18A0B8", markeredgewidth=0.9, linewidth=0),
    Line2D([0], [0], marker="*", color="w",
           label="S1PR1-specific neighbors in IBD",
           markerfacecolor="#A07BD0", markersize=14,
           markeredgecolor="#6E4FA0", markeredgewidth=0.9, linewidth=0),
    Line2D([0], [0], marker="^", color="w",
           label="SPHK2-specific neighbors in IBD",
           markerfacecolor="#A2A2A2", markersize=11,
           markeredgecolor="#6E6E6E", markeredgewidth=0.9, linewidth=0),
]

ax.legend(
    handles        = legend_elements,
    loc            = "upper right",
    fontsize       = 9.5,
    frameon        = False,
    handletextpad  = 0.7,
    labelspacing   = 0.65,
    borderpad      = 0.5,
    title          = "GWAS genes",
    title_fontsize = 11,
)

# ── 7f. Axis limits with padding ─────────────────────────────────────────────
pad_x = (max(all_x) - min(all_x)) * 0.08
pad_y = (max(all_y) - min(all_y)) * 0.12
ax.set_xlim(min(all_x) - pad_x, max(all_x) + pad_x)
ax.set_ylim(min(all_y) - pad_y, max(all_y) + pad_y + 50)

plt.tight_layout(pad=0.3)

plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
print(f"Saved PNG : {OUT_PNG}")

plt.savefig(OUT_SVG, format="svg", bbox_inches="tight", facecolor="white",
            metadata={"Creator": "SPHK1_network_revision2.py"})
print(f"Saved SVG : {OUT_SVG}")

plt.close()
