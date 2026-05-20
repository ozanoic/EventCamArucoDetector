"""
Generates the figures used by the algorithm-explanation PDF.

Each figure is saved as a PNG into docs/figures/.
Run from the repo root:
    python docs/make_figures.py
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D

FIG_DIR = Path(__file__).parent / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def save(fig, name: str) -> None:
    path = FIG_DIR / name
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {path}")


def boxed(ax, x, y, w, h, text, fc="#dfe7fd", ec="#3b4cca", fontsize=9):
    ax.add_patch(
        FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.05,rounding_size=0.10",
            linewidth=1.2, facecolor=fc, edgecolor=ec,
        )
    )
    ax.text(x + w / 2, y + h / 2, text,
            ha="center", va="center", fontsize=fontsize, wrap=True)


def arrow(ax, x1, y1, x2, y2, color="#444"):
    ax.annotate(
        "", xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=1.4,
                        shrinkA=2, shrinkB=4),
    )


# ---------------------------------------------------------------------------
# 1. End-to-end pipeline flowchart
# ---------------------------------------------------------------------------

def fig_pipeline():
    fig, ax = plt.subplots(figsize=(9, 10))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 14)
    ax.axis("off")

    nodes = [
        # (x, y, w, h, text, fill)
        (1.0, 12.5, 8.0, 1.0, "Input .mat file\n[x, y, polarity, timestamp]  (Nx4 events)", "#ffe1c2"),
        (1.0, 10.9, 8.0, 1.0, "Sort events by timestamp; choose multi-window set\n"
                              "windowDurations_ms = [150 200 ... 750]; tickStep = 1 ms", "#dfe7fd"),
        (1.0, 9.3, 8.0, 1.0, "Tick loop  (parallel parfor or sequential)\n"
                             "for tNow = tStart : 1 ms : tEnd ", "#dfe7fd"),
        (1.0, 7.7, 8.0, 1.0, "For each window dt: take events in [tNow-dt, tNow]\n"
                             "accumarray --> event count image (HxW)", "#dfe7fd"),
        (1.0, 6.1, 8.0, 1.0, "detectQuadBlob: imfill+watershed AND convex hull on edges\n"
                             "-> candidate 4-corner quads", "#dfe7fd"),
        (1.0, 4.5, 8.0, 1.0, "fitgeotrans('projective') + imwarp\n"
                             "warp 160x160 px square; canonical unwarped marker", "#dfe7fd"),
        (1.0, 2.9, 8.0, 1.0, "Scanline transition decode (vert+horiz+vote)\n"
                             "-> 6x6 inner code; rotate/flip/invert; ARUCO_MIP_36h12 lookup", "#dfe7fd"),
        (1.0, 1.3, 8.0, 1.0, "Filter by requestedMarkerIds; record win_Xms = marker ID or -1\n"
                             "(per-tick, per-window)", "#c8eed1"),
    ]
    for x, y, w, h, text, fc in nodes:
        boxed(ax, x, y, w, h, text, fc=fc, fontsize=10)

    for i in range(len(nodes) - 1):
        x = 5
        y1 = nodes[i][1]
        y2 = nodes[i + 1][1] + nodes[i + 1][3]
        arrow(ax, x, y1, x, y2)

    ax.text(5, 13.7, "EventCamArucoDetector  ·  end-to-end pipeline",
            ha="center", fontsize=13, fontweight="bold")
    save(fig, "01_pipeline.png")


# ---------------------------------------------------------------------------
# 2. Multi-window sliding tick scheme
# ---------------------------------------------------------------------------

def fig_sliding_tick():
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.4, 5.2)
    ax.axis("off")

    # event timeline ticks
    rng = np.random.default_rng(42)
    ev_t = np.sort(rng.uniform(0.2, 9.8, 350))
    ax.scatter(ev_t, np.zeros_like(ev_t) + 0.25, s=4, c="#888", alpha=0.55)
    ax.text(0, 0.25, "events  ", ha="right", va="center", fontsize=10)
    ax.plot([0, 10], [0.25, 0.25], color="#bbb", lw=0.6)

    # tNow marker
    tNow = 7.5
    ax.axvline(tNow, color="#cc1f1a", lw=1.6, ls="--")
    ax.text(tNow, 5.0, "tNow", ha="center", color="#cc1f1a", fontsize=11,
            fontweight="bold")

    # windows
    windows = [(0.8, 1.2, "150 ms"),
               (1.7, 2.0, "300 ms"),
               (2.6, 2.8, "500 ms"),
               (3.5, 3.6, "750 ms")]
    cmap = plt.cm.viridis(np.linspace(0.2, 0.9, len(windows)))
    for (left, y, lbl), color in zip(windows, cmap):
        ax.add_patch(Rectangle((tNow - left, y - 0.18), left, 0.36,
                               linewidth=1.2, edgecolor=color, facecolor=color,
                               alpha=0.30))
        ax.plot([tNow - left, tNow], [y, y], color=color, lw=2)
        ax.text(tNow - left - 0.05, y, lbl, ha="right", va="center",
                fontsize=9, color="black")

    ax.text(5, 4.7, "At every 1 ms tick, the detector tries N lookback windows in parallel.",
            ha="center", fontsize=10)
    ax.text(5, 4.3, "Each window collects the events inside [tNow - dt, tNow] only.",
            ha="center", fontsize=10, color="#444")

    ax.text(0.2, 0.7, "event timeline", fontsize=9, color="#666")
    save(fig, "02_sliding_tick.png")


# ---------------------------------------------------------------------------
# 3. Event accumulation -> event-count image
# ---------------------------------------------------------------------------

def fig_event_accumulation():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4),
                             gridspec_kw={"width_ratios": [1.1, 1]})

    # Left: raw events scatter inside the chosen window
    ax = axes[0]
    rng = np.random.default_rng(1)
    # simulate a marker quad with events on the edges
    H, W = 60, 80
    edges = []
    for t in np.linspace(0, 1, 220):
        # rotating square edges
        s = 20
        cx, cy = 40 + 4 * np.cos(2 * np.pi * t), 30 + 3 * np.sin(2 * np.pi * t)
        corners = np.array([[cx - s, cy - s], [cx + s, cy - s],
                            [cx + s, cy + s], [cx - s, cy + s]])
        for i in range(4):
            p1, p2 = corners[i], corners[(i + 1) % 4]
            n = 14
            xs = np.linspace(p1[0], p2[0], n) + rng.normal(0, 0.4, n)
            ys = np.linspace(p1[1], p2[1], n) + rng.normal(0, 0.4, n)
            edges.append(np.column_stack([xs, ys]))
    pts = np.vstack(edges)
    ax.scatter(pts[:, 0], pts[:, 1], s=4, c="#1a73e8", alpha=0.4)
    ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.set_aspect("equal")
    ax.set_title("(a) events in window\n(x, y) of every event in [tNow-dt, tNow]",
                 fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])

    # Right: countImg via accumarray
    img = np.zeros((H, W))
    xi = np.clip(pts[:, 0].astype(int), 0, W - 1)
    yi = np.clip(pts[:, 1].astype(int), 0, H - 1)
    np.add.at(img, (yi, xi), 1.0)
    ax = axes[1]
    ax.imshow(img, cmap="gray_r", interpolation="nearest")
    ax.set_title("(b) countImg = accumarray([y,x], 1, [H,W])\n"
                 "then activeMask = countImg > 0", fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])

    arrow_ax = fig.add_axes([0.49, 0.45, 0.02, 0.1])
    arrow_ax.axis("off")
    arrow_ax.annotate("", xy=(1, 0.5), xytext=(0, 0.5),
                      arrowprops=dict(arrowstyle="-|>", color="#222", lw=1.6))

    save(fig, "03_event_accumulation.png")


# ---------------------------------------------------------------------------
# 4. Quad detection - method A (fill+watershed) vs method B (convex hull on edges)
# ---------------------------------------------------------------------------

def fig_quad_detection():
    H, W = 60, 80
    rng = np.random.default_rng(3)

    # synthesize a square outline of events with one gap
    cx, cy, s = 40, 30, 18
    corners = np.array([[cx - s, cy - s], [cx + s, cy - s],
                        [cx + s, cy + s], [cx - s, cy + s]])
    edge_pts = []
    for i in range(4):
        p1, p2 = corners[i], corners[(i + 1) % 4]
        n = 90
        if i == 1:
            n = 30  # gap on the right edge
        xs = np.linspace(p1[0], p2[0], n) + rng.normal(0, 0.3, n)
        ys = np.linspace(p1[1], p2[1], n) + rng.normal(0, 0.3, n)
        edge_pts.append(np.column_stack([xs, ys]))
    pts = np.vstack(edge_pts)

    img = np.zeros((H, W), bool)
    xi = np.clip(pts[:, 0].astype(int), 0, W - 1)
    yi = np.clip(pts[:, 1].astype(int), 0, H - 1)
    img[yi, xi] = True

    # method A: fill the holes
    from scipy.ndimage import binary_fill_holes  # ships with matplotlib? not always
    try:
        filled = binary_fill_holes(img)
    except Exception:
        filled = img.copy()

    fig, axes = plt.subplots(1, 3, figsize=(11, 4))

    axes[0].imshow(img, cmap="gray_r")
    axes[0].set_title("activeMask\n(thin event edges, with a gap)", fontsize=10)

    axes[1].imshow(filled, cmap="gray_r")
    axes[1].set_title("Method A: imfill + watershed\n-> filled blob, fit min-area rect",
                      fontsize=10)
    # overlay rect
    axes[1].add_patch(Rectangle((cx - s, cy - s), 2 * s, 2 * s,
                                 fill=False, edgecolor="#cc1f1a", lw=2))

    axes[2].imshow(img, cmap="gray_r")
    axes[2].set_title("Method B: convex hull of raw edge component\n"
                       "-> robust to gaps", fontsize=10)
    hull_x = [cx - s, cx + s, cx + s, cx - s, cx - s]
    hull_y = [cy - s, cy - s, cy + s, cy + s, cy - s]
    axes[2].plot(hull_x, hull_y, color="#1a73e8", lw=2)
    axes[2].scatter([cx - s, cx + s, cx + s, cx - s],
                    [cy - s, cy - s, cy + s, cy + s],
                    s=42, color="#1a73e8", zorder=5)

    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect("equal")

    save(fig, "04_quad_detection.png")


# ---------------------------------------------------------------------------
# 5. Perspective unwarp
# ---------------------------------------------------------------------------

def fig_unwarp():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))

    # Left: distorted quad on event image
    ax = axes[0]
    ax.set_xlim(0, 80); ax.set_ylim(60, 0)
    ax.set_aspect("equal")
    quad = np.array([[18, 12], [62, 22], [70, 52], [10, 46]])
    poly = plt.Polygon(quad, fill=False, edgecolor="#cc1f1a", lw=2)
    ax.add_patch(poly)
    ax.scatter(quad[:, 0], quad[:, 1], s=46, color="#cc1f1a", zorder=5)
    for (x, y), label in zip(quad, ["TL", "TR", "BR", "BL"]):
        ax.annotate(label, (x, y), textcoords="offset points",
                    xytext=(8, -2), fontsize=9, color="#cc1f1a")
    # draw fake event speckle inside
    rng = np.random.default_rng(7)
    for _ in range(180):
        u = rng.random(); v = rng.random()
        p = (1 - u) * (1 - v) * quad[0] + u * (1 - v) * quad[1] + \
            u * v * quad[2] + (1 - u) * v * quad[3]
        if rng.random() < 0.5:
            ax.plot(p[0], p[1], '.', color="#444", markersize=2)
    ax.set_title("(a) Detected quad on event image\nsrcCorners (TL, TR, BR, BL)",
                 fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])

    # Right: warped square
    ax = axes[1]
    side = 160
    ax.set_xlim(-10, side + 10); ax.set_ylim(side + 10, -10)
    ax.set_aspect("equal")
    ax.add_patch(Rectangle((0, 0), side, side, fill=False,
                           edgecolor="#1a73e8", lw=2))
    # 8x8 grid lines
    for k in range(9):
        ax.plot([0, side], [k * 20, k * 20], color="#aac", lw=0.6)
        ax.plot([k * 20, k * 20], [0, side], color="#aac", lw=0.6)
    # corner labels
    for (x, y), label in zip(
        [(0, 0), (side - 1, 0), (side - 1, side - 1), (0, side - 1)],
        ["TL", "TR", "BR", "BL"]):
        ax.annotate(label, (x, y), textcoords="offset points",
                    xytext=(6, 14), fontsize=9, color="#1a73e8")
    ax.set_title("(b) After fitgeotrans(src, dst, 'projective') + imwarp\n"
                 "160x160 px canonical marker (8x8 cells of 20 px)",
                 fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])

    # arrow between
    fig.add_artist(FancyArrowPatch((0.487, 0.50), (0.518, 0.50),
                                    transform=fig.transFigure,
                                    arrowstyle="-|>", color="#222", lw=1.8,
                                    mutation_scale=18))
    save(fig, "05_unwarp.png")


# ---------------------------------------------------------------------------
# 6. Scanline transition decoding
# ---------------------------------------------------------------------------

def fig_scanline_decode():
    # Build a synthetic 8x8 ARUCO marker with a black border, 6x6 random inner
    rng = np.random.default_rng(11)
    inner = rng.integers(0, 2, size=(6, 6))
    grid = np.ones((8, 8))               # 1 = white
    grid[0, :] = 0; grid[-1, :] = 0      # black border
    grid[:, 0] = 0; grid[:, -1] = 0
    grid[1:-1, 1:-1] = inner             # inner code (1=white, 0=black)

    img = np.kron(grid, np.ones((20, 20)))  # 160x160 px

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.6))
    ax = axes[0]
    ax.imshow(img, cmap="gray", vmin=0, vmax=1)
    # mark cell boundaries
    for k in range(9):
        ax.axvline(k * 20 - 0.5, color="#cc1f1a", lw=0.5, alpha=0.6)
        ax.axhline(k * 20 - 0.5, color="#cc1f1a", lw=0.5, alpha=0.6)
    ax.set_title("(a) Unwarped marker (8x8 cells)\n"
                 "red lines = cell boundaries we sample", fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])

    # vertical scanline highlight
    ax = axes[1]
    ax.imshow(img, cmap="gray", vmin=0, vmax=1)
    for col in range(8):
        ax.axvline(col * 20 + 10, color="#1a73e8", lw=1.4, alpha=0.8)
    for row in range(1, 8):
        ax.axhline(row * 20 - 0.5, color="#cc1f1a", lw=1.0, alpha=0.7)
    ax.set_title("(b) Vertical scanlines (blue)\n"
                 "flip current colour at every boundary that crosses threshold",
                 fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])

    # Combined codes
    ax = axes[2]
    ax.imshow(grid, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    ax.set_title("(c) Decoded 8x8 grid\ncode_V, code_H, majority(V,H) all tried",
                 fontsize=10)
    for r in range(8):
        for c in range(8):
            ax.text(c, r, "1" if grid[r, c] else "0", color="#cc1f1a"
                    if grid[r, c] else "#ffeb3b", ha="center", va="center",
                    fontsize=8)
    ax.set_xticks([]); ax.set_yticks([])

    save(fig, "06_scanline_decode.png")


# ---------------------------------------------------------------------------
# 7. ARUCO_MIP_36h12 lookup with rotation/flip/invert
# ---------------------------------------------------------------------------

def fig_dictionary_lookup():
    fig, ax = plt.subplots(figsize=(9, 5.4))
    ax.set_xlim(0, 10); ax.set_ylim(0, 6)
    ax.axis("off")

    boxed(ax, 0.5, 4.5, 2.4, 1.0, "Inner 6x6 bits\n(testCode)", fc="#ffe1c2")
    boxed(ax, 3.6, 5.2, 2.4, 0.6, "inv = 0 or 1\n(invert colours)", fc="#dfe7fd")
    boxed(ax, 6.4, 5.2, 2.4, 0.6, "flip = 0 or 1\n(fliplr)", fc="#dfe7fd")
    boxed(ax, 6.4, 4.3, 2.4, 0.6, "rot in {0, 90, 180, 270}", fc="#dfe7fd")
    boxed(ax, 3.6, 2.7, 5.0, 1.0,
          "Pack 36 bits into a uint64 code\n"
          "search sorted dictionary (250 codes) via binary search",
          fc="#fff3cd")
    boxed(ax, 0.5, 2.7, 2.4, 1.0,
          "candidates =\n{ codeV, codeH,\nmajority(V,H) }",
          fc="#ffe1c2")

    boxed(ax, 3.6, 0.9, 5.0, 1.2,
          "match found  -->  return marker ID (0..249)\n"
          "no match in all 2*2*4*3 = 48 variants  -->  ID = -1",
          fc="#c8eed1")

    arrow(ax, 2.9, 5.0, 3.6, 5.5)
    arrow(ax, 6.0, 5.5, 6.4, 5.5)
    arrow(ax, 6.0, 4.6, 6.4, 4.6)
    arrow(ax, 2.9, 3.2, 3.6, 3.2)
    arrow(ax, 7.6, 4.3, 7.6, 3.7)
    arrow(ax, 6.1, 2.7, 6.1, 2.1)

    ax.set_title("Decoder tries 48 variants per candidate before giving up",
                 fontsize=11, fontweight="bold")
    save(fig, "07_dictionary_lookup.png")


# ---------------------------------------------------------------------------
# 8. Merge: union timeline across two runs
# ---------------------------------------------------------------------------

def fig_merge_union():
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.set_xlim(0, 10); ax.set_ylim(-0.3, 4.2)
    ax.axis("off")

    # v1 timeline
    ax.text(-0.05, 3.0, "v1 (windows 5..150 ms)\n4826 ticks",
            ha="right", va="center", fontsize=9)
    ax.add_patch(Rectangle((0.15, 2.85), 8.5, 0.30, facecolor="#dfe7fd",
                           edgecolor="#3b4cca"))

    # v2 timeline (starts later because max window grew to 750 ms)
    ax.text(-0.05, 2.0, "v2 (windows 150..750 ms)\n4226 ticks",
            ha="right", va="center", fontsize=9)
    ax.add_patch(Rectangle((0.75, 1.85), 7.9, 0.30, facecolor="#dfe7fd",
                           edgecolor="#3b4cca"))

    # union timeline
    ax.text(-0.05, 1.0, "union (mergeResults)\n4826 ticks",
            ha="right", va="center", fontsize=9)
    ax.add_patch(Rectangle((0.15, 0.85), 8.5, 0.30, facecolor="#c8eed1",
                           edgecolor="#286f43"))

    # callout: -1 fill in v2 columns for the first 600 ticks
    ax.annotate("", xy=(0.45, 0.85), xytext=(0.45, 1.55),
                arrowprops=dict(arrowstyle="-|>", color="#cc1f1a", lw=1.4))
    ax.text(0.55, 1.6,
            "for these 600 ticks, v2's columns\nare marked -1 (not-attempted)",
            fontsize=8, color="#cc1f1a")

    ax.text(5, 3.8, "Merging files with different window sets keeps every tick.",
            ha="center", fontsize=11, fontweight="bold")
    ax.text(5, 0.4,
            "Per-window rate uses attemptedPerWindow as denominator;\n"
            "AnyDetectPct uses the full union length.",
            ha="center", fontsize=9, color="#444")

    save(fig, "08_merge_union.png")


# ---------------------------------------------------------------------------
# 9. Viewer mock - tab layout schematic
# ---------------------------------------------------------------------------

def fig_viewer_mock():
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.5); ax.axis("off")

    # window border
    ax.add_patch(Rectangle((0.1, 0.1), 9.8, 6.3, fill=False, edgecolor="#222",
                            lw=1.4))

    # toolbar
    ax.add_patch(Rectangle((0.1, 5.7), 9.8, 0.7, facecolor="#eee",
                            edgecolor="#222"))
    boxed(ax, 0.3, 5.85, 0.6, 0.4, "Marker:", fc="#fff", ec="#aaa", fontsize=8)
    boxed(ax, 1.0, 5.85, 1.4, 0.4, "All markers v", fc="#fff", ec="#aaa", fontsize=8)
    boxed(ax, 2.7, 5.85, 1.7, 0.4, "Save Current Tab", fc="#fff", ec="#aaa", fontsize=8)
    boxed(ax, 4.6, 5.85, 1.8, 0.4, "Save Whole Window", fc="#fff", ec="#aaa", fontsize=8)
    boxed(ax, 6.6, 5.85, 2.2, 0.4, "Save All Tabs (PDF)", fc="#fff", ec="#aaa", fontsize=8)

    # tabs row
    tabs = ["Summary", "cross low", "cross med", "cross high",
            "linear low", "linear med", "...", "zoom high"]
    x = 0.3
    for t in tabs:
        w = 0.95
        fc = "#fff" if t != "Summary" else "#dfe7fd"
        ax.add_patch(Rectangle((x, 5.30), w, 0.35, facecolor=fc,
                                edgecolor="#aaa"))
        ax.text(x + w / 2, 5.475, t, ha="center", va="center", fontsize=8)
        x += w + 0.04

    # content area  (Summary tab visible)
    boxed(ax, 0.3, 2.6, 4.7, 2.5,
          "Per-Dataset Summary  (table)\n"
          "Dataset | Ticks | AnyDetectPct | BestWindow | BestRatePct",
          fc="#fff", ec="#aaa", fontsize=9)
    boxed(ax, 5.1, 2.6, 4.6, 2.5,
          "Overall Detection Rate per Dataset\n(bar chart)",
          fc="#fff", ec="#aaa", fontsize=9)
    boxed(ax, 0.3, 0.3, 9.4, 2.1,
          "Detection Rate (%) - Dataset x Window  (heatmap)\n"
          "rows = datasets sorted low->med->high, cols = window durations",
          fc="#fff", ec="#aaa", fontsize=9)

    ax.text(5, 6.15, "viewAllResults('Data')  -  tabbed GUI",
            ha="center", fontsize=11, fontweight="bold")

    save(fig, "09_viewer_mock.png")


# ---------------------------------------------------------------------------
# 10. Detail-tab mock
# ---------------------------------------------------------------------------

def fig_detail_mock():
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.set_xlim(0, 10); ax.set_ylim(0, 5.8); ax.axis("off")
    ax.add_patch(Rectangle((0.1, 0.1), 9.8, 5.6, fill=False, edgecolor="#222",
                            lw=1.4))

    ax.text(5, 5.4,
            "zoom_med  |  marker ID 3  |  4826 ticks (4.83 s)  |  "
            "Any-window: 4599 (95.3%)  |  Best: 500 ms at 100.0%",
            ha="center", fontsize=9, fontweight="bold")

    # raster
    ax.add_patch(Rectangle((0.3, 2.2), 9.4, 2.8, facecolor="#fafafa",
                            edgecolor="#aaa"))
    ax.text(5, 4.85, "Detection Timeline (per window)", ha="center", fontsize=9)
    rng = np.random.default_rng(13)
    for wi, y in enumerate(np.linspace(2.4, 4.6, 13)):
        n = int(450 - wi * 6 + rng.normal(0, 10))
        xs = rng.uniform(0.5, 9.5, n)
        ax.plot(xs, np.full_like(xs, y), '|', markersize=5,
                color=plt.cm.hsv(wi / 13))
    ax.text(0.5, 4.55, "150ms", fontsize=7)
    ax.text(0.5, 2.45, "750ms", fontsize=7)

    # bar
    ax.add_patch(Rectangle((0.3, 0.3), 9.4, 1.7, facecolor="#fafafa",
                            edgecolor="#aaa"))
    ax.text(5, 1.85, "Detection Rate per Window  (rate = detections / attempts)",
            ha="center", fontsize=9)
    wins = np.arange(13)
    rates = 70 + 25 * (1 - np.exp(-0.5 * wins)) + rng.normal(0, 1.2, 13)
    rates = np.clip(rates, 0, 100)
    ax.bar(0.5 + wins * 0.7, rates / 60, width=0.55,
           color=plt.cm.hsv(wins / 13), bottom=0.5)
    for w, r in zip(wins, rates):
        ax.text(0.5 + w * 0.7, 0.5 + r / 60 + 0.04, f"{r:.0f}%",
                ha="center", fontsize=6)
    save(fig, "10_detail_mock.png")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print(f"Writing figures into {FIG_DIR}")
    fig_pipeline()
    fig_sliding_tick()
    fig_event_accumulation()
    fig_quad_detection()
    fig_unwarp()
    fig_scanline_decode()
    fig_dictionary_lookup()
    fig_merge_union()
    fig_viewer_mock()
    fig_detail_mock()
    print("done.")
