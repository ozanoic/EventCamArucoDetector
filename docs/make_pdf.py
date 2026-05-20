"""
Builds docs/algorithm_explanation.pdf -- a detailed walk-through of the
EventCamArucoDetector pipeline with diagrams.

Run from the repo root after make_figures.py:
    python docs/make_pdf.py
"""

from __future__ import annotations

from pathlib import Path

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT, TA_JUSTIFY
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak,
    KeepTogether, Table, TableStyle,
)
from reportlab.pdfgen import canvas as _canvas

DOC_DIR = Path(__file__).parent
FIG_DIR = DOC_DIR / "figures"
OUT_PDF = DOC_DIR / "algorithm_explanation.pdf"


# ---------------------------------------------------------------------------
# Styles
# ---------------------------------------------------------------------------

styles = getSampleStyleSheet()

H1 = ParagraphStyle("H1",
                    parent=styles["Heading1"],
                    fontSize=18, spaceBefore=8, spaceAfter=14,
                    textColor=colors.HexColor("#1a237e"))

H2 = ParagraphStyle("H2",
                    parent=styles["Heading2"],
                    fontSize=14, spaceBefore=14, spaceAfter=8,
                    textColor=colors.HexColor("#283593"))

H3 = ParagraphStyle("H3",
                    parent=styles["Heading3"],
                    fontSize=11.5, spaceBefore=10, spaceAfter=6,
                    textColor=colors.HexColor("#37474f"))

BODY = ParagraphStyle("Body",
                      parent=styles["BodyText"],
                      fontSize=10.5, leading=15,
                      spaceBefore=2, spaceAfter=6,
                      alignment=TA_JUSTIFY)

BULLET = ParagraphStyle("Bullet",
                        parent=BODY,
                        leftIndent=18, bulletIndent=6,
                        spaceBefore=0, spaceAfter=2)

CODE = ParagraphStyle("Code",
                      parent=styles["Code"],
                      fontName="Courier", fontSize=8.5,
                      leftIndent=10, leading=11,
                      backColor=colors.HexColor("#f4f5f7"),
                      borderColor=colors.HexColor("#cfd8dc"),
                      borderWidth=0.5, borderPadding=4,
                      spaceBefore=4, spaceAfter=6)

CAPTION = ParagraphStyle("Caption",
                         parent=BODY,
                         fontSize=9.0, leading=11,
                         textColor=colors.HexColor("#444"),
                         alignment=TA_LEFT,
                         spaceBefore=2, spaceAfter=12)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def fig(name: str, width_cm: float, caption: str | None = None):
    """Image flowable scaled to width_cm; optional caption beneath."""
    path = FIG_DIR / name
    img = Image(str(path), width=width_cm * cm, height=None)
    # preserve aspect ratio: compute height from native size
    iw, ih = img.imageWidth, img.imageHeight
    img.drawWidth = width_cm * cm
    img.drawHeight = (ih / iw) * width_cm * cm
    if caption:
        return KeepTogether([img, Paragraph(caption, CAPTION)])
    return img


def code(text: str):
    safe = (text.replace("&", "&amp;")
                .replace("<", "&lt;")
                .replace(">", "&gt;")
                .replace("\n", "<br/>")
                .replace(" ", "&nbsp;"))
    return Paragraph(safe, CODE)


def p(text: str):
    return Paragraph(text, BODY)


def b(text: str):
    return Paragraph(text, BULLET, bulletText="•")


def h1(text: str):
    return Paragraph(text, H1)


def h2(text: str):
    return Paragraph(text, H2)


def h3(text: str):
    return Paragraph(text, H3)


# ---------------------------------------------------------------------------
# Page footer with page numbers
# ---------------------------------------------------------------------------

def _draw_footer(c: _canvas.Canvas, doc):
    c.saveState()
    c.setFont("Helvetica", 8)
    c.setFillColorRGB(0.45, 0.45, 0.45)
    c.drawString(2 * cm, 1.2 * cm, "EventCamArucoDetector  ·  algorithm explainer")
    c.drawRightString(A4[0] - 2 * cm, 1.2 * cm, f"page {doc.page}")
    c.restoreState()


# ---------------------------------------------------------------------------
# Story
# ---------------------------------------------------------------------------

def build_story() -> list:
    s: list = []

    # --- Title block ---
    title = Paragraph(
        "EventCamArucoDetector",
        ParagraphStyle("T", parent=H1, fontSize=26, spaceAfter=4,
                       textColor=colors.HexColor("#0d1b5a"),
                       alignment=1))
    subtitle = Paragraph(
        "Detecting ArUco markers from raw event-camera streams",
        ParagraphStyle("ST", parent=H2, fontSize=14,
                       textColor=colors.HexColor("#3b4cca"),
                       alignment=1))
    s += [Spacer(1, 4 * cm), title, subtitle, Spacer(1, 0.7 * cm)]

    s.append(Paragraph(
        "This document walks through the full pipeline used to find "
        "ArUco markers in the asynchronous output of an event camera. "
        "It covers every stage end-to-end: data layout, the multi-window "
        "sliding-tick scheme, quad detection, perspective rectification, "
        "scanline transition decoding, dictionary lookup, the optional "
        "marker-ID filter, result merging across runs, and the tabbed "
        "viewer that visualises the result tables.",
        ParagraphStyle("Lead", parent=BODY, fontSize=11.5,
                       textColor=colors.HexColor("#222"),
                       alignment=TA_JUSTIFY)))
    s.append(Spacer(1, 1.0 * cm))
    s.append(fig("01_pipeline.png", 14,
                 "Figure 1 - End-to-end pipeline. Each box maps directly "
                 "to a code unit in the repository."))
    s.append(PageBreak())

    # =====================================================================
    s.append(h1("1. Why event cameras need a different detector"))
    s.append(p(
        "A standard ArUco detector expects a global-shutter intensity image: "
        "find an external contour, fit a quadrilateral, threshold, decode. "
        "Event cameras don't produce frames at all. They emit an "
        "asynchronous stream of 4-tuples - "
        "<b>(x, y, polarity, timestamp)</b> - one per pixel whenever the "
        "log-intensity at that pixel changes by more than a contrast "
        "threshold. There is no concept of \"the picture at time t\"; you "
        "only ever see what changed."))
    s.append(p(
        "Static markers therefore produce <i>no</i> events. We rely on "
        "either the marker moving in the scene or the camera moving "
        "relative to it. When that happens, the high-contrast border of "
        "the marker triggers events at every crossing, which gives us "
        "edges to work with - but no fill, no greyscale and no "
        "guaranteed continuity."))

    s.append(h2("1.1  Input data layout"))
    s.append(p(
        "Each dataset is stored as a MATLAB <tt>.mat</tt> file containing "
        "a single Nx4 array called <tt>events</tt>:"))
    s.append(code(
        "events(:,1) = x        (column, 0-based)\n"
        "events(:,2) = y        (row,    0-based)\n"
        "events(:,3) = polarity (+1 on, 0/-1 off; ignored by the detector)\n"
        "events(:,4) = t        (timestamp in microseconds)"))
    s.append(p(
        "The DVS we used produces events on a 320x240 sensor "
        "(<tt>sensorSize = [240, 320]</tt>). At typical motion speeds this "
        "is around 1-5 million events per recording, sparse in space and "
        "dense in time."))

    # =====================================================================
    s.append(h1("2. Multi-window sliding-tick scheme"))
    s.append(p(
        "We do <i>not</i> try to track individual events. Instead we treat "
        "the stream as a piecewise time-integration problem. At every "
        "<b>tick</b> - by default every 1 ms - we ask the same question "
        "for several different lookback durations:"))
    s.append(b("How many events occurred in the window <tt>[tNow - dt, tNow]</tt>?"))
    s.append(b("Treat those events as a static binary image."))
    s.append(b("Run the rest of the pipeline on that image."))

    s.append(fig("02_sliding_tick.png", 15,
                 "Figure 2 - At each 1 ms tick we run several windows in "
                 "parallel. The set of windows is configured by "
                 "<tt>params.windowDurations_ms</tt>."))
    s.append(p(
        "Why several windows? Because the best integration time depends on "
        "the motion: a slow rotation needs hundreds of milliseconds to "
        "accumulate a clean edge, but for fast translation the same window "
        "will smear the marker into a blur. By trying multiple windows at "
        "the same tick and recording which (if any) succeeded, we get a "
        "robust detector and we can later analyse which window is best "
        "for which motion. That analysis is exactly what the result viewer "
        "exposes."))

    s.append(h2("2.1  Where do we start and stop?"))
    s.append(p(
        "Setting up the tick loop is one of the few places this code is "
        "fussy. We have to skip the first <i>max(windowDurations)</i> of "
        "the recording, because before that point no window has enough "
        "history. The relevant lines from <tt>detectAruco.m</tt>:"))
    s.append(code(
        "tMin   = evT(1);\n"
        "tMax   = evT(end);\n"
        "tStart = tMin + max(windowDurations_us);\n"
        "tEnd   = tMax;\n"
        "numTicks = floor((tEnd - tStart) / tickStep_us) + 1;"))
    s.append(p(
        "Note that running the same recording with a different window "
        "set therefore produces a different number of ticks - if we add a "
        "750 ms window we lose the first 750 ms of the recording. That is "
        "why <tt>mergeResults</tt> has to be careful about timelines."))

    # =====================================================================
    s.append(h1("3. From events to a binary image"))
    s.append(p(
        "Inside the window we have a list of <tt>(x, y)</tt> pixels that "
        "fired at least once. Two binary-search calls extract the right "
        "slice of the sorted <tt>evT</tt> array; <tt>accumarray</tt> then "
        "builds the per-pixel count image in one vectorised step:"))
    s.append(code(
        "iStart = bsearchLeft(evT, tNow - dt);\n"
        "iEnd   = bsearchRight(evT, tNow);\n"
        "wX = evX(iStart:iEnd);  wY = evY(iStart:iEnd);\n"
        "countImg   = accumarray([wY, wX], 1, [H, W]);\n"
        "activeMask = countImg > 0;"))
    s.append(p(
        "<tt>activeMask</tt> is what the quad detector consumes. "
        "<tt>countImg</tt> survives for later: when we want to <i>decode</i> "
        "a candidate marker we need greyscale-like intensity, so we use "
        "<tt>uint8(countImg / max(countImg) * 255)</tt> instead of the "
        "binary mask."))
    s.append(fig("03_event_accumulation.png", 15,
                 "Figure 3 - Left: the raw (x, y) of every event in the "
                 "window. Right: the same events accumulated into an HxW "
                 "count image."))

    s.append(h2("3.1  Sub-threshold rejection"))
    s.append(p(
        "If a window contains fewer than 10 events we skip it altogether - "
        "any quad we could fit would be noise. This is the cheapest "
        "early-exit in the loop."))

    # =====================================================================
    s.append(h1("4. Quad detection - two methods, then dedup"))
    s.append(p(
        "Event cameras give us only the marker's edges, never its fill. "
        "Sometimes the edges form a closed loop, sometimes there is a "
        "missing segment because the camera was momentarily still. We "
        "therefore run two complementary detectors on every frame and "
        "merge their hits."))

    s.append(fig("04_quad_detection.png", 16,
                 "Figure 4 - Method A fills closed boundaries and fits a "
                 "minimum-area rectangle to the filled blob. Method B "
                 "takes the convex hull of the raw edge component, which "
                 "survives a missing side."))

    s.append(h3("Method A - fill + watershed (closed-edge markers)"))
    s.append(b("<tt>imfill(activeMask, 'holes')</tt> closes the marker interior."))
    s.append(b("<tt>bwareaopen(bw, 100)</tt> drops noise blobs."))
    s.append(b("Distance transform + extended-maxima + watershed splits two "
               "markers that have touched."))
    s.append(b("<tt>findQuadCandidates</tt> fits a minimum-area rectangle to "
               "the convex hull of each surviving blob and accepts it if "
               "area, aspect, and rectangularity all pass."))

    s.append(h3("Method B - convex hull on raw edges (markers with gaps)"))
    s.append(b("Run <tt>bwconncomp</tt> directly on <tt>activeMask</tt> with "
               "4-connectivity (so two pixels need to share a side, not "
               "just a corner)."))
    s.append(b("For each component, take the convex hull and fit a "
               "minimum-area rectangle as before."))
    s.append(b("Because we work on the convex hull, a missing side does not "
               "invalidate the candidate."))

    s.append(h3("Deduplication"))
    s.append(p(
        "Both methods often return the same physical marker. "
        "<tt>deduplicateQuads</tt> walks the combined list and drops any "
        "quad whose centre is within half a side-length of an earlier "
        "one. The first-found wins."))

    s.append(h2("4.1  The filters that actually do the work"))
    s.append(p("Inside <tt>findQuadCandidates</tt> three predicates decide "
               "whether a region looks like an ArUco frame:"))
    s.append(b("<b>Area filter</b> - the <i>convex-hull</i> area must be "
               "between <tt>params.minArea</tt> (default 625 px) and "
               "<tt>params.maxArea</tt> (default 40% of the image)."))
    s.append(b("<b>Rectangularity</b> = convexArea / minRectArea must "
               "exceed 0.6. Genuine squares score close to 1; a noisy "
               "blob with a long tail scores much lower."))
    s.append(b("<b>Aspect ratio</b> of the minimum-area rectangle must "
               "be below 3.0. Elongated rectangles are not markers."))

    # =====================================================================
    s.append(h1("5. Perspective rectification"))
    s.append(p(
        "Each surviving quad gives us four image-space corners. We need "
        "to flatten that quad into a canonical 160 x 160 px square so the "
        "rest of the decoder can assume axis-aligned cells. We do this "
        "with a projective homography:"))
    s.append(code(
        "tform = fitgeotrans(srcCorners, dstCorners+1, 'projective');\n"
        "warpedImg = imwarp(countU8, tform, 'OutputView', ...);"))
    s.append(fig("05_unwarp.png", 15,
                 "Figure 5 - A distorted detected quad is mapped onto a "
                 "fixed 8x8 grid of 20-pixel cells (160x160 px)."))

    s.append(h2("5.1  Corner ordering"))
    s.append(p(
        "<tt>fitgeotrans</tt> only behaves correctly if the source and "
        "destination corners are in matching order. "
        "<tt>orderCornersForUnwarp_local</tt> sorts the quad by angle "
        "around its centroid, then rotates the list so the corner with "
        "the smallest x+y is first (top-left), and finally flips the "
        "order if the cross product disagrees with the canonical "
        "(TL, TR, BR, BL) winding. This makes the unwarp deterministic "
        "regardless of which corner the detector returned first."))

    # =====================================================================
    s.append(h1("6. Scanline transition decoding"))
    s.append(p(
        "Now we have a 160x160 image of a marker. A normal photo would "
        "let us threshold each cell directly to read 0/1. Event-derived "
        "images don't have a meaningful interior - the only signal is "
        "<i>at</i> the cell boundaries, where the polarity flips. "
        "We decode by following the scanlines and toggling colour at "
        "every boundary that crosses a learned threshold."))

    s.append(fig("06_scanline_decode.png", 17,
                 "Figure 6 - Cell boundaries (red) are where event "
                 "intensity peaks. Vertical and horizontal scanlines "
                 "(blue) sample column / row centres and flip the "
                 "current colour at each boundary that exceeds the "
                 "auto-tuned threshold."))

    s.append(h2("6.1  Boundary intensity, then Otsu threshold"))
    s.append(p(
        "For every cell boundary we sample the mean intensity in a "
        "narrow 11-pixel band (parameter <tt>boundaryHalfWidth = 5</tt>). "
        "The full set of boundary samples is then thresholded with "
        "<tt>graythresh</tt> (Otsu's method) to find the per-marker "
        "<tt>transThresh</tt>. This adapts to whatever brightness the "
        "current accumulated count image happens to have."))

    s.append(h2("6.2  Three candidate codes, six geometric variants"))
    s.append(p(
        "Reading the marker as a single vertical scan or a single "
        "horizontal scan can fail when one direction has a weak "
        "boundary; we therefore build <i>three</i> candidate 8x8 codes:"))
    s.append(b("<tt>codeV</tt> - scan each column top-to-bottom, flipping "
               "colour at every horizontal boundary above threshold."))
    s.append(b("<tt>codeH</tt> - scan each row left-to-right, flipping "
               "colour at every vertical boundary above threshold."))
    s.append(b("<tt>majority</tt> - bit-wise majority of <tt>codeV</tt> "
               "and <tt>codeH</tt>, which heals individual misreads."))

    s.append(fig("07_dictionary_lookup.png", 15,
                 "Figure 7 - For each of the three candidate codes, the "
                 "decoder tries all 2x2x4 = 16 combinations of invert, "
                 "horizontal flip, and 90-degree rotation before giving "
                 "up. Including the three candidates that is 48 lookups "
                 "per quad, each a single binary search."))

    s.append(h2("6.3  Inner code and dictionary lookup"))
    s.append(p(
        "After dropping the outer black border (rows/cols 1 and 8) we "
        "have a 6x6 grid of bits. We pack those 36 bits into a "
        "<tt>uint64</tt> with the dictionary's canonical bit order and "
        "binary-search the ARUCO_MIP_36h12 dictionary (250 codes), which "
        "is pre-sorted at startup. A hit returns the marker ID "
        "(0 to 249); a miss makes us try the next variant."))

    # =====================================================================
    s.append(h1("7. The requestedMarkerIds filter"))
    s.append(p(
        "Most experiments only care about <i>one</i> marker - the one "
        "printed in the recording. We added <tt>params.requestedMarkerIds</tt> "
        "so that decodes outside that set are treated as no-detection:"))
    s.append(code(
        "params.requestedMarkerIds = 3;          % only accept marker ID 3\n"
        "params.requestedMarkerIds = [3 7 12];   % accept any of these\n"
        "params.requestedMarkerIds = [];          % accept any decoded ID"))
    s.append(p(
        "Concretely, inside the quad loop the check changes from"))
    s.append(code(
        "if mid &gt;= 0\n"
        "    bestID = mid; break;\n"
        "end"))
    s.append(p("to"))
    s.append(code(
        "if mid &gt;= 0 &amp;&amp; (isempty(requestedMarkerIds) || ...\n"
        "                       any(requestedMarkerIds == int32(mid)))\n"
        "    bestID = mid; break;\n"
        "end"))
    s.append(p(
        "The decoded ID is still stored per window (in the <tt>win_Xms</tt> "
        "columns) so that the viewer can slice the results by marker after "
        "the fact - even when several IDs are recorded inside the same "
        "run, you can still ask \"how often did marker 7 show up?\""))

    # =====================================================================
    s.append(h1("8. Output schema"))
    s.append(p("Every detection run writes a struct with these fields:"))
    tbl = Table([
        ["field", "shape", "meaning"],
        ["tNow_us",              "Nx1 double", "timestamp of each tick (microseconds)"],
        ["anyDetected",          "Nx1 double", "1 if any window detected this tick, else 0"],
        ["win_Xms",              "Nx1 double", "decoded marker ID for window X ms (-1 = none)"],
        ["windowDurations_ms",   "1xW double", "the window set used"],
        ["detectionsPerWindow",  "1xW double", "total hits per window across all ticks"],
        ["requestedMarkerIds",   "1xK double", "filter that was applied ([] = any)"],
        ["attemptedPerWindow",   "1xW double", "(after merge) how many ticks each window was actually run on"],
    ], colWidths=[4 * cm, 3 * cm, 8.5 * cm])
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dfe7fd")),
        ("FONT", (0, 0), (-1, 0), "Helvetica-Bold", 10),
        ("FONT", (0, 1), (-1, -1), "Helvetica", 9.5),
        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#aab")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.white, colors.HexColor("#f6f8fc")]),
    ]))
    s.append(tbl)
    s.append(Spacer(1, 0.4 * cm))

    s.append(h2("8.1  Merging multiple runs over the same input"))
    s.append(p(
        "Running detection twice on the same recording with different "
        "window sets gives us two result files with <i>different</i> "
        "tick counts (because <tt>tStart</tt> depends on "
        "<tt>max(windowDurations)</tt>). <tt>mergeResults</tt> joins "
        "them on the union of timestamps, fills missing entries with "
        "-1, and records <tt>attemptedPerWindow</tt> so we can report "
        "an honest per-window rate."))
    s.append(fig("08_merge_union.png", 15,
                 "Figure 8 - Union-of-timestamps merge. v1 covers the "
                 "full recording but only with small windows; v2 skips "
                 "the first 600 ticks but adds large windows. The merged "
                 "file keeps every tick and is honest about which window "
                 "was attempted where."))

    s.append(p(
        "There are two natural denominators after a merge, and "
        "<tt>viewAllResults</tt> uses both:"))
    s.append(b("<b>AnyDetectPct</b> = #ticks with any successful window / "
               "<i>union length</i>. This is the right number for \"how "
               "often did we see the marker at all?\""))
    s.append(b("<b>BestRatePct</b> = #detections of that window / "
               "<i>attemptedPerWindow</i> of that window. This is the "
               "right number for \"when this window was tried, how often "
               "did it succeed?\""))
    s.append(p(
        "Because the denominators differ, you can legitimately see "
        "AnyDetectPct = 95.3% next to a 500 ms window at 100% - the "
        "first 600 ticks were only covered by smaller windows in v1, "
        "which missed some, while the 500 ms window from v2 nailed "
        "everything it was actually tried on."))

    # =====================================================================
    s.append(h1("9. The result viewer"))
    s.append(p(
        "<tt>viewAllResults('Data')</tt> opens a single window with one "
        "summary tab and one detail tab per dataset. A marker-filter "
        "dropdown lets you recompute every tab for a specific decoded "
        "ID without re-running detection. Three save buttons export the "
        "current tab, the current window, or all tabs as a multi-page "
        "PDF."))
    s.append(fig("09_viewer_mock.png", 16,
                 "Figure 9 - The summary tab. The dropdown switches the "
                 "view between \"all markers\" and any specific ID found "
                 "in the data."))
    s.append(fig("10_detail_mock.png", 16,
                 "Figure 10 - A per-dataset detail tab. The raster shows "
                 "when each window fired; the bar chart shows the "
                 "per-window success rate using attemptedPerWindow as "
                 "denominator."))

    # =====================================================================
    s.append(h1("10. Parallelism and performance notes"))
    s.append(b("The tick loop is <b>embarrassingly parallel</b>; each tick "
               "is independent. With Parallel Computing Toolbox and "
               "<tt>params.useParallel = true</tt> we use <tt>parfor</tt> "
               "and a <tt>parallel.pool.DataQueue</tt> for live "
               "progress."))
    s.append(b("Without the toolbox (or when "
               "<tt>params.useParallel = false</tt>) the same loop runs "
               "sequentially with a 1%-resolution progress print."))
    s.append(b("Dictionary lookup is <tt>O(log 250)</tt> per quad per "
               "variant - <tt>buildDictionaryArrays</tt> sorts the codes "
               "once so the worker copies are immutable."))
    s.append(b("<tt>accumarray</tt> is the only big-O hotspot inside the "
               "inner loop and is fully vectorised."))

    s.append(h2("10.1  What's the bottleneck?"))
    s.append(p(
        "Empirically the dominant cost is <tt>imwarp</tt>, not decoding. "
        "If a recording has many false-positive quads, we pay one "
        "perspective warp per quad per window per tick. The "
        "<tt>minRectangularity</tt> and <tt>maxAspect</tt> filters in "
        "<tt>findQuadCandidates</tt> are the main lever for keeping that "
        "in check."))

    # =====================================================================
    s.append(h1("11. Where the bodies are buried"))
    s.append(p("Things worth remembering when you read or modify the code:"))
    s.append(b("<tt>events(:,3)</tt> (polarity) is ignored. The detector "
               "only uses positions and timestamps."))
    s.append(b("Marker decoding tries 48 variants per candidate. If you "
               "speed-tune this, that is the loop to look at."))
    s.append(b("<tt>tStart</tt> depends on the <i>largest</i> requested "
               "window. Changing the window set shifts the start, so "
               "tick counts change between runs and have to be merged on "
               "the union of timestamps."))
    s.append(b("With a non-empty <tt>requestedMarkerIds</tt>, decodes "
               "outside the set are stored as -1 in the relevant cell. "
               "This is silent: nothing in the result file tells you "
               "\"by the way, I rejected three valid marker-9 hits "
               "here\". If that matters, log them in "
               "<tt>detectAruco</tt>."))
    s.append(b("Merging an unfiltered v1 with a filtered v2 produces a "
               "file whose semantics differ across windows. The viewer's "
               "marker dropdown will surface IDs that exist only in some "
               "windows. Keep the filter consistent across runs you "
               "intend to merge."))

    # =====================================================================
    s.append(h1("12. File map"))
    file_tbl = Table([
        ["file", "responsibility"],
        ["main.m",                "configure inputs / parameters, drive the per-dataset loop, kick off merge + viewer"],
        ["detectAruco.m",         "tick loop, event accumulation, decoder, the parfor/sequential split"],
        ["detectQuadBlob.m",      "Method A (fill+watershed) + Method B (convex hull on edges) + dedup"],
        ["findQuadCandidates.m",  "area / rectangularity / aspect filters, minimum-area rect fit, corner ordering"],
        ["Utils/mergeResults.m",   "union-timeline merge of two or more result files"],
        ["Utils/mergeAllResults.m","batch-merge every v1/v2 pair under Data/"],
        ["viewAllResults.m",      "tabbed GUI, marker-filter dropdown, PDF exports"],
        ["analyzeResults.m",      "single-result deep-dive: raster, bars, heatmap, gap distribution, rolling rate"],
    ], colWidths=[5 * cm, 11 * cm])
    file_tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dfe7fd")),
        ("FONT", (0, 0), (-1, 0), "Helvetica-Bold", 10),
        ("FONT", (0, 1), (-1, -1), "Helvetica", 9.5),
        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#aab")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.white, colors.HexColor("#f6f8fc")]),
    ]))
    s.append(file_tbl)

    return s


def main():
    doc = SimpleDocTemplate(
        str(OUT_PDF),
        pagesize=A4,
        leftMargin=2.0 * cm, rightMargin=2.0 * cm,
        topMargin=1.6 * cm, bottomMargin=1.8 * cm,
        title="EventCamArucoDetector - Algorithm Explanation",
        author="EventCamArucoDetector",
    )
    story = build_story()
    doc.build(story, onFirstPage=_draw_footer, onLaterPages=_draw_footer)
    print(f"wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
