#!/usr/bin/env python3
"""Generate the NeuroRDViz setup report as a PDF."""

from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_LEFT
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Preformatted,
)

OUTPUT = "/Users/lg002/NeuroRDViz/NeuroRDViz_Setup_Report.pdf"

doc = SimpleDocTemplate(
    OUTPUT,
    pagesize=letter,
    topMargin=0.75 * inch,
    bottomMargin=0.75 * inch,
    leftMargin=0.75 * inch,
    rightMargin=0.75 * inch,
)

styles = getSampleStyleSheet()
title_style = styles["Title"]
h1 = styles["Heading1"]
h2 = styles["Heading2"]
body = styles["BodyText"]
body_bold = ParagraphStyle("BodyBold", parent=body, fontName="Helvetica-Bold")
code_style = ParagraphStyle(
    "Code",
    parent=body,
    fontName="Courier",
    fontSize=8.5,
    leading=11,
    leftIndent=18,
    backColor=HexColor("#f0f0f0"),
    borderPadding=4,
)
bullet = ParagraphStyle("Bullet", parent=body, leftIndent=24, bulletIndent=12)

story = []

# ── Title ────────────────────────────────────────────────────────────
story.append(Paragraph("NeuroRDViz Modernization — Setup Report", title_style))
story.append(Spacer(1, 12))

# ── 1  Repo & Code ──────────────────────────────────────────────────
story.append(Paragraph("1. Repository Cloned and Updated", h1))
story.append(Paragraph(
    "Cloned the existing repo from <b>https://github.com/neurord/NeuroRDViz.git</b> "
    "into <font face='Courier'>/Users/lg002/NeuroRDViz</font>.", body))
story.append(Paragraph(
    "Replaced the old PyQt4/Python 2 <font face='Courier'>NeuroRDViz.py</font> with the "
    "modernized PySide6 version (from <font face='Courier'>NeuroRDViz_Working_Morphology.txt</font>).", body))
story.append(Paragraph(
    "Copied the test data file <font face='Courier'>Model_mglur_spine-CaN.h5</font> into the project directory.", body))
story.append(Spacer(1, 8))

# ── 2  Key Code Changes ─────────────────────────────────────────────
story.append(Paragraph("2. Key Changes in the Modernized Code", h1))

changes = [
    ("<b>Qt framework</b>: PyQt4 → PySide6 (Qt6). Compatibility aliases map old "
     "<font face='Courier'>QtGui.*</font> widget references to their new "
     "<font face='Courier'>QtWidgets.*</font> locations."),
    ("<b>HDF5 bytes handling</b>: All molecule name lookups normalize bytes→str consistently, "
     "since h5py returns <font face='Courier'>bytes</font> for string datasets. "
     "The old code compared raw bytes directly, which breaks under Python 3."),
    ("<b>Output set iteration</b>: <font face='Courier'>get_mol_info()</font> now iterates "
     "all output sets (not just <font face='Courier'>outputsets[1:]</font> and hardcoded keys), "
     "handles varying time-step counts across sets, and picks the maximum sample count for "
     "molecules that appear in multiple output sets."),
    ("<b>Voxel volume extraction</b>: <font face='Courier'>get_voxel_volumes()</font> robustly "
     "handles the HDF5 grid — tries structured array field 'volume' first, then falls back "
     "to a separate dataset, then to uniform 1.0."),
    ("<b>Time-step mismatch handling</b>: <font face='Courier'>get_voxel_molecule_conc()</font> "
     "defensively handles cases where output sets have different numbers of time points "
     "(truncates or pads)."),
    ("<b>Animation loop</b>: The <font face='Courier'>anim()</font> generator now updates via "
     "<font face='Courier'>surf.mlab_source.set(scalars=...)</font> (preferred Mayavi update path) "
     "with a fallback to direct unstructured grid modification. Slider updates block signals to "
     "prevent feedback loops."),
    ("<b>Signal connections</b>: <font face='Courier'>comboBox.activated[str]</font> replaced with "
     "<font face='Courier'>comboBox.textActivated</font> (PySide6-native signal)."),
    ("<b>Dialog API updates</b>: <font face='Courier'>QMessageBox.question()</font> uses "
     "<font face='Courier'>StandardButton</font> enum flags instead of positional string arguments."),
]
for c in changes:
    story.append(Paragraph(c, bullet, bulletText="\u2022"))
story.append(Spacer(1, 8))

# ── 3  Environment ───────────────────────────────────────────────────
story.append(Paragraph("3. Environment Setup", h1))
story.append(Paragraph(
    "Creating a working Python environment on macOS ARM64 (Apple Silicon) required "
    "navigating several dependency conflicts:", body))
story.append(Spacer(1, 6))

table_data = [
    ["Attempt", "Approach", "Failure"],
    ["1", "pip install everything",
     "Mayavi wheel build crashes — pip's VTK wheels trigger an IOKit "
     "ENABLE_FIELD_RECOGNITION error on macOS ARM during TVTK class generation"],
    ["2", "Conda-forge for everything at once",
     "Solver took 10+ minutes and hung"],
    ["3", "Conda-forge VTK + pip Mayavi\n+ pip PySide6",
     "Qt library conflict — conda VTK ships its own libQt6Core.dylib, "
     "clashing with pip PySide6's bundled copy"],
    ["4", "Conda-forge VTK + conda-forge\nPySide6 + pip Mayavi",
     "TVTK/VTK version mismatch — Mayavi 4.8.3's generated tvtk wrapper "
     "classes don't handle VTK 9.5.2's new properties"],
    ["5 \u2713", "All from conda-forge:\nmayavi, h5py, pyside6",
     "Worked, but required pinning numpy<2 (see below)"],
]
t = Table(table_data, colWidths=[0.6 * inch, 2.0 * inch, 4.2 * inch])
t.setStyle(TableStyle([
    ("BACKGROUND", (0, 0), (-1, 0), HexColor("#4472C4")),
    ("TEXTCOLOR", (0, 0), (-1, 0), HexColor("#FFFFFF")),
    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
    ("FONTSIZE", (0, 0), (-1, -1), 8),
    ("LEADING", (0, 0), (-1, -1), 10),
    ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ("GRID", (0, 0), (-1, -1), 0.5, HexColor("#CCCCCC")),
    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [HexColor("#FFFFFF"), HexColor("#F2F2F2")]),
    ("BACKGROUND", (0, 5), (-1, 5), HexColor("#E2EFDA")),
    ("LEFTPADDING", (0, 0), (-1, -1), 4),
    ("RIGHTPADDING", (0, 0), (-1, -1), 4),
    ("TOPPADDING", (0, 0), (-1, -1), 3),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
]))
story.append(t)
story.append(Spacer(1, 8))

story.append(Paragraph(
    "<b>The numpy&lt;2 pin</b>: VTK 9.4.2's <font face='Courier'>numpy_interface/algorithms.py</font> "
    "references <font face='Courier'>numpy.in1d</font>, which was removed in numpy 2.0. "
    "Pinning to numpy 1.26.4 resolves this. This is a known upstream issue; future VTK releases will fix it.",
    body))
story.append(Spacer(1, 6))

story.append(Paragraph("<b>Final working environment:</b>", body))
story.append(Preformatted(
    "conda create -n neurordviz --strict-channel-priority \\\n"
    "    -c conda-forge python=3.11 mayavi h5py pyside6\n"
    "conda install -n neurordviz -c conda-forge \"numpy<2\"",
    code_style))
story.append(Spacer(1, 6))

ver_data = [
    ["Package", "Version"],
    ["Python", "3.11"],
    ["VTK", "9.4.2"],
    ["Mayavi", "4.8.3"],
    ["PySide6", "6.9.2"],
    ["numpy", "1.26.4"],
    ["h5py", "3.15.1"],
    ["traits", "7.1.0"],
    ["traitsui", "8.0.0"],
]
vt = Table(ver_data, colWidths=[1.5 * inch, 1.5 * inch])
vt.setStyle(TableStyle([
    ("BACKGROUND", (0, 0), (-1, 0), HexColor("#4472C4")),
    ("TEXTCOLOR", (0, 0), (-1, 0), HexColor("#FFFFFF")),
    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
    ("FONTSIZE", (0, 0), (-1, -1), 9),
    ("GRID", (0, 0), (-1, -1), 0.5, HexColor("#CCCCCC")),
    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [HexColor("#FFFFFF"), HexColor("#F2F2F2")]),
    ("LEFTPADDING", (0, 0), (-1, -1), 6),
    ("TOPPADDING", (0, 0), (-1, -1), 3),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
]))
story.append(vt)
story.append(Spacer(1, 8))

# ── 4  Verified ──────────────────────────────────────────────────────
story.append(Paragraph("4. Verified Working", h1))
story.append(Preformatted(
    "conda activate neurordviz\n"
    "cd /Users/lg002/NeuroRDViz\n"
    "python NeuroRDViz.py Model_mglur_spine-CaN.h5",
    code_style))
story.append(Spacer(1, 6))
story.append(Paragraph(
    "The application launches, renders the 36-voxel spine morphology as hexahedral mesh elements "
    "(432 edges), populates the molecule dropdown with all 20 species, and is ready for "
    "concentration animation on molecule selection.", body))
story.append(Spacer(1, 8))

# ── 5  New UI Features ──────────────────────────────────────────────
story.append(Paragraph("5. New UI Features", h1))

story.append(Paragraph("<b>Global Colorbar Range</b>", body_bold))
story.append(Paragraph(
    "The concentration colorbar now shows the global min/max across the <i>entire</i> animation, "
    "not just the current frame. This is locked before the first frame renders and re-locked after "
    "every frame update so the color scale remains stable throughout playback.", body))
story.append(Spacer(1, 6))

story.append(Paragraph("<b>Multi-Viewer Architecture</b>", body_bold))
story.append(Paragraph(
    "Multiple viewers can run independent animations simultaneously. Each viewer has its own "
    "Mayavi scene, molecule dropdown, start/stop buttons, progress bar, and scrub slider.", body))
story.append(Spacer(1, 4))
story.append(Paragraph(
    "Getting this working required solving a deep Mayavi scoping bug. Both the global "
    "<font face='Courier'>mlab</font> module and the scene-specific "
    "<font face='Courier'>viewer.visualization.scene.mlab</font> route pipeline operations "
    "through Mayavi's global engine, which adds new objects to whichever scene was last activated — "
    "always the newest viewer. This caused all animation surfaces to render in the most recently "
    "created viewer regardless of which viewer's controls were used.", body))
story.append(Spacer(1, 4))
story.append(Paragraph(
    "The fix: <b>never create new pipeline objects during animation</b>. "
    "<font face='Courier'>Visualization.update_plot()</font> (fired on "
    "<font face='Courier'>scene.activated</font>, when the scene IS the active one) now stores "
    "<font face='Courier'>self.ug</font> and <font face='Courier'>self.surf</font>. Animation "
    "setup configures the <i>existing</i> surface's LUT (colormap, range, scalar bar) and sets "
    "concentration scalars directly on the existing unstructured grid. Each frame update modifies "
    "<font face='Courier'>ug.point_data.scalars</font> in-place and calls "
    "<font face='Courier'>scene.render()</font> on the correct viewer's scene. A per-viewer "
    "<font face='Courier'>QTimer</font> drives each animation independently, replacing the old "
    "<font face='Courier'>@mlab.animate</font> decorator (which also used global figure routing).", body))
story.append(Spacer(1, 6))

story.append(Paragraph("<b>Animation Controls</b>", body_bold))
ui_features = [
    "<b>Per-viewer Start/Stop buttons</b>: Start or stop each viewer's animation independently. "
    "Starting viewers at different times produces desynchronized animations.",
    "<b>Start All / Stop All</b>: Global buttons (also in File menu and toolbar) that launch or halt "
    "all viewer animations simultaneously. 'Start All' resets all frames to 0 first, ensuring perfect sync.",
    "<b>Per-viewer progress bar + scrub slider</b>: Each viewer shows its own progress bar, time label, "
    "and horizontal slider beneath the 3D view. Dragging the slider seeks that viewer independently.",
    "<b>Global scrub slider</b>: A slider at the bottom of the window that, when dragged, moves "
    "<i>all</i> viewers to the same percentage point — useful for comparing molecules at the same time offset.",
    "<b>+ Add Viewer button</b>: Adds a new viewer from the bottom of the window without using the menu.",
]
for f in ui_features:
    story.append(Paragraph(f, bullet, bulletText="\u2022"))

# ── Build ────────────────────────────────────────────────────────────
doc.build(story)
print(f"PDF written to {OUTPUT}")
