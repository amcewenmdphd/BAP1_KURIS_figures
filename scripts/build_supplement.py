#!/usr/bin/env python3
"""Assemble the final, numbered supplement.

For each supplementary figure, compose a US-Letter PDF with the figure placed at
double-column width at the top of the page and its caption below — on the same
page when there is room, otherwise on the following page (Genetics in Medicine
convention). Writes supplement/Supplementary_Figure_N.pdf.

The numbered, captioned supplementary *tables* are written by the individual
table scripts (via bap1figs.supplement.render_supp_table); this script also
writes the reader's guide (supplement/SUPPLEMENT_GUIDE.md) that lists every
figure and table with its caption.

Run the table scripts (or scripts/run_all.py) first so the source figures and
supplement tables exist.
"""

import sys
from pathlib import Path

import fitz  # PyMuPDF

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bap1figs import config, supplement, richtext

LETTER = (612.0, 792.0)      # points
MARGIN = 54.0                # 0.75 inch
CONTENT_W = LETTER[0] - 2 * MARGIN
CONTENT_H = LETTER[1] - 2 * MARGIN

TITLE_FS, CAP_FS, LEAD = 9.5, 9.0, 1.32
FIG_GAP = 16.0               # points between figure and caption
# Captions use DejaVu Sans (the figures' own font) via bap1figs.richtext, which
# also italicizes gene symbols and "et al." and renders the Unicode (minus sign,
# curly quotes, ≤ / ≥) that base-14 Helvetica can't.


def _caption_height(fig):
    title_lines = richtext.layout(richtext.normalize(fig.title), "bld", TITLE_FS, CONTENT_W)
    cap_lines = richtext.layout(richtext.normalize(fig.caption), "reg", CAP_FS, CONTENT_W)
    h = len(title_lines) * TITLE_FS * LEAD + 4.0 + len(cap_lines) * CAP_FS * LEAD
    return h, title_lines, cap_lines


def _draw_caption(page, top, title_lines, cap_lines):
    tw_title = fitz.TextWriter(page.rect)
    tw_body = fitz.TextWriter(page.rect)
    y = top
    for line in title_lines:
        y += TITLE_FS * LEAD
        for s, it, x in line:
            tw_title.append((MARGIN + x, y), s, font=richtext.font("bld", it), fontsize=TITLE_FS)
    y += 4.0
    for line in cap_lines:
        y += CAP_FS * LEAD
        for s, it, x in line:
            tw_body.append((MARGIN + x, y), s, font=richtext.font("reg", it), fontsize=CAP_FS)
    tw_title.write_text(page, color=(0, 0, 0))
    tw_body.write_text(page, color=(0.15, 0.15, 0.15))


def compose_figure(fig):
    src_path = config.FIGURES / f"{fig.stem}.pdf"
    if not src_path.exists():
        print(f"WARN: missing {src_path.relative_to(config.REPO)}; skipping {fig.label}")
        return None
    src = fitz.open(src_path)
    srect = src[0].rect
    aspect = srect.width / srect.height           # w / h

    # Fit to double-column width, but never taller than the content area.
    w = CONTENT_W
    h = w / aspect
    if h > CONTENT_H:
        h = CONTENT_H
        w = h * aspect
    fx = MARGIN + (CONTENT_W - w) / 2.0            # centre horizontally

    cap_h, title_lines, cap_lines = _caption_height(fig)
    same_page = (h + FIG_GAP + cap_h) <= CONTENT_H

    out = fitz.open()
    page = out.new_page(width=LETTER[0], height=LETTER[1])
    fig_rect = fitz.Rect(fx, MARGIN, fx + w, MARGIN + h)
    page.show_pdf_page(fig_rect, src, 0)          # vector-embed, stays crisp
    if same_page:
        _draw_caption(page, MARGIN + h + FIG_GAP, title_lines, cap_lines)
    else:
        cpage = out.new_page(width=LETTER[0], height=LETTER[1])
        _draw_caption(cpage, MARGIN, title_lines, cap_lines)

    out_path = config.SUPPLEMENT / f"{fig.out_name}.pdf"
    out.save(out_path)
    out.close(); src.close()
    tag = "caption on same page" if same_page else "caption on following page"
    print(f"INFO: wrote {out_path.relative_to(config.REPO)} ({tag})")
    return out_path


def write_guide():
    lines = [
        "# Supplementary material — reader's guide",
        "",
        "BAP1 saturation genome editing functional evidence for KURIS/NDD variant "
        "interpretation. Figures and tables are numbered in the order they are cited; "
        "each file is a self-contained US-Letter PDF (figures carry their caption, "
        "tables carry their title and caption). Working copies of every table are also "
        "provided as `.tsv` and `.xlsx` in `tables/`.",
        "",
        "## Supplementary figures",
        "",
    ]
    for f in supplement.SUPP_FIGURES:
        lines += [f"**{f.title}**", "", f"{f.caption}", "",
                  f"_File:_ `supplement/{f.out_name}.pdf`", ""]
    lines += ["## Supplementary tables", ""]
    for t in supplement.SUPP_TABLES:
        lines += [f"**{t.title}**", "", f"{t.caption}", "",
                  f"_File:_ `supplement/{t.out_name}.pdf` "
                  f"(data: `tables/{t.key}.tsv`, `tables/{t.key}.xlsx`)", ""]
    guide = config.SUPPLEMENT / "SUPPLEMENT_GUIDE.md"
    guide.write_text("\n".join(lines))
    print(f"INFO: wrote {guide.relative_to(config.REPO)}")


def main():
    for fig in supplement.SUPP_FIGURES:
        compose_figure(fig)
    write_guide()


if __name__ == "__main__":
    main()
