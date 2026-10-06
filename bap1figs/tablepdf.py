"""Render a DataFrame to a publication-style PDF table (one helper shared by
every supplementary-table script, so they all look identical).

Booktabs styling: three horizontal rules (top / below header / bottom), no
vertical lines, no shading — the convention in most journals. Columns are fitted
to their content; wide free-text columns wrap; long tables paginate. Pages are
real US-Letter, portrait unless a table is too wide to stay legible. Optional
light rules separate row groups. Footnotes sit below the bottom rule.
"""

from functools import lru_cache

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath

from . import config, richtext

RULE = "#1A1A1A"          # top/bottom/header rules
GROUP_RULE = "#BFBFBF"    # faint separators between row groups
TEXT = "#1A1A1A"
DETAIL = "#333333"        # full-width detail (clinical occurrence) text
FOOT = "#444444"

LETTER = (8.5, 11.0)
MARGIN = 0.55           # inches
CELL_PAD = 0.22         # inches of horizontal padding per column (total)
MIN_COL = 0.5           # inches, smallest a column may be squeezed to
MAX_COL = 2.6           # inches, widest before a column is wrapped
LEADING = 1.65          # line height as a multiple of font size (row airiness)


def _fmt(v):
    if v is None or (isinstance(v, float) and v != v):     # None / NaN
        return ""
    if isinstance(v, float):
        return f"{v:.4f}" if abs(v) < 1000 else f"{v:.1f}"
    return str(v)


_FONT = {(bold, italic): FontProperties(
             weight="bold" if bold else "normal",
             style="italic" if italic else "normal")
         for bold in (False, True) for italic in (False, True)}


@lru_cache(maxsize=None)
def _glyph_w(tok, fontsize, bold, italic=False):
    """Ink-bbox width (inches) of a single space-free token — the real glyph
    measurement (matplotlib TextPath), memoized since the same tokens (ClinVar
    boilerplate words, amino-acid codes, HGVS fragments) recur constantly
    across a table, and re-building a font glyph path per call is the
    expensive part of this. Returns 0 for an empty token."""
    if not tok:
        return 0.0
    prop = _FONT[(bold, italic)]
    return TextPath((0, 0), tok, size=fontsize, prop=prop).get_extents().width / 72.0


@lru_cache(maxsize=None)
def _space_w(fontsize, bold, italic=False):
    """Width (inches) of a single space at this fontsize/weight. A space has no
    ink of its own, so TextPath's ink-bbox can't measure it directly (a
    whitespace-only string reports an empty/-inf extent); derive it instead
    from the advance-width difference between 'i i' and two 'i's."""
    return max(0.0, _glyph_w("i i", fontsize, bold, italic)
               - 2 * _glyph_w("i", fontsize, bold, italic))


def _text_w(s, fontsize, bold=False, italic=False):
    """Rendered width of s, in inches — real glyph metrics (matplotlib TextPath),
    not approximated by character count, so wrap points and underlines land
    exactly where the text actually ends. Decomposed into cached per-word
    widths (+ the derived space width) rather than one TextPath call on the
    whole string, so repeated words/strings across a table are measured once."""
    if not s:
        return 0.0
    if " " not in s:
        return _glyph_w(s, fontsize, bold, italic)
    words = s.split(" ")
    return (sum(_glyph_w(w, fontsize, bold, italic) for w in words)
            + _space_w(fontsize, bold, italic) * (len(words) - 1))


_TOKEN_BREAK_CHARS = (":", ">")

# Tolerance (inches) for "does this text fit in this width" checks. A width
# computed via one arithmetic path (e.g. natural-width-minus-padding, later
# uniformly rescaled) and a width re-measured via a fresh _text_w call at
# render time are mathematically equal but not always bit-identical —
# sub-ULP float noise can otherwise make a value compare as "1 point wider
# than its own column" and trigger a pointless hard split (e.g. '184' ->
# '18'/'4'). This is many orders of magnitude below anything visible.
_FIT_EPS = 1e-6


def _split_token(tok, width_in, fontsize, bold=False):
    """Break a single unbreakable token (no spaces) wider than width_in: prefer
    splitting right after the last natural delimiter (':' / '>', as in HGVS
    strings like 'NM_004656.4:c.687C>A') whose left fragment still fits,
    keeping the delimiter with that fragment; fall back to a hard character
    split (longest prefix that fits) only if no such delimiter fits. Without
    this, a token wider than its (possibly squeezed) column renders unbroken
    and bleeds into the next column instead of wrapping."""
    if _text_w(tok, fontsize, bold) <= width_in + _FIT_EPS:
        return [tok]
    cut = None
    for i, ch in enumerate(tok):
        if ch in _TOKEN_BREAK_CHARS and _text_w(tok[:i + 1], fontsize, bold) <= width_in + _FIT_EPS:
            cut = i + 1
    if cut and cut < len(tok):
        return [tok[:cut]] + _split_token(tok[cut:], width_in, fontsize, bold)
    lo, hi = 1, len(tok)                       # longest hard-split prefix that fits
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if _text_w(tok[:mid], fontsize, bold) <= width_in + _FIT_EPS:
            lo = mid
        else:
            hi = mid - 1
    lo = max(1, lo)
    rest = _split_token(tok[lo:], width_in, fontsize, bold) if lo < len(tok) else []
    return [tok[:lo]] + rest


def _wrap(s, width_in, fontsize, bold=False):
    """Word-wrap s to fit within width_in (inches) at fontsize, measuring each
    candidate line's real rendered width rather than approximating by
    character count — so lines use the full available width (no wasted
    whitespace) and a squeezed column never overflows into its neighbor. A
    single word wider than width_in is broken via _split_token."""
    words = s.split(" ")
    lines, cur = [], ""
    for w in words:
        cand = w if not cur else f"{cur} {w}"
        if not cur or _text_w(cand, fontsize, bold) <= width_in + _FIT_EPS:
            cur = cand
        else:
            lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)
    out = []
    for line in lines:
        out.extend(_split_token(line, width_in, fontsize, bold)
                    if _text_w(line, fontsize, bold) > width_in + _FIT_EPS else [line])
    return out


def _layout_rich(text, width_in, fontsize, bold=False):
    """Word-wrap text to width_in (inches) at fontsize, italicizing gene
    symbols (bap1figs.richtext.word_runs), measured with _text_w — the same
    matplotlib glyph metrics used to draw it, so positions land exactly where
    the words end instead of drifting (see the note on richtext.word_runs).
    Returns a list of lines; each line is a list of (subtext, italic,
    x_offset_in) with x_offset_in the left offset of the run within the line,
    in inches."""
    space_w = _space_w(fontsize, bold)
    lines, cur, x = [], [], 0.0
    for word in richtext.word_runs(text):
        ww = sum(_text_w(s, fontsize, bold=bold, italic=it) for s, it in word)
        if cur and x + ww > width_in + _FIT_EPS:
            lines.append(cur)
            cur, x = [], 0.0
        wx = x
        for s, it in word:
            cur.append((s, it, wx))
            wx += _text_w(s, fontsize, bold=bold, italic=it)
        x = wx + space_w
    if cur:
        lines.append(cur)
    return lines


def _word_floor(cells_list, fontsize, bold=False):
    """Width (inches) of the longest single space-delimited word among a
    column's header + cell values, + CELL_PAD — the true floor that column
    can be squeezed to without _wrap's hard-split fallback ever kicking in.
    A value with no space at all (a formatted CI like '(0.0028-0.0441)', a
    camelCase method name like 'StandardizedClass') has no safe squeeze point
    smaller than its own full width."""
    w = 0.0
    for s in cells_list:
        if not s:
            continue
        for word in s.split(" "):
            w = max(w, _text_w(word, fontsize, bold=bold))
    return w + CELL_PAD


def _fit(cols, natural, floor, page, fontsize, no_squeeze=frozenset()):
    """Fit columns to a page: wrap the widest columns, then shrink the font only
    as a last resort. No column is squeezed past its own word-floor (see
    _word_floor) or past MIN_COL, whichever is larger; no_squeeze columns are
    never chosen at all (their whole content, not just one word, must stay
    unbroken — e.g. because it's meant to read as one line). The width a
    protected column would have given up is taken from other columns instead,
    or from the final uniform font shrink, which keeps every column's
    proportions (so nothing that fit before stops fitting). Returns (widths,
    wrap_width_in_by_col, fontsize)."""
    avail_w = page[0] - 2 * MARGIN
    widths = list(natural)
    # wrapw always holds "current width minus padding" for every column, not
    # just squeezed ones — both so the final uniform shrink (below) can rescale
    # it correctly (CELL_PAD doesn't shrink on its own, so a column that was
    # never squeeze-selected would otherwise keep its pre-shrink, too-generous
    # wrap width and never re-check whether short header text like "FN" still
    # fits at the smaller fontsize) and so header wrapping can reuse it
    # directly instead of recomputing "width - CELL_PAD" with a stale constant.
    wrapw = {c: w - CELL_PAD for c, w in zip(cols, widths)}
    while sum(widths) > avail_w:
        idx = max((i for i in range(len(cols))
                    if widths[i] > max(MIN_COL, floor[i]) + 1e-6
                    and cols[i] not in no_squeeze),
                  key=lambda i: widths[i], default=None)
        if idx is None:
            break
        widths[idx] = max(MIN_COL, floor[idx], widths[idx] * 0.85)
        wrapw[cols[idx]] = widths[idx] - CELL_PAD
    if sum(widths) > avail_w:
        scale = avail_w / sum(widths)
        fontsize *= scale
        widths = [w * scale for w in widths]
        wrapw = {c: w * scale for c, w in wrapw.items()}
    return widths, wrapw, fontsize


def render_table_pdf(df, out_name, aligns=None, star_note=None, note=None,
                     group_cols=None, orient=None, detail_cols=None, fontsize=7.5,
                     title=None, caption=None, directory=None, rich_cols=None,
                     no_squeeze=None):
    """Write `df` to <directory>/<out_name>.pdf, booktabs-styled, on Letter pages.

    aligns      per-column 'left'/'center'/'right' (default: left; numbers centre nicely)
    note        method footnote(s), always shown: str or list of lines
    star_note   clarification of the '*' marker; shown only if a '*' appears
    group_cols  columns whose change marks a group; a faint rule is drawn between groups
    orient      force 'portrait'/'landscape'; default auto (portrait unless illegible)
    detail_cols columns rendered as a full-width, labelled second row under each entry's
                metadata row (for long free text); they leave the grid. PDF only.
    rich_cols   grid columns whose cell values are split by gene symbol
                (bap1figs.richtext.split_runs) and drawn with the symbol
                italicized (e.g. 'BAP1; HGNC:950') instead of as plain text.
                Values must be short enough not to need wrapping.
    no_squeeze  grid columns exempted from column-width squeezing, so a value
                with no good word-break point (a formatted CI like
                '(0.0028-0.0441)') never gets hard-split mid-number — other
                columns give up the width instead. Use for short, unbreakable,
                information-dense columns; overusing it just pushes the
                squeeze onto whatever's left.
    title       bold heading above the table on the first page (e.g. 'Supplementary Table 1. ...')
    caption     descriptive caption paragraph(s) under the title (str or list of lines)
    directory   output directory (default config.TABLES)
    """
    config.set_rcparams()
    directory = directory or config.TABLES
    all_cols = list(df.columns)
    detail_cols = [c for c in (detail_cols or []) if c in all_cols]
    cols = [c for c in all_cols if c not in detail_cols]         # grid columns
    rich_cols = set(rich_cols or []) & set(cols)
    no_squeeze = set(no_squeeze or []) & set(cols)
    aligns = aligns or {}
    cells = {c: [_fmt(v) for v in df[c]] for c in all_cols}
    natural = [max([_text_w(c, fontsize, bold=True)]
                   + [_text_w(v, fontsize) for v in cells[c]]) + CELL_PAD for c in cols]
    floor = [max(_word_floor([c], fontsize, bold=True),
                 _word_floor(cells[c], fontsize)) for c in cols]

    portrait, landscape = LETTER, (LETTER[1], LETTER[0])

    def _lines(wrapw, col, s, fs):
        if wrapw[col] and s:
            return _wrap(s, wrapw[col], fs)
        return [s]

    def _plan(page, fontsize):
        widths, wrapw, fs = _fit(cols, natural, floor, page, fontsize, no_squeeze)
        rows = [[_lines(wrapw, c, cells[c][r], fs) for c in cols] for r in range(len(df))]
        nls = [max(len(w) for w in row) for row in rows]
        line_h = fs * LEADING / 72.0
        avail = (page[1] - 2 * MARGIN) - line_h * 1.9
        # greedy page count
        pages, used = 1, 0.0
        for nl in nls:
            h = nl * line_h
            if used + h > avail and used > 0:
                pages += 1; used = 0.0
            used += h
        return dict(page=page, widths=widths, wrapw=wrapw, fontsize=fs,
                    rows=rows, nls=nls, pages=pages, total_lines=sum(nls))

    if orient == "portrait":
        plan = _plan(portrait, fontsize)
    elif orient == "landscape":
        plan = _plan(landscape, fontsize)
    else:
        pp, lp = _plan(portrait, fontsize), _plan(landscape, fontsize)
        # fewer pages wins; then less wrapping (rules out cramped same-page layouts);
        # portrait wins a true tie (so narrow/tall tables stay upright).
        plan = pp if (pp["pages"], pp["total_lines"]) <= (lp["pages"], lp["total_lines"]) else lp
    page, widths, fontsize = plan["page"], plan["widths"], plan["fontsize"]
    wrapw = plan["wrapw"]
    body = list(zip(plan["rows"], plan["nls"]))
    table_w = sum(widths)

    # Group boundaries: rule drawn *above* a row whose group key differs from the previous.
    breaks = set()
    if group_cols:
        keys = ["".join(str(df.iloc[r][c]) for c in group_cols) for r in range(len(df))]
        breaks = {r for r in range(1, len(df)) if keys[r] != keys[r - 1]}

    line_h = fontsize * LEADING / 72.0
    # Wrap header labels to their column width (bold text is a touch wider).
    hlines = [_wrap(cols[j], wrapw[cols[j]], fontsize, bold=True) for j in range(len(cols))]
    header_h = (max(len(h) for h in hlines) + 0.9) * line_h
    raw_notes = []
    if note:
        raw_notes += [note] if isinstance(note, str) else list(note)
    if star_note and any("*" in cells[c][r] for c in cols for r in range(len(df))):
        raw_notes.append(star_note)
    foot_fs = max(5.5, fontsize - 1.5)
    foot_width = page[0] - 2 * MARGIN
    footnotes = [ln for fn in raw_notes for ln in (_wrap(fn, foot_width, foot_fs) or [""])]
    note_h = (0.6 + len(footnotes)) * line_h if footnotes else 0.0

    # Title + caption block (first page only): bold heading, then a caption
    # paragraph. Rich text (gene symbols + "et al." italic) laid out with
    # _text_w — the same matplotlib glyph metrics used to draw it — so the
    # word spacing here matches the rest of the table instead of drifting
    # (laying out with one font's metrics and drawing with another's is what
    # produced unevenly-gapped title/caption text before this).
    title_fs = fontsize + 2.0
    cap_fs = max(6.0, fontsize - 0.5)
    cap_line_h = cap_fs * 1.35 / 72.0
    content_w_in = page[0] - 2 * MARGIN
    title_lines = _layout_rich(richtext.normalize(title), content_w_in, title_fs, bold=True) if title else []
    raw_cap = ([caption] if isinstance(caption, str) else list(caption)) if caption else []
    caption_lines = [ln for c in raw_cap
                     for ln in _layout_rich(richtext.normalize(c), content_w_in, cap_fs, bold=False)]
    title_block_h = 0.0
    if title_lines:
        title_block_h += len(title_lines) * (title_fs * 1.25 / 72.0) + 0.06
    if caption_lines:
        title_block_h += len(caption_lines) * cap_line_h + 0.10
    if title_block_h:
        title_block_h += 0.12                                   # gap before the table

    avail_h = page[1] - 2 * MARGIN
    x0 = (page[0] - table_w) / 2.0
    xedges = [x0]
    for w in widths:
        xedges.append(xedges[-1] + w)

    def cell_x(j):
        al = aligns.get(cols[j], "left")
        if al == "left":
            return xedges[j] + CELL_PAD / 2, "left"
        if al == "right":
            return xedges[j + 1] - CELL_PAD / 2, "right"
        return (xedges[j] + xedges[j + 1]) / 2, "center"

    def _draw_rich(ax, x, y, text, ha, fs):
        """Draw text as separate runs (bap1figs.richtext.split_runs), italicizing
        gene symbols, at the same (x, ha) a plain cell would use. For a short,
        single-line, non-wrapping value only (e.g. 'BAP1; HGNC:950')."""
        runs = richtext.split_runs(text)
        total = sum(_text_w(s, fs, italic=it) for s, it in runs)
        x0_run = x - total / 2 if ha == "center" else (x - total if ha == "right" else x)
        for s, it in runs:
            ax.text(x0_run, y, s, ha="left", va="center", color=TEXT, fontsize=fs,
                    fontstyle="italic" if it else "normal", zorder=3)
            x0_run += _text_w(s, fs, italic=it)

    # Full-width detail lines (long free-text columns), wrapped to the table width.
    # Each detail line is (text, underline_width_in): the former column header (the
    # label at the start of a field's first line) is underlined.
    detail_indent = 0.16
    detail_width = table_w - detail_indent - CELL_PAD
    detail_body = []
    for r in range(len(df)):
        dl = []
        for dc in detail_cols:
            v = cells[dc][r]
            if v:
                wrapped = _wrap(f"{dc}: {v}", detail_width, fontsize)
                for idx, ln in enumerate(wrapped):
                    dl.append((ln, _text_w(dc, fontsize) if idx == 0 else None))
        detail_body.append(dl)

    # Flatten to visual lines so tall entries flow across pages without clipping.
    # Each visual line: (row_index, kind 'grid'/'detail', line_index, is_row_start).
    vlines = []
    for r, (wrapped, nl) in enumerate(body):
        for k in range(nl):
            vlines.append((r, "grid", k, k == 0))
        for k in range(len(detail_body[r])):
            vlines.append((r, "detail", k, False))

    def draw_header(ax, top):
        ax.plot([x0, x0 + table_w], [top, top], color=RULE, lw=1.2, zorder=3)
        for j in range(len(cols)):
            xt, ha = cell_x(j)
            ax.text(xt, top - header_h / 2, "\n".join(hlines[j]), ha=ha, va="center",
                    color=TEXT, fontsize=fontsize, fontweight="bold", zorder=3, linespacing=1.2)
        ax.plot([x0, x0 + table_w], [top - header_h, top - header_h], color=RULE, lw=0.7, zorder=3)

    def draw_title_block(ax, top):
        """Draw the title + caption at the top of the first page; return the new top."""
        y = top
        xl = MARGIN
        for line in title_lines:
            y -= title_fs * 1.25 / 72.0
            for s, it, xin in line:
                ax.text(xl + xin, y, s, ha="left", va="baseline", color=TEXT,
                        fontsize=title_fs, fontweight="bold",
                        fontstyle="italic" if it else "normal", zorder=3)
        if title_lines:
            y -= 0.06
        for line in caption_lines:
            y -= cap_line_h
            for s, it, xin in line:
                ax.text(xl + xin, y, s, ha="left", va="baseline", color=FOOT,
                        fontsize=cap_fs, fontstyle="italic" if it else "normal", zorder=3)
        if caption_lines:
            y -= 0.10
        return y - 0.12 if title_block_h else y

    cap = max(1, int((avail_h - header_h - note_h) / line_h))            # lines/page (rest)
    cap_first = max(1, int((avail_h - title_block_h - header_h - note_h) / line_h))

    out = directory / f"{out_name}.pdf"
    with PdfPages(out) as pdf:
        i, page_no = 0, 0
        while i < len(vlines):
            this_cap = cap_first if page_no == 0 else cap
            chunk = vlines[i:i + this_cap]
            i += this_cap
            page_no += 1
            last_page = i >= len(vlines)

            fig = plt.figure(figsize=page)
            ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, page[0]); ax.set_ylim(0, page[1])
            ax.axis("off")
            top = page[1] - MARGIN
            if page_no == 1 and title_block_h:
                top = draw_title_block(ax, top)
            draw_header(ax, top)

            y = top - header_h
            for pos, (r, kind, k, is_start) in enumerate(chunk):
                if is_start and r in breaks and pos != 0:
                    ax.plot([x0, x0 + table_w], [y, y], color=GROUP_RULE, lw=0.4, zorder=2)
                if kind == "grid":
                    wrapped = body[r][0]
                    for j in range(len(cols)):
                        if k < len(wrapped[j]):
                            xt, ha = cell_x(j)
                            if cols[j] in rich_cols:
                                _draw_rich(ax, xt, y - line_h / 2, wrapped[j][k], ha, fontsize)
                            else:
                                ax.text(xt, y - line_h / 2, wrapped[j][k], ha=ha, va="center",
                                        color=TEXT, fontsize=fontsize, zorder=3)
                else:
                    text, ul = detail_body[r][k]
                    xd = x0 + detail_indent
                    ax.text(xd, y - line_h / 2, text, ha="left", va="center",
                            color=DETAIL, fontsize=fontsize, zorder=3)
                    if ul:
                        yb = y - line_h / 2 - fontsize * 0.52 / 72.0
                        ax.plot([xd, xd + ul], [yb, yb], color=DETAIL, lw=0.5, zorder=3)
                y -= line_h

            ax.plot([x0, x0 + table_w], [y, y], color=RULE, lw=1.2, zorder=3)        # bottom rule
            if footnotes and last_page:
                for kk, fn in enumerate(footnotes):
                    ax.text(x0, y - line_h * (0.9 + kk * 0.95), fn, ha="left", va="top",
                            fontsize=foot_fs, style="italic", color=FOOT)
            pdf.savefig(fig); plt.close(fig)
    print(f"INFO: wrote {out.relative_to(config.REPO)}")
    return out
