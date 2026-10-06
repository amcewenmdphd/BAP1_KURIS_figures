"""Inline rich text for the supplement captions and titles: italicize gene
symbols and the Latin abbreviation "et al.".

A caption string is normalized (Küry spelling, "et al." punctuation), split into
whitespace words, each word into (text, italic) sub-runs (word_runs), and
word-wrapped to a width. Two renderers consume this: the PyMuPDF figure
composer (build_supplement) draws with fitz/DejaVu Sans and measures with the
matching layout() here; the matplotlib table renderer (tablepdf) draws with
whatever font matplotlib's rcParams actually resolve to (Arial on this system)
and does its own layout with word_runs() + its own matplotlib-metric text
measurement, so the two never mix a layout computed in one font with glyphs
drawn in another.
"""

import re

import fitz
import matplotlib.font_manager as fm

GENES = {"BAP1"}
_GENE_RE = re.compile(r"(" + "|".join(sorted(GENES, key=len, reverse=True)) + r")")

# fitz.Font per (weight, italic) — used for both measuring and (in the composer)
# drawing. The table renderer measures with these and draws with matplotlib in
# the same family, so positions line up.
_STYLE = {
    ("reg", False): dict(style="normal", weight="normal"),
    ("reg", True):  dict(style="italic", weight="normal"),
    ("bld", False): dict(style="normal", weight="bold"),
    ("bld", True):  dict(style="italic", weight="bold"),
}
_FONT = {k: fitz.Font(fontfile=fm.findfont(fm.FontProperties(family="DejaVu Sans", **v)))
         for k, v in _STYLE.items()}


def font(weight, italic):
    """The fitz.Font for a (weight in {'reg','bld'}, italic bool) combination."""
    return _FONT[(weight, bool(italic))]


def normalize(text):
    """Unify the Küry spelling and give 'et al.' its trailing period."""
    text = re.sub(r"\bKury\b", "Küry", text)
    text = re.sub(r"\bet al\b(?!\.)", "et al.", text)
    return text


def _split_word(w):
    """Split one space-free word into (subtext, italic) runs by gene symbol."""
    return [(p, bool(_GENE_RE.fullmatch(p))) for p in _GENE_RE.split(w) if p]


def split_runs(text):
    """Split any text into (substring, italic) runs around gene-symbol
    occurrences (e.g. 'BAP1; HGNC:950' -> [('BAP1', True), ('; HGNC:950',
    False)]) — for drawing a short, non-wrapping cell value with an inline
    italic gene symbol outside the caption/title word-wrap path."""
    return [(p, bool(_GENE_RE.fullmatch(p))) for p in _GENE_RE.split(text) if p]


def word_runs(text):
    """List of words; each word is a list of (subtext, italic) sub-runs.
    'et al.' (or 'et al') is emitted as two italic words. Pure text analysis,
    no font metrics — callers measure/position runs with their own renderer's
    glyph metrics (layout() below does this with fitz for the PyMuPDF figure
    composer; tablepdf does its own equivalent with matplotlib's metrics, since
    on this system matplotlib actually renders with Arial, not DejaVu Sans —
    laying out with one font's spacing and drawing with another's is what
    produced the unevenly-gapped title/caption text this was fixed for)."""
    toks = text.split()
    out, i = [], 0
    while i < len(toks):
        w = toks[i]
        if w == "et" and i + 1 < len(toks) and toks[i + 1] in ("al.", "al"):
            out.append([("et", True)])
            out.append([("al.", True)])
            i += 2
            continue
        out.append(_split_word(w))
        i += 1
    return out


def layout(text, weight, fontsize, width_pt):
    """Greedy word-wrap `text` to `width_pt` points at `fontsize`, measured with
    fitz/DejaVu Sans metrics — correct for the PyMuPDF figure composer, which
    also draws with fitz/DejaVu Sans. Returns a list of lines; each line is a
    list of (subtext, italic, x_pt) with x_pt the left offset of the sub-run
    within the block, in points."""
    sp = font(weight, False).text_length(" ", fontsize=fontsize)
    lines, cur, x = [], [], 0.0
    for word in word_runs(text):
        ww = sum(font(weight, it).text_length(s, fontsize=fontsize) for s, it in word)
        if cur and x + ww > width_pt:
            lines.append(cur)
            cur, x = [], 0.0
        wx = x
        for s, it in word:
            cur.append((s, it, wx))
            wx += font(weight, it).text_length(s, fontsize=fontsize)
        x = wx + sp
    if cur:
        lines.append(cur)
    return lines
