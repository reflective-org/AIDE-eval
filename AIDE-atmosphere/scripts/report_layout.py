"""
Minimal text-page layout engine for building a multi-page PDF with matplotlib.

No LaTeX, no reportlab - matplotlib's PdfPages is the only dependency available
on this machine, so text pages are laid out as figures with a y-cursor.

All pages are A4 landscape (11.69 x 8.27 in) so text and figures interleave
without the page size jumping around.
"""

import textwrap
import matplotlib.pyplot as plt

PAGE_W, PAGE_H = 11.69, 8.27

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#3d3c3a"
MUTED = "#6f6e6a"
RULE = "#d9d8d3"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
HILITE = "#eef4fc"
WARN = "#fdf0ea"

L, R = 0.055, 0.955          # left / right margins in figure fraction
TOP, BOT = 0.925, 0.055


class Page:
    """One text page. Call the block methods in order; they advance a cursor."""

    def __init__(self, pdf, title=None, kicker=None, page_no=None):
        self.pdf = pdf
        self.fig = plt.figure(figsize=(PAGE_W, PAGE_H), facecolor=SURFACE)
        self.y = TOP
        self.page_no = page_no
        if kicker:
            self.fig.text(L, self.y + 0.035, " ".join(kicker.upper()),
                          fontsize=6.8, color=MUTED, ha="left")
        if title:
            self.fig.text(L, self.y, title, fontsize=16, color=INK, ha="left",
                          va="top", weight="bold")
            self.y -= 0.075
            self.fig.add_artist(
                plt.Line2D([L, R], [self.y + 0.018, self.y + 0.018],
                           transform=self.fig.transFigure, color=RULE, lw=1.0))
            self.y -= 0.012

    # ------------------------------------------------------------ blocks
    def heading(self, text, gap=0.030):
        self.y -= gap
        self.fig.text(L, self.y, text, fontsize=11.5, color=INK, ha="left",
                      va="top", weight="bold")
        self.y -= 0.033

    def sub(self, text, gap=0.018):
        self.y -= gap
        self.fig.text(L, self.y, text, fontsize=9.5, color=INK2, ha="left",
                      va="top", weight="bold")
        self.y -= 0.028

    def para(self, text, size=9.0, color=INK2, width=132, gap=0.012, indent=0.0):
        self.y -= gap
        for line in textwrap.wrap(text, width=width) or [""]:
            self.fig.text(L + indent, self.y, line, fontsize=size, color=color,
                          ha="left", va="top")
            self.y -= 0.0225
        return self

    def bullets(self, items, size=9.0, width=126, gap=0.010):
        self.y -= gap
        for it in items:
            lines = textwrap.wrap(it, width=width)
            for j, line in enumerate(lines):
                if j == 0:
                    self.fig.text(L + 0.008, self.y, "•", fontsize=size,
                                  color=BLUE, ha="left", va="top")
                self.fig.text(L + 0.022, self.y, line, fontsize=size, color=INK2,
                              ha="left", va="top")
                self.y -= 0.0225
            self.y -= 0.004

    def eq(self, text, gap=0.014, color=INK):
        """A display equation / formula, monospace, indented."""
        self.y -= gap
        self.fig.text(L + 0.030, self.y, text, fontsize=9.5, color=color,
                      ha="left", va="top", family="DejaVu Sans Mono")
        self.y -= 0.030

    def code(self, lines, gap=0.010, size=8.2):
        self.y -= gap
        for line in lines:
            self.fig.text(L + 0.022, self.y, line, fontsize=size, color=INK2,
                          ha="left", va="top", family="DejaVu Sans Mono")
            self.y -= 0.0195
        self.y -= 0.004

    def callout(self, text, color=BLUE, bg=HILITE, width=124, gap=0.016,
                size=9.0):
        """A boxed remark that should not be missed."""
        self.y -= gap
        lines = textwrap.wrap(text, width=width)
        h = 0.0225 * len(lines) + 0.020
        self.fig.add_artist(
            plt.Rectangle((L - 0.010, self.y - h + 0.016), (R - L) + 0.020, h,
                          transform=self.fig.transFigure, facecolor=bg,
                          edgecolor="none", zorder=0))
        self.fig.add_artist(
            plt.Rectangle((L - 0.010, self.y - h + 0.016), 0.0035, h,
                          transform=self.fig.transFigure, facecolor=color,
                          edgecolor="none", zorder=1))
        self.y -= 0.004
        for line in lines:
            self.fig.text(L + 0.004, self.y, line, fontsize=size, color=INK,
                          ha="left", va="top", zorder=2)
            self.y -= 0.0225
        self.y -= 0.008

    def table(self, headers, rows, widths, size=8.2, gap=0.018, colors=None,
              header_size=8.0):
        """widths are fractions of the text column, summing to <= 1."""
        self.y -= gap
        span = R - L
        xs = [L]
        for w in widths[:-1]:
            xs.append(xs[-1] + w * span)
        for x, htxt in zip(xs, headers):
            self.fig.text(x, self.y, htxt, fontsize=header_size, color=MUTED,
                          ha="left", va="top", weight="bold")
        self.y -= 0.024
        self.fig.add_artist(
            plt.Line2D([L, R], [self.y + 0.010, self.y + 0.010],
                       transform=self.fig.transFigure, color=RULE, lw=0.8))
        self.y -= 0.006
        for i, row in enumerate(rows):
            wrapped, nlines = [], 1
            for j, cell in enumerate(row):
                cw = max(int(widths[j] * 168) - 2, 8)
                wr = textwrap.wrap(str(cell), width=cw) or [""]
                wrapped.append(wr)
                nlines = max(nlines, len(wr))
            rc = (colors[i] if colors and colors[i] else INK2)
            for k in range(nlines):
                for j, wr in enumerate(wrapped):
                    if k < len(wr):
                        self.fig.text(xs[j], self.y, wr[k], fontsize=size,
                                      color=rc, ha="left", va="top")
                self.y -= 0.0205
            self.y -= 0.005
        return self

    def spacer(self, h=0.02):
        self.y -= h

    def footer(self, text=""):
        self.fig.text(L, 0.030, text, fontsize=7.2, color=MUTED, ha="left")
        if self.page_no is not None:
            self.fig.text(R, 0.030, str(self.page_no), fontsize=7.2,
                          color=MUTED, ha="right")

    def close(self, footer_text=""):
        self.footer(footer_text)
        self.pdf.savefig(self.fig, facecolor=SURFACE)
        plt.close(self.fig)


def figure_page(pdf, fig, page_no=None, footer_text=""):
    """Save an already-built figure as one page, stamping the page number."""
    if page_no is not None:
        fig.text(R, 0.020, str(page_no), fontsize=7.2, color=MUTED, ha="right")
    if footer_text:
        fig.text(L, 0.020, footer_text, fontsize=7.2, color=MUTED, ha="left")
    pdf.savefig(fig, facecolor=SURFACE)
    plt.close(fig)
