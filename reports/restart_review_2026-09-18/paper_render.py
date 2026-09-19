"""Typeset the report on A4 pages using local, installed plotting libraries.

Individual figure PDFs retain vector graphics. The integrated reading PDF uses
300-dpi figure previews, serif text, numbered sections, real tables and captions.
"""
from __future__ import annotations

import re
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextToPath
from PIL import Image


def plain(text):
    text = re.sub(r'\[([^]]+)\]\(([^)]+)\)', r'\1', text)
    return text.replace('**', '').replace('`', '').replace('\\_', '_')


class Paper:
    width, height = 8.27, 11.69
    left, right, top, bottom = .78, .67, .85, .72

    def __init__(self, pdf):
        self.pdf = pdf
        self.page = 0
        self.fig = None
        self.y = 0
        self.metrics = TextToPath()
        self.new_page()

    @property
    def available_width(self):
        return self.width - self.left - self.right

    def finish_page(self):
        if self.fig is not None:
            self.fig.text(.5, .033, str(self.page), ha='center', fontsize=9, family='DejaVu Serif')
            self.pdf.savefig(self.fig)
            plt.close(self.fig)

    def new_page(self):
        self.finish_page()
        self.page += 1
        self.fig = plt.figure(figsize=(self.width, self.height))
        self.fig.text(self.left/self.width, .965, 'TMC1A: continuum and ice-band modelling',
                      fontsize=8, color='.35', family='DejaVu Serif')
        self.fig.text(1-self.right/self.width, .965, 'Local results review · 18 September 2026',
                      ha='right', fontsize=7.5, color='.35', family='DejaVu Serif')
        self.y = self.height - self.top

    def ensure(self, inches):
        if self.y - inches < self.bottom:
            self.new_page()

    def wrap(self, text, size=10, width=None, weight='normal'):
        width = (width or self.available_width) * 72
        prop = FontProperties(family='DejaVu Serif', size=size, weight=weight)
        lines, line = [], ''
        for word in plain(text).split():
            trial = (line + ' ' + word).strip()
            measured = self.metrics.get_text_width_height_descent(trial, prop, ismath=False)[0]
            if measured > width and line:
                lines.append(line)
                line = word
            else:
                line = trial
        if line:
            lines.append(line)
        return lines or ['']

    def lines(self, lines, size=10, weight='normal', color='black', gap=.1, line_height=None):
        step = line_height or size * 1.38 / 72
        self.ensure(min(len(lines), 3)*step + gap)
        for index, line in enumerate(lines):
            # Keep the last two lines together where possible.
            self.ensure((2 if len(lines)-index == 2 else 1)*step)
            self.fig.text(self.left/self.width, self.y/self.height, line, va='top',
                          fontsize=size, fontweight=weight, color=color, family='DejaVu Serif')
            self.y -= step
        self.y -= gap

    def paragraph(self, text, size=10, weight='normal', gap=.10):
        self.lines(self.wrap(text, size, weight=weight), size, weight, gap=gap)

    def heading(self, text, level):
        size = {1: 18, 2: 12.5, 3: 11}.get(level, 10.5)
        lines = self.wrap(text, size, weight='bold')
        self.ensure(len(lines)*size*1.4/72 + .6)
        self.y -= .10 if level > 1 else 0
        self.lines(lines, size, weight='bold', gap=.13)

    def table(self, rows):
        if not rows:
            return
        n = len(rows[0])
        ratios = [0.24, .47, .29] if n == 3 else [1/n]*n
        widths = [self.available_width*r for r in ratios]
        size, step = 8.1, .153

        def row_height(row):
            wrapped = [self.wrap(cell, size, widths[i]-.10) for i, cell in enumerate(row)]
            return wrapped, max(map(len, wrapped))*step+.14

        def draw(row, header=False):
            wrapped, h = row_height(row)
            if header:
                self.fig.add_artist(plt.Rectangle((self.left/self.width, (self.y-h)/self.height),
                    self.available_width/self.width, h/self.height, transform=self.fig.transFigure,
                    facecolor='#edf0f2', edgecolor='none', zorder=0))
            x = self.left
            for i, lines in enumerate(wrapped):
                for j, line in enumerate(lines):
                    self.fig.text((x+.04)/self.width, (self.y-.065-j*step)/self.height, line,
                                  va='top', fontsize=size, fontweight='bold' if header else 'normal',
                                  family='DejaVu Serif')
                x += widths[i]
            self.y -= h
            self.fig.add_artist(plt.Line2D([self.left/self.width, 1-self.right/self.width],
                [self.y/self.height]*2, transform=self.fig.transFigure, color='.7', lw=.45))

        self.ensure(row_height(rows[0])[1] + row_height(rows[1] if len(rows)>1 else rows[0])[1])
        draw(rows[0], True)
        for row in rows[1:]:
            h = row_height(row)[1]
            if self.y-h < self.bottom:
                self.new_page()
                draw(rows[0], True)
            draw(row)
        self.y -= .16

    def figure(self, path, caption):
        with Image.open(path) as source:
            aspect = source.height/source.width
        image_height = min(self.available_width*aspect, 6.8)
        lines = self.wrap(caption, 9)
        caption_height = len(lines)*9*1.38/72+.16
        self.ensure(image_height+caption_height+.10)
        image_width = image_height/aspect
        x = self.left+(self.available_width-image_width)/2
        ax = self.fig.add_axes([x/self.width, (self.y-image_height)/self.height,
                               image_width/self.width, image_height/self.height])
        ax.imshow(plt.imread(path)); ax.axis('off')
        self.y -= image_height+.08
        self.lines(lines, size=9, gap=.18)


def build_paper_pdf(here: Path):
    lines = (here/'REPORT.md').read_text().splitlines()
    with PdfPages(here/'REPORT.pdf', metadata={
        'Title': 'Radiative-transfer modelling of TMC1A: continuum and ice-band constraints',
        'Subject': 'Local modelling results and limitations; no new simulations',
        'Keywords': 'TMC1A, JWST, MCFOST, radiative transfer, dust, ice'}) as pdf:
        paper = Paper(pdf)
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line:
                i += 1; continue
            if line.startswith('|'):
                rows = []
                while i < len(lines) and lines[i].strip().startswith('|'):
                    cells = [x.strip() for x in lines[i].strip().strip('|').split('|')]
                    if not all(set(c) <= set('-: ') for c in cells):
                        rows.append(cells)
                    i += 1
                paper.table(rows); continue
            image = re.fullmatch(r'!\[([^]]*)\]\(([^)]+)\)', line)
            if image:
                caption = image[1]
                next_index = i+1
                while next_index < len(lines) and not lines[next_index].strip():
                    next_index += 1
                if next_index < len(lines) and re.match(r'(?:\*\*)?(?:Figure|Fig\.|Appendix figure)\s', lines[next_index]):
                    caption = lines[next_index]
                    i = next_index
                paper.figure(here/image[2], plain(caption))
            elif line.startswith('#'):
                level = len(line)-len(line.lstrip('#'))
                paper.heading(line.lstrip('# '), level)
            elif line.startswith('- '):
                paper.paragraph('• '+line[2:])
            else:
                paper.paragraph(line)
            i += 1
        paper.finish_page()
