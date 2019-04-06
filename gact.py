"""
Tahmid F. Mehdi
Moses Lab, University of Toronto
Motif logo Helper
Feb 23, 2018

This script is from ImportanceOfBeingErnest's answer to:
(https://stackoverflow.com/questions/42615527/sequence-logos-in-matplotlib-aligning-xticks)
"""

import matplotlib as mpl
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties

fp = FontProperties(family="Arial", weight="bold")
globscale = 1.33
LETTERS = { "T" : TextPath((-0.305, 0), "T", size=1, prop=fp),
            "G" : TextPath((-0.384, 0), "G", size=1, prop=fp),
            "A" : TextPath((-0.35, 0), "A", size=1, prop=fp),
            "C" : TextPath((-0.356, 0), "C", size=1, prop=fp) }
COLOR_SCHEME = {'G': 'orange',
                'A': 'green',
                'C': 'blue',
                'T': 'red'}


def letterAt(letter, x, y, yscale=1, ax=None):
    text = LETTERS[letter]

    t = mpl.transforms.Affine2D().scale(1*globscale, yscale*globscale) + \
        mpl.transforms.Affine2D().translate(x,y) + ax.transData
    p = PathPatch(text, lw=0, fc=COLOR_SCHEME[letter],  transform=t)
    if ax != None:
        ax.add_artist(p)
    return p
