"""
Tahmid F. Mehdi
Moses Lab, University of Toronto
Motif logo drawer
Feb 23, 2018

Copyright 2019 Tahmid Mehdi
This file is part of svime.

svime is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

svime is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with svime.  If not, see <http://www.gnu.org/licenses/>.
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import gact
import numpy as np
import pandas as pd
import sys


def motif_logo(pwdFile, n, img):
    """
    Generates motif logo images

    :param string pwdFile: name of variationalPWD file
    :param int n: number of occurrences of the motif
    :param img: output image name
    """
    # load concentration parameters for the motif
    psm = pd.read_table(pwdFile, skiprows=[0], header=None, sep="\t")
    psm.drop(psm.columns[len(psm.columns)-1], axis=1, inplace=True)  # drop last column
    window = len(psm.columns)  # window size
    # normalize psm into expected PWM
    windowSums = psm.sum(axis=0)
    pwm = psm/windowSums
    all_scores = []
    for pos in range(window):
        f = pwm[pwm.columns[pos]]  # load weights for corresponding position
        # calculate information content (IC) for each base in the position
        entropy = -np.sum(f*np.log2(f))
        correction = 3/(2*n*np.log(2))
        IC = 2 - (entropy + correction)
        height = IC*f
        base = [('A', height[0]), ('C', height[1]), ('G', height[2]), ('T', height[3])]
        all_scores.append(base)
    # draw logos
    fig, ax = plt.subplots(figsize=(10, 3))

    x = 1
    maxi = 0
    for scores in all_scores:
        y = 0
        sortedScores = sorted(scores, key=lambda item: item[1])
        for base, score in sortedScores:
            gact.letterAt(base, x, y, score, ax)
            y += score
        x += 1
        maxi = max(maxi, y)

    plt.xticks(range(1,x))
    plt.xlim((0, x))
    plt.ylim((0, maxi))
    plt.ylabel('bits', fontsize=20)

    # Tufte-like axes
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_position('zero')
    ax.spines['bottom'].set_color('white')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    plt.tight_layout()
    plt.savefig(img, dpi=300)


outDir = sys.argv[1]
# count number of motifs discovered
clusterAssignments = pd.read_csv(outDir+'/results/clusters.csv', header=None)
clusters = np.unique(clusterAssignments[1])
for c in clusters:
    n = len(np.where(clusterAssignments[1] == c)[0])  # count number of occurrences of motif c
    motif_logo(outDir+'/results/variationalPWD_motif'+str(c)+".txt", n, outDir+'/results/cluster'+str(c)+'.png')
