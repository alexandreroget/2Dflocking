#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 16:27:31 2019

@author: aroget
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import *

def readFile(measures,n):
    x1 = np.zeros(n)
    x2 = np.zeros(n)
    v1 = np.zeros(n)
    v2 = np.zeros(n)
    with open(measures, "r") as fichier:
        for i in range(n):
            line = fichier.readline()
            splited = line.split()
            x1[i] = float(splited[0])
            x2[i] = float(splited[2])
            v1[i] = float(splited[4])
            v2[i] = float(splited[6])

    return x1, x2, v1, v2

n_files = 100
n_particles = 1000
a = 5  # Size of the box

axes = gca()
line, = axes.plot([],[],'rx')
xlim(0.0, 2 * a)
ylim(0.0, 2 * a)

for i in range(n_files):
    filename = "data/" + str(i) + ".txt"
    x1,x2,v1,v2 = readFile(filename,n_particles)

    # Plot the data
    line.set_xdata(x1)
    line.set_ydata(x2)
    draw()

    # and display the mean velocity
    title(
        r"V1 = {0:.5f}".format(sum(v1) / n_particles) + r" V2 = {0:.5f}".format(sum(v2) / n_particles)
    )

    pause(0.05)

show()
