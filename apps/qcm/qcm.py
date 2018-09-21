#!/usr/bin/python
#coding: utf-8

import sys
import matplotlib.pyplot as plt
import math

#Usage: ./qcm.py input.txt input.tiff output.tiff

assert(len(sys.argv) > 3)

# Import note
# f = open(sys.argv[1])
# Rrep = [ line.strip() for line in f ]
# f.close()
Rrep = [ "bc", "a", "bd", "ac", "c", "b", "c", "c", "ad", "b",
         "abd", "abd", "b", "c", "a", "a", "ac", "ac", "b", "a",
         "c", "a", "c", "c", "ab", "ab", "ad", "a", "c", "c" ]


# Import student
f = open(sys.argv[1])
login = f.readline().strip()
Srep = [ line.strip() for line in f ]
f.close()

#print Srep, len(Rrep)
# Get the note
notes = [ 1 if u == v else 0 if v == "" else -.25 for u,v in zip(Rrep, Srep) ]

# Print
note = 25.0 * sum(notes) / len(notes)
note = math.ceil(note * 2) / 2  # arrondi à 0.5
print "%s,%.01f,%s" % (login, note, sys.argv[2])

###############################
##  Write the image         ###
###############################

img = plt.imread(sys.argv[2])
ypixels, xpixels, _ = img.shape

dpi = 200.
xinch = xpixels / dpi
yinch = ypixels / dpi

fig = plt.figure(figsize=(xinch,yinch))
ax = plt.axes([0., 0., 1., 1.], frameon=False, xticks=[],yticks=[])
ax.imshow(img)

x = 450
y = 1750
for k,n in enumerate(notes):
    i = k / 10
    j = k % 10
    plt.text(x + 481 * i, 1750 + j * 100, str(n), color="red")

plt.text(180,1000, login + " : " +  str(note), color="red", fontsize=15)
#plt.show()
plt.savefig(sys.argv[3], dpi=dpi)
