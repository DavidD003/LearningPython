import numpy as np
import csv
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
import math as m
###########################################################################
####  Author: David De Sa						                        ###
#### Contact: daviddesa03@gmail.com					                    ###
#### Support: https://www.patreon.com/DavidD003				            ###
####    More: https://github.com/DavidD003 				                ###
####          https://www.youtube.com/channel/UCr8HDUNoTaT2QIfLktXk3Jg  ###
###########################################################################

# The purpose of this code is to take in a data set, and output a visualization of the
# data variance.

# Source Data
filepath = r'Filepath\FileName.csv'
#   Data should be a .csv file, with 3 columns containing the data in the x, y, z axes

# Labels
xlab = 'x'
ylab = 'y'
zlab = 'z'

# Print Data points?
printpoints = True

# Toggle data points?
tgl = True
#   If false, points will always be present. If True, points will only appear after the 'flower' has grown

# Print principle Component Axes?
printprinciplecomponents = True

# Print Rectangular axis?
printrectaxes = True

# Number of radial samples?
radsamp = 50
#   The circular slices of the plot will have this many points along the edge

# Slices?
slices = 25
#   Determines how many slices will be formed in each direction along the primary axis
#   (For log spacing, equal number appear before and after 1 on the axis)

# Slice Spread
spread = 0
#   Determines whether the spacing of slices will be linear(0) or logarithmic (1)

# Scale
scl = 15
#   Only used for logarithmic slice spacing. The greater the value, the more condensed the slices will be near the origin

# Draw Speed
ps = 0.01
#   The number of seconds to elapse between slices being drawn

# Oscillating?
osc = True
#   Determines whether or not the 'flower' grows and shrinks cyclically instead of just growing

################################################
Program:

pts = []
with open(filepath) as f:
    while True:
        try:
            reader = csv.reader(f)
            row = next(reader)
            pts.append(row[0:3])
        except StopIteration:
            break
pts = np.asarray(pts).astype(np.float).T
xavg = np.sum(pts[0][:]) / pts.shape[1]
yavg = np.sum(pts[1][:]) / pts.shape[1]
zavg = np.sum(pts[2][:]) / pts.shape[1]
xbar = np.array([xavg, yavg, zavg]).reshape(3, 1)
xstdv = (sum([(x-xavg)**2 for x in pts[0][:]])/pts.shape[1])**(1/2)
ystdv = (sum([(y-yavg)**2 for y in pts[1][:]])/pts.shape[1])**(1/2)
zstdv = (sum([(z-zavg)**2 for z in pts[2][:]])/pts.shape[1])**(1/2)
xz = np.asarray([(x - xavg) / xstdv for x in pts[0][:]])
yz = np.asarray([(y - yavg) / ystdv for y in pts[1][:]])
zz = np.asarray([(z - zavg) / zstdv for z in pts[2][:]])
pts = np.vstack((xz, yz, zz))
covm = np.cov(pts)
evals, evects = np.linalg.eig(covm)
e1, ev1 = evects[:, 0], evals[0]
e2, ev2 = evects[:, 1], evals[1]
e3, ev3 = evects[:, 2], evals[2]
e1, e2, e3 = e1/np.linalg.norm(e1), e2 / \
    np.linalg.norm(e2), e3/np.linalg.norm(e3)
maxval = sorted(list(evals), reverse=True)[0]
eigenindex = list(evals).index(maxval)
emaxval = evects[:, eigenindex]
bss = [emaxval]
for i in range(3):
    if i != eigenindex:
        bss.append(evects[:, i])
evs = np.hstack(tuple(bss)).reshape(3, 3, order='F')
invevs = np.linalg.inv(evs)
eigens = (e1, e2, e3)
x = []
y = []
z = []
spr = []
tstpt = []
i = -1
for e in eigens:
    i = i+1
    if i == eigenindex:
        prcp = e
    else:
        spr.append(e)
for i in np.linspace(0, 361, radsamp):
    mypt = list(np.sin((i * m.pi) / 180) *
                spr[0] + np.cos((i * m.pi) / 180) * spr[1])
    mypt.extend([0, 0, 0, 0])
    tstpt.append(mypt)
tstpt1 = []
if spread == 0:
    for i in np.linspace(0.01, maxval, slices):
        for pt in tstpt:
            mypt = list(pt[0:3] + i * prcp)
            mypt.extend([0, 0, 0, i])
            tstpt1.append(mypt)
            mypt = list(pt[0:3] - i * prcp)
            mypt.extend([0, 0, 0, -i])
            tstpt1.append(mypt)
else:
    ls = [i for i in reversed(
        list(np.logspace((1), (1 / scl), num=slices, base=maxval * 1000) / 1000))]
    for i in ls:
        for pt in tstpt:
            mypt = list(pt[0:3] + i * prcp)
            mypt.extend([0, 0, 0, i])
            tstpt1.append(mypt)
            mypt = list(pt[0:3] - i * prcp)
            mypt.extend([0, 0, 0, -i])
            tstpt1.append(mypt)
tstpt.extend(tstpt1)
allpts = []
allpts.append([0, 0, 0, 0, 0, 0, 0])
for pt in tstpt:
    curr = pt[0:3]
    curr = curr / np.linalg.norm(curr)
    new = (covm @ curr).tolist()
    new1 = list(new)
    new1.extend(pt[3:])
    allpts.append(new1)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel(xlab)
ax.set_ylabel(ylab)
ax.set_zlabel(zlab)
if printrectaxes == True:
    axpts = [[[-1, 1], [0, 0], [0, 0]],
             [[0, 0], [-1, 1], [0, 0]], [[0, 0], [0, 0], [-1, 1]]]
    for axpt in axpts:
        ax.plot3D(axpt[0], axpt[1], axpt[2], color='black')
if printpoints == True and tgl == False:
    dta = [ax.scatter3D(pts[0][:], pts[1][:], pts[2][:], color='black')]
if printprinciplecomponents == True:
    ax.plot3D([-ev1*e1[0], ev1*e1[0]], [-ev1*e1[1], ev1*e1[1]],
              [-ev1*e1[2], ev1*e1[2]], color='blue')
    ax.plot3D([-ev2*e2[0], ev2*e2[0]], [-ev2*e2[1], ev2*e2[1]],
              [-ev2*e2[2], ev2*e2[2]], color='blue')
    ax.plot3D([-ev3*e3[0], ev3*e3[0]], [-ev3*e3[1], ev3*e3[1]],
              [-ev3*e3[2], ev3*e3[2]], color='blue')
plt.gca().set_xlim(left=-4, right=4)
plt.gca().set_ylim(bottom=-4, top=4)
plt.gca().set_zlim(zmin=-4, zmax=4)
allx = [x[0] for x in allpts if x[len(x) - 1] == 0]
ally = [x[1] for x in allpts if x[len(x) - 1] == 0]
allz = [x[2] for x in allpts if x[len(x) - 1] == 0]
ax.plot_trisurf(allx, ally, allz,
                cmap='viridis', edgecolor='none')
plt.draw()
plt.pause(ps)
sct = []
if spread == 0:
    while True:
        for i in np.linspace(0.01, maxval, slices):
            boole = [i, -i]
            for j in range(2):
                allx = [x[0] for x in allpts if x[len(x) - 1] == boole[j]]
                allx.extend([x[3]
                             for x in allpts if x[len(x) - 1] == boole[j]])
                ally = [x[1] for x in allpts if x[len(x)-1] == boole[j]]
                ally.extend([x[4] for x in allpts if x[len(x)-1] == boole[j]])
                allz = [x[2] for x in allpts if x[len(x)-1] == boole[j]]
                allz.extend([x[5] for x in allpts if x[len(x)-1] == boole[j]])
                sct.append(ax.plot_trisurf(allx, ally, allz,
                                           cmap='viridis', edgecolor='none'))
            plt.draw()
            plt.pause(ps)
        if printpoints == True and tgl == True:
            dta = [ax.scatter3D(pts[0][:], pts[1][:],
                                pts[2][:], color='black')]
            plt.draw()
        if osc != True:
            break
        else:
            plt.pause(3)
            if printpoints == True and tgl == True:
                dta[0].remove()
                plt.draw()
                plt.pause(0.5)
            j = 0
            while len(sct) > 0:
                j = j+1
                sct[len(sct)-1].remove()
                sct.remove(sct[len(sct)-1])
                if j % 2 == 0:
                    plt.draw()
                    plt.pause(ps)
else:
    while True:
        for i in ls:
            boole = [i, -i]
            for j in range(2):
                allx = [x[0] for x in allpts if x[len(x) - 1] == boole[j]]
                allx.extend([x[3]
                             for x in allpts if x[len(x) - 1] == boole[j]])
                ally = [x[1] for x in allpts if x[len(x)-1] == boole[j]]
                ally.extend([x[4] for x in allpts if x[len(x)-1] == boole[j]])
                allz = [x[2] for x in allpts if x[len(x)-1] == boole[j]]
                allz.extend([x[5] for x in allpts if x[len(x)-1] == boole[j]])
                sct.append(ax.plot_trisurf(allx, ally, allz,
                                           cmap='viridis', edgecolor='none'))
            plt.draw()
            plt.pause(ps)
        if printpoints == True and tgl == True:
            dta = [ax.scatter3D(pts[0][:], pts[1][:],
                                pts[2][:], color='black')]
            plt.draw()
        if osc != True:
            break
        else:
            plt.pause(3)
            if printpoints == True and tgl == True:
                dta[0].remove()
                plt.draw()
                plt.pause(0.5)
            j = 0
            while len(sct) > 0:
                j = j+1
                sct[len(sct)-1].remove()
                sct.remove(sct[len(sct)-1])
                if j % 2 == 0:
                    plt.draw()
                    plt.pause(ps)
plt.show()
