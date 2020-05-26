###########################################################################
####  Author: David De Sa						                        ###
#### Contact: daviddesa03@gmail.com					                    ###
#### Support: https://www.patreon.com/DavidD003				            ###
####    More: https://github.com/DavidD003 				                ###
####          https://www.youtube.com/channel/UCr8HDUNoTaT2QIfLktXk3Jg  ###
###########################################################################

import csv
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt

# This code takes in 3 Dimensional data, and yields an animation of the line created by the past few historical data points
# illustrating the pathway that data took through the data space
filepath = r'FilePath\filename.csv'
# Data should be .csv, with the data  in the first 3 columns, and figure labels in the 4th column

# Decision Variables:

#   Axis labels:
FromHeader = False
#       Set to True to pull from header, if false, manually enter below
xlab = 'x'
ylab = 'y'
zlab = 'z'

#       The data columns are read in x-y-z format from left to right

#   Caption?
capt = True
#       If true, data from the 4th column is taken as text to be printed out on the screen

#   How many data points long should the snake be?
length = 3

#   How many tail (historical) points to draw?
tail = 10000000

#   Snake growth granularity (>=1)
grain = 20
#       Enter a number thats one or larger.
#       This will determine how coarsely the snake grows
#       (1)  yields the most coarse animation

#   Frames per Second
fps = 100
#       Animation speed: Excessively high values will not show an image.

#   Lower Dimensional projections?
xy = True
xz = True
yz = True
x = True
y = True
z = True
#       Determines whether or not data will be projected onto the planes formed by axis intercepts, and axes

#   Create surface on completion?
surf = False
#       Recommend low surface density for large datasets, so enter a high number below:
#   Surface Density
dens = 1
#       One of every n points will be taken to construct a surface
#       Using '1' will use every point

tails = []
sct = []
cpts = []
hist = [[], [], []]
dta = -1
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
with open(filepath) as f:
    reader = csv.reader(f)
    if FromHeader == True:
        row = next(reader)
        xlab = row[0]
        ylab = row[1]
        zlab = row[2]
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_zlabel(zlab)
    row = next(reader)
    x1 = [float(row[0])]
    y1 = [float(row[1])]
    z1 = [float(row[2])]
    while True:
        try:
            dta = dta+1
            row = next(reader)
            x1.append(float(row[0]))
            y1.append(float(row[1]))
            z1.append(float(row[2]))
            if capt == True:
                capttxt = row[3]
            lead = [float(row[0]), float(row[1]), float(row[2])]
            if dta % dens == 0:
                hist[0].append(lead[0])
                hist[1].append(lead[1])
                hist[2].append(lead[2])
            x2 = x1.copy()
            y2 = y1.copy()
            z2 = z1.copy()
            tails.append([x2[0], y2[0], z2[0]])
            if len(tails) > tail:
                tails.remove(tails[0])
            if len(x1) > length:
                state = 2
                x1.remove(x1[0])
                y1.remove(y1[0])
                z1.remove(z1[0])
            else:
                state = 1
        except StopIteration:
            if len(x1) == 1:
                break
            state = 3
            x2 = x1.copy()
            y2 = y1.copy()
            z2 = z1.copy()
            x1.remove(x1[0])
            y1.remove(y1[0])
            z1.remove(z1[0])
            tails.append([x2[0], y2[0], z2[0]])
            if len(tails) > tail:
                tails.remove(tails[0])
        grow = [(x2[len(x2) - 1] - x2[len(x2) - 2])/grain, (y2[len(y2) - 1] -
                                                            y2[len(y2) - 2])/grain, (z2[len(z2) - 1] - z2[len(z2) - 2])/grain]
        shrink = [(x2[1]-x2[0])/grain, (y2[1]-y2[0]) /
                  grain, (z2[1]-z2[0])/grain]
        for i in range(grain):
            if state == 1:
                x2[len(x2)-1] = x2[len(x2)-2] + i*grow[0]
                y2[len(y2)-1] = y2[len(y2)-2] + i*grow[1]
                z2[len(z2)-1] = z2[len(z2)-2] + i*grow[2]
            elif state == 2:
                x2[len(x2)-1] = x2[len(x2)-2] + i*grow[0]
                y2[len(y2)-1] = y2[len(y2)-2] + i*grow[1]
                z2[len(z2)-1] = z2[len(z2)-2] + i*grow[2]
                x2[0] = x2[0] + shrink[0]
                y2[0] = y2[0] + shrink[1]
                z2[0] = z2[0] + shrink[2]
            else:
                x2[0] = x2[0] + shrink[0]
                y2[0] = y2[0] + shrink[1]
                z2[0] = z2[0] + shrink[2]
            ax.lines = []
            while len(sct) > 0:
                sct[0].remove()
                sct.remove(sct[0])
            sx = [i[0] for i in tails]
            sx.append(lead[0])
            sy = [i[1] for i in tails]
            sy.append(lead[1])
            sz = [i[2] for i in tails]
            sz.append(lead[2])
            ax.plot3D(x2, y2, z2, color='black')
            sct.append(ax.scatter3D(sx, sy, sz, color='black'))
            zeros = [0] * len(x2)
            xmin = [plt.gca().get_xlim()[0]]*len(x2)
            ymin = [plt.gca().get_ylim()[0]]*len(x2)
            zmin = [plt.gca().get_zlim()[0]] * len(x2)
            xmax = [plt.gca().get_xlim()[1]] * len(x2)
            ymax = [plt.gca().get_ylim()[1]] * len(x2)
            zmax = [plt.gca().get_zlim()[1]] * len(x2)
            if capt == True:
                while len(cpts) > 0:
                    cpts[0].remove()
                    cpts.remove(cpts[0])
                cpts.append(plt.figtext(
                    0.1, 0.9, capttxt, fontsize='xx-large'))
            if xy == True:
                ax.plot3D(x2, y2, zmin, color='red')
            if yz == True:
                ax.plot3D(xmin, y2, z2, color='blue')
            if xz == True:
                ax.plot3D(x2, ymax, z2, color='green')
            if x == True:
                sct.append(ax.scatter3D(
                    x2[len(x2)-1], ymin, zmin, color='red'))
            if y == True:
                sct.append(ax.scatter3D(
                    xmax, y2[len(y2)-1], zmin, color='blue'))
            if z == True:
                sct.append(ax.scatter3D(
                    xmax, ymax, z2[len(z2)-1], color='green'))
            plt.draw()
            plt.pause(fps ** -1)
if surf == True:
    ax.plot_trisurf(hist[0], hist[1], hist[2],
                    cmap='viridis', edgecolor='none')
plt.show()
