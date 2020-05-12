# -*- coding: utf-8 -*-
"""Compression Gradient.ipynb
# Partial Image Compressor
The code in this notebook takes in an image, and then applies the Single value Decomposition method to create progressively better approximations for the image. 
The approximations are overlayed to create the appearance of a gradient of deconstruction of an image, centered around a specified area which stays at a high resolution.
# How To Use The Program
In the code block below, modify the necessary paramaters. The required paramters are listed here, Within the code cell, all paramteres have comments explaining how to use them and their purpose.
Once the parameters are set, run the code cell. The image output will be saved to the file pane on the left as 'Gradient Output' which can then be downloaded.
**Filename** : 
Enter a file name, including the filetype extension
**Focus Area** : Enter the coordinates for the top left and bottom right corner of the area of the image that should stay in focus. The coordinates should identify pixels in (row,column) form.
**Rank** : Select the range of image compression ranks that will be used to re create the image. From lower to higher values, compresssion decreases and fidelity increases with exponentially diminishing returns.
**Gradient Type** : Select the type of gradient to apply to the image around the specified Focus Area
Examples of different Gradient Types are found at the bottom of the workbook.
"""
import os
import numpy as np
from matplotlib.image import imsave
from matplotlib.image import imread
import matplotlib.pyplot as plt
# Upload the file you want to use and then enter its name ("use quotes") here:
# Files can be uploaded by dragging them from your file viewer onto the file pane at the left
filename = r'Secret.jpg'

#########################################################
## All parameters entered below must be integer values ##
#########################################################
# Identify the region that you want to remain at high resolution. It is defined by the top left and bottom right corners.
# The format is (row, column). Microsoft Paint can be used to identify the pixel coordinates of a point.
FocusArea = [((95, 368), (192, 433)), ((13, 21), (123, 222))]
# Monkey Island[((95, 368), (192, 433)), ((13, 21), (123, 222))]
# Grad Photo   ((715, 801), (1432, 1257))
# Kawhi ((385, 133), (496, 249)), ((192, 486), (256, 532))
#             r1,c1      r2,c2       #Sepearate distinct focus areas by a comma within the square brackets. Gradient Type must be equal to 4 to manage multiple areas

# Specify the rank of compression. This number corresponds to the number of layers that, together, will construct the gradient
# Larger numbers take longer times to compute. At Rank 80, images usually approximate the original.
MinRank, MaxRank = 1, 9   # Recommended values: 1,15

# Specify the step-value for incrementing rank from one gradient band to the next
RankStep = -1  # Recommended value: -1
#   There is a greater change in image fidelity among lower Ranks. As rank increases, there are diminishing returns.
#   Leave RankStep as -1 to use a non-linear step function that strives to make the rate of compression across the gradient linear.


# Indicate the gradient type as a number
GradientType = 4
# 1 : Square. The gradient will increase compression as it moves away from the focus area
# 2 : Horizontal. The gradient will increase compression as it moves away from the focus area in the vertical direction, creating horizontal bands
# 3 : Vertical. The gradient will increase compression as it moves away from the focus area in the horizontal direction, creating vertical bands
# 4 : Multi-Focus Blending. The gradient will increase compression as it moves away from the focus areas. Gradients from different focus areas will be blended according to their distance from their focus area
#                           To specify multiple focus areas, duplicate the original format using a single comma to separate the areas
# 5 : Cross. The Focus Region entered will be extended upward and across the image, creating a cross shape of focus, with 4 quadrants of vertical gradient decomposition

# Specify how to compute gradient band width as a function of Rank.
BandWidth = 1  # Recommended Value: 1
#   There is a greater change in image fidelity among lower Ranks. As rank increases, there are diminishing returns.
#   The rank is raised to the power of band width to compute a band weight.
#   -Entering 0 yields equal bands
#   -A positive number has band width correspond to rank
#   -A negative number has band width negatively correlate with rank

# Specify Sample Spacing for Multi-Focus Blending. The smaller the spacing, the smoother the gradient blend, at the cost of speed
samplespacing = 1  # Recommended value 1


def rSVD(X, r, q, p):
    ny = X.shape[1]
    P = np.random.randn(ny, r+p)
    Z = X @ P
    for k in range(q):
        Z = X @ (X.T @ Z)
    Q, R = np.linalg.qr(Z, mode='reduced')
    Y = Q.T @ X
    UY, S, VT = np.linalg.svd(Y, full_matrices=0)
    U = Q @ UY
    return U, S, VT


A = imread(filename)
X = np.mean(A, axis=2)
plt.set_cmap('gray')
hires = FocusArea
ranklim = MaxRank
if RankStep == -1:
    ranks = [i**2 for i in range(MinRank, MaxRank)]
    rkcopy = ranks.copy()
    for ch in rkcopy:
        if rkcopy.index(ch) != len(rkcopy)-1:
            if int((ch + ranks[ranks.index(ch) + 1]) / 2) < 400:
                ranks.insert(ranks.index(ch) + 1,
                             int((ch + ranks[ranks.index(ch) + 1]) / 2))
            else:
                ranks = ranks[:ranks.index(ch)]
                break
    ranks.sort()
elif RankStep > 0:
    ranks = [i for i in range(MinRank, MaxRank, RankStep)]
else:
    print('Only non-negative integer rank step values are permitted')
    import sys
    sys.exit()
modcolln = []
for i in ranks:
    U, S, VT = np.linalg.svd(X, full_matrices=0)
    r = i
    XSVD = U[:, :(r+1)] @ np.diag(S[:(r+1)]
                                  ) @ VT[:(r+1), :]
    modcolln.append(XSVD)


def evalfun(rankno):
    return rankno**BandWidth


evals = [evalfun(i) for i in ranks]
wts = [i / sum(evals) for i in evals]
imghgt, imgwdth = len(X), len(X[0])
if GradientType == 4:
    bounds = [[] for i in list(ranks)]
for fcs in hires:
    i = -1
    lastbound = []
    FirstPass = True
    for wgt in wts:
        i = i + 1
        if GradientType == 1:
            if FirstPass == True:
                FirstPass = False
                bounds = list(ranks)
                lastbound = [0, 0, 100000000000, 10000000000]
            lastbound = [fcs[0][0]*wgt+lastbound[0], fcs[0][1]*wgt+lastbound[1],
                         min(imghgt, lastbound[2])-((imghgt-fcs[1][0])*wgt), min(imgwdth, lastbound[3])-((imgwdth-fcs[1][1])*wgt)]
            bounds[i] = (int(lastbound[0]), int(lastbound[1]), int(lastbound[2]),
                         int(lastbound[3]))
        elif GradientType == 2:
            if FirstPass == True:
                FirstPass = False
                bounds = list(ranks)
                lastbound = [0, 0, 1000000000, 0]
            lastbound = [fcs[0][0]*wgt+lastbound[0], 0,
                         min(imghgt, lastbound[2])-((imghgt-fcs[1][0])*wgt), imgwdth-1]
            bounds[i] = (int(lastbound[0]), int(lastbound[1]), int(lastbound[2]),
                         int(lastbound[3]))
        elif GradientType == 3 or GradientType == 5:
            if FirstPass == True:
                FirstPass = False
                bounds = list(ranks)
                lastbound = [0, 0, 0, 1000000000000]
            lastbound = [0, fcs[0][1]*wgt+lastbound[1],
                         imghgt-1, min(imgwdth, lastbound[3])-((imgwdth-fcs[1][1])*wgt)]
            bounds[i] = (int(lastbound[0]), int(lastbound[1]), int(lastbound[2]),
                         int(lastbound[3]))
        elif GradientType == 4:
            if FirstPass == True:
                FirstPass = False
                lastbound = [0, 0, 1000000000000, 1000000000000]
            lastbound = [fcs[0][0]*wgt+lastbound[0], fcs[0][1]*wgt+lastbound[1],
                         min(imghgt, lastbound[2])-((imghgt-fcs[1][0])*wgt), min(imgwdth, lastbound[3])-((imgwdth-fcs[1][1])*wgt)]
            bounds[i].append((int(lastbound[0]), int(lastbound[1]), int(lastbound[2]),
                              int(lastbound[3]), wgt))
if GradientType == 4:
    widthsamps = range(0, imgwdth, samplespacing)
    hghtsamps = range(0, imghgt, samplespacing)
    samps = []
    for w in widthsamps:
        for h in hghtsamps:
            samps.append((h, w))
    for fa in hires:
        samps[:] = [smp for smp in samps if True != (
            smp[0] > fa[0][0] and smp[0] < fa[1][0] and smp[1] > fa[0][1] and smp[1] < fa[1][1])]
    printtiles = []
    while len(samps) > 0:
        point = samps[0]
        TrueNearestFound = [0, 0, 1000000000, 1000000000]
        TrueNearestFoundDelta = [1000000000,
                                 1000000000, -1000000000, -1000000000]
        printweights = []
        printranks = []
        for fa in range(0, len(hires)):
            NearestFound = [10000000000000, 0, wts[0]]
            for rindex in range(0, len(ranks)):
                if True != (point[0] > hires[fa][1][0] and point[1] > hires[fa][0][1] and point[1] < hires[fa][1][1]):
                    if (point[0]-bounds[rindex][fa][0] < point[0]-NearestFound[0]) and (point[0]-bounds[rindex][fa][0] > 0):
                        NearestFound = [point[0]-bounds[rindex][fa]
                                        [0], rindex+1, bounds[rindex+1][fa][4]]
                else:
                    if (bounds[rindex][fa][0]-point[0] < NearestFound[0]-point[0]) and (bounds[rindex][fa][0]-point[0] > 0):
                        NearestFound = [bounds[rindex][fa][0] -
                                        point[0], rindex+1, bounds[rindex+1][fa][4]]
                if (point[0]-bounds[rindex][fa][0] > 0) and (point[0]-bounds[rindex][fa][0] < TrueNearestFoundDelta[0]):
                    TrueNearestFound[0] = bounds[rindex][fa][0]
                    TrueNearestFoundDelta[0] = point[0]-bounds[rindex][fa][0]
                if (point[0]-bounds[rindex][fa][2] > 0) and (point[0]-bounds[rindex][fa][2] < TrueNearestFoundDelta[0]):
                    TrueNearestFound[0] = bounds[rindex][fa][2]
                    TrueNearestFoundDelta[0] = point[0]-bounds[rindex][fa][2]
                if (point[0]-bounds[rindex][fa][0] < 0) and (point[0]-bounds[rindex][fa][0] > TrueNearestFoundDelta[2]):
                    TrueNearestFound[2] = bounds[rindex][fa][0]
                    TrueNearestFoundDelta[2] = point[0]-bounds[rindex][fa][0]
                if (point[0]-bounds[rindex][fa][2] < 0) and (point[0]-bounds[rindex][fa][2] > TrueNearestFoundDelta[2]):
                    TrueNearestFound[2] = bounds[rindex][fa][2]
                    TrueNearestFoundDelta[2] = point[0]-bounds[rindex][fa][2]
                if (point[1]-bounds[rindex][fa][1] > 0) and (point[1]-bounds[rindex][fa][1] < TrueNearestFoundDelta[1]):
                    TrueNearestFound[1] = bounds[rindex][fa][1]
                    TrueNearestFoundDelta[1] = point[1]-bounds[rindex][fa][1]
                if (point[1]-bounds[rindex][fa][3] > 0) and (point[1]-bounds[rindex][fa][3] < TrueNearestFoundDelta[1]):
                    TrueNearestFound[1] = bounds[rindex][fa][3]
                    TrueNearestFoundDelta[1] = point[1]-bounds[rindex][fa][3]
                if (point[1]-bounds[rindex][fa][1] < 0) and (point[1]-bounds[rindex][fa][1] > TrueNearestFoundDelta[3]):
                    TrueNearestFound[3] = bounds[rindex][fa][1]
                    TrueNearestFoundDelta[3] = point[1]-bounds[rindex][fa][1]
                if (point[1]-bounds[rindex][fa][3] < 0) and (point[1]-bounds[rindex][fa][3] > TrueNearestFoundDelta[3]):
                    TrueNearestFound[3] = bounds[rindex][fa][3]
                    TrueNearestFoundDelta[3] = point[1]-bounds[rindex][fa][3]
            printweights.append(NearestFound[2])
            printranks.append(NearestFound[1])
        printtiles.append([((TrueNearestFound[0], TrueNearestFound[1]),
                            (TrueNearestFound[2], TrueNearestFound[3])), printranks, printweights])
        samps[:] = [smp for smp in samps if True != (
            smp[0] >= TrueNearestFound[0] and smp[0] <= TrueNearestFound[2] and smp[1] >= TrueNearestFound[1] and smp[1] <= TrueNearestFound[3])]
WIP = modcolln.pop(0)
i = -1
for bnd in bounds:
    i = i+1
    if GradientType == 1:
        if i < len(modcolln):
            WIP[bnd[0]:bnd[2], bnd[1]:bnd[3]
                ] = modcolln[i][bnd[0]:bnd[2], bnd[1]:bnd[3]]
        else:
            WIP[hires[0][0][0]:hires[0][1][0], hires[0][0][1]:hires[0][1][1]
                ] = X[hires[0][0][0]:hires[0][1][0], hires[0][0][1]:hires[0][1][1]]
    elif GradientType == 2:
        if i < len(modcolln):
            WIP[bnd[0]:bnd[2], bnd[1]:bnd[3]
                ] = modcolln[i][bnd[0]:bnd[2], bnd[1]:bnd[3]]
        else:
            WIP[hires[0][0][0]:hires[0][1][0],
                :] = X[hires[0][0][0]:hires[0][1][0], :]
    elif GradientType == 3:
        if i < len(modcolln):
            WIP[bnd[0]:bnd[2], bnd[1]:bnd[3]
                ] = modcolln[i][bnd[0]:bnd[2], bnd[1]:bnd[3]]
        else:
            WIP[:, hires[0][0][1]:hires[0][1][1]] = X[:,
                                                      hires[0][0][1]:hires[0][1][1]]
    elif GradientType == 4:
        WIP = modcolln[int(len(modcolln)/2)]
        for tile in printtiles:
            constituents = [modcolln[tile[1][i]-1] *
                            (tile[2][i]/sum(tile[2])) for i in range(len(tile[1]))]
            WIP[tile[0][0][0]:tile[0][1][0], tile[0][0][1]:tile[0][1][1]] = sum(
                constituents)[tile[0][0][0]:tile[0][1][0], tile[0][0][1]:tile[0][1][1]]
        for org in hires:
            WIP[org[0][0]:org[1][0], org[0][1]:org[1][1]
                ] = X[org[0][0]:org[1][0], org[0][1]:org[1][1]]
        break
    elif GradientType == 5:
        if i < len(modcolln):
            WIP[bnd[0]:bnd[2], bnd[1]:bnd[3]
                ] = modcolln[i][bnd[0]:bnd[2], bnd[1]:bnd[3]]
        else:
            WIP[hires[0][0][0]:hires[0][1][0], :] = X[hires[0]
                                                      [0][0]:hires[0][1][0], :]
            WIP[:, hires[0][0][1]:hires[0][1][1]] = X[:,
                                                      hires[0][0][1]:hires[0][1][1]]
imsave('Gradient Output.jpg', WIP)
imsave('OrigForGradientOverlay.jpg', X)
plt.rcParams['figure.figsize'] = [16, 6]
plt.rcParams.update({'font.size': 18})
fig, axs = plt.subplots(1, 2)
plt.set_cmap('gray')
axs[0].imshow(X)
axs[0].axis('off')
axs[1].imshow(WIP)
axs[1].axis('off')
plt.show()
