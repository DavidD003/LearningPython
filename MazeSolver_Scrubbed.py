import os
from matplotlib.image import imsave
from matplotlib.image import imread
import matplotlib.pyplot as plt
import numpy as np
import PIL
from copy import copy
from math import floor
from math import ceil
import math
from time import time
import csv
from scipy import stats
import imageio
from collections import Counter
filename = os.path.join(r'C:\Users\xx-_-\Pictures\Mazes', 'sphinx.png')
sttime = time()
zonesize = 150
tRun = True  
endfreezeframes = 10
gifsteps = 200  
frmdrn = 0.5  
stepblocks = 3
smallestgap = 10  
dg = 5  
colourthreshold = 0.25
stampsize = 0.5  
steppropn = 0.75
proxcutoff = 0.4
CustomStepSize = 0  
nPrSt = 1000  
nPrTr = 100  
PathColour = 'white'
BoundaryColour = 'black'
if PathColour == 'white':
    pthclr = [255, 255, 255]
    grypth = 255  
else:
    pthclr = PathColour
    grypth = 255  
if BoundaryColour == 'black':
    bndclr = [0, 0, 0]
    grybnd = 0  
else:
    bndclr = BoundaryColour
    grybnd = 0  
img = np.asarray(PIL.Image.open(filename))
tempsaved = False  
if len(img[0][0]) != 3:  
    rgba_image = PIL.Image.open(filename)
    rgba_image.load()
    background = PIL.Image.new("RGB", rgba_image.size, (255, 255, 255))
    background.paste(rgba_image, mask=rgba_image.split()
                     [3])  
    img = np.asarray(background)
def cutoffdist(stepdist):
    if stepdist == smallestgap:
        return smallestgap+1
    else:
        return proxcutoff*stepdist
def showWIP():
    for st in stamps:
        for rw in range(max(0, floor(st.coord[0] - round(zones[st.zoneIndex].stampsize / 2))), min(len(origimg) - 1, floor(st.coord[0] + round(zones[st.zoneIndex].stampsize / 2))), 1):
            for cl in range(max(0, floor(st.coord[1] - round(zones[st.zoneIndex].stampsize / 2))), min(len(origimg[1]) - 1, floor(st.coord[1] + round(zones[st.zoneIndex].stampsize / 2))), 1):
                if img[rw][cl] != grybnd:
                    if st.status == 0:  
                        dispimg[rw][cl] = 60
                    elif st.status == 1:  
                        dispimg[rw][cl] = 125
                    elif st.status == 2:  
                        dispimg[rw][cl] = 200
    imgplot = plt.imshow(dispimg)
    plt.show()
def clearPath(coords1, coords2, myImg):
    y1, x1 = coords1[0], coords1[1]  
    y2, x2 = coords2[0], coords2[1]
    if myImg[min(len(myImg)-1, max(0, y2))][min(len(myImg[0])-1, max(0, x2))] == grybnd or myImg[min(len(myImg)-1, max(0, y1))][min(len(myImg[0])-1, max(0, x1))] == grybnd:
        return False  
    if x1 == x2:
        dy = y2 - y1
        if dy < 0:
            ystepper = -1
        else:
            ystepper = 1
        for i in range(0, dy+2*ystepper, ystepper):
            if i != dy+ystepper:
                if myImg[min(len(myImg)-1, max(0, int(y1+i)))][int(x1)] == grybnd:
                    return False
        return True
    elif y1 == y2:
        dx = x2 - x1
        if dx < 0:
            xstepper = -1
        else:
            xstepper = 1
        for i in range(0, dx+2*xstepper, xstepper):
            if i != dx+xstepper:
                if myImg[y1][min(len(myImg[0])-1, max(0, int(x1+i)))] == grybnd:
                    return False
        return True
    else:
        sqtr1 = True
        sqtr2 = True
        ditr = True
        dy = y2-y1
        dx = x2-x1
        if dy < 0:
            ystepper = -1
        else:
            ystepper = 1
        if dx < 0:
            xstepper = -1
        else:
            xstepper = 1
        for xadj in [-1, 0, 1]:
            for yadj in [-1, 0, 1]:
                st1y, st1x, st2y, st2x = y1+yadj, x1+xadj, y2+yadj, x2+xadj
                xdisp = st2x - st1x
                ydisp = st2y - st1y
                if ydisp < 0:
                    stepper = -1
                else:
                    stepper = 1
                xPery = xdisp / ydisp
                for obl in range(0, ydisp+stepper, stepper):
                    ystep = obl
                    for xstep in [floor(ystep*xPery), ceil(ystep*xPery)]:
                        if img[min(len(img) - 1, max(0, st1y + ystep))][min(len(img[1]) - 1, max(0, st1x + xstep))] == grybnd:
                            ditr = False
                            break
                    if ditr == False:
                        break
                if xdisp < 0:
                    stepper = -1
                else:
                    stepper = 1
                yPerx = ydisp / xdisp
                for obl in range(0, xdisp+stepper, stepper):
                    xstep = obl
                    for ystep in [floor(xstep*yPerx), ceil(xstep*yPerx)]:
                        ystep = round(xstep*yPerx)
                        if img[min(len(img) - 1, max(0, st1y + ystep))][min(len(img[1]) - 1, max(0, st1x + xstep))] == grybnd:
                            ditr = False
                            break
                    if ditr == False:
                        break
                if ditr == False:
                    break
            if ditr == False:
                break
        if ditr == True:
            return True
        for i in range(0, dy+2*ystepper, ystepper):
            if i != dy+ystepper:
                if myImg[min(len(myImg)-1, max(0, int(y1+i)))][int(x1)] == grybnd:
                    sqtr1 = False
                    break  
        if sqtr1 == True:  
            for i in range(0, dx+2*xstepper, xstepper):
                if i != dx+xstepper:
                    if myImg[min(len(myImg)-1, max(0, int(y1+dy)))][min(len(myImg[0])-1, max(0, int(x1+i)))] == grybnd:
                        sqtr1 = False
                        break  
        for i in range(0, dx+2*xstepper, xstepper):
            if i != dx+xstepper:
                if myImg[min(len(myImg)-1, max(0, int(y1)))][min(len(myImg[0])-1, max(0, int(x1+i)))] == grybnd:
                    sqtr2 = False
                    break  
        if sqtr2 == True:  
            for i in range(0, dy+2*ystepper, ystepper):
                if i != dy+ystepper:
                    if myImg[min(len(myImg)-1, max(0, int(y1+i)))][min(len(myImg[0])-1, max(0, int(x1+dx)))] == grybnd:
                        sqtr2 = False
                        break  
        if (sqtr1 == True or sqtr2 == True):
            return True
        else:
            return False
def findStartEnd(myimg):
    myimg = copy(myimg)
    stpts = []
    endpts = []
    for i in range(len(myimg)):
        for j in range(len(myimg[i])):
            if myimg[i][j][0] > 200 and (myimg[i][j][1] < 100 and myimg[i][j][2] < 100):
                stpts.append([i, j])  
                myimg[i][j] = pthclr  
            elif myimg[i][j][1] > 200 and (myimg[i][j][0] < 160 and myimg[i][j][2] < 160):
                endpts.append([i, j])  
                myimg[i][j] = pthclr  
    if len(stpts) == 0:
        print('No start points could be found')
        quit()
    if len(endpts) == 0:
        print('No end points could be found')
        quit()
    start = [round(np.mean([x[0] for x in stpts])),
             round(np.mean([x[1] for x in stpts]))]
    end = [round(np.mean([x[0] for x in endpts])),
           round(np.mean([x[1] for x in endpts]))]
    return (stpts, endpts, start, end, myimg)
origimg = copy(img)
def most_frequent(myList):
    occurence_count = Counter(myList)
    return occurence_count.most_common(1)[0][0]
class zone:
    def __init__(self):
        self.stepsize = 1  
        self.stampsize = 5  
        self.stamps = []  
zones = []
zpr = ceil(len(img[0]) / zonesize)
zpc = ceil(len(img)/zonesize)
for i in range(zpr * zpc):
    zones.append(zone())
if (PathColour == 'white' and BoundaryColour == 'black') or (PathColour == 'black' and BoundaryColour == 'white'):
    stpts, endpts, start, end, img = findStartEnd(img)
    img = np.mean(img, axis=2)  
    di = len(img) * len(img[0])
    sizes = [[] for x in zones]
    ysize = [0 for x in range(len(img[0])+1)]
    for i in range(len(img)):
        xsize = 0
        for j in range(len(img[i])):
            zId = min(floor(j/zonesize), zpr-1)+zpr*floor(i/zonesize)
            if img[i][j] > round(255*colourthreshold):
                img[i][j] = 255
                xsize = xsize + 1
                ysize[j] = ysize[j]+1
            else:
                img[i][j] = 0
                if xsize > smallestgap:
                    sizes[zId].append(xsize)
                    if xsize < di and xsize > 2:
                        di = xsize
                xsize = 0  
                if ysize[j] > smallestgap:  
                    sizes[zId].append(ysize[j])
                    if ysize[j] < di:
                        di = ysize[j]
                ysize[j] = 0  
            if len(sizes[zId]) == 0:
                sizes[zId] = [smallestgap]
    for z in sizes:
        if CustomStepSize < 1:
            di = most_frequent(z)
        else:
            di = CustomStepSize
        di = ceil(di)
        sz = floor(di * stampsize)
        zones[sizes.index(z)].stepsize = di
        zones[sizes.index(z)].stampsize = ceil(di*stampsize)
    plt.set_cmap('gray')  
else:
    stpts, endpts, start, end, img = findStartEnd(img)
def stepsz():
    if prestamp == True:
        zId = min(floor(start[1] / zonesize), zpr - 1) +            zpr * floor(start[0] / zonesize)
    else:
        zId = min(floor(sourcestamp.coord[1]/zonesize),
                  zpr-1)+zpr*floor(sourcestamp.coord[0]/zonesize)
    return round(steppropn*zones[zId].stepsize)
class stamp:  
    def __init__(self, coord, stepcount=0, prev=0, direction=0, grandpID=0, ssz=0):
        self.ID = id(self)
        self.coord = [int(round(x)) for x in coord]  
        self.prevstamp = prev  
        self.stepcount = stepcount + 1
        self.TgtDistNumSteps = floor(
            (((self.coord[0] - end[0]) ** 2 + (self.coord[1] - end[1]) ** 2)**(0.5)) / stepsz())
        self.score = self.TgtDistNumSteps
        self.status = 1
        self.dirn = direction
        self.gpID = grandpID
        self.zoneIndex = min(
            floor(self.coord[1]/zonesize), zpr-1)+zpr*floor(self.coord[0]/zonesize)
        zones[min(floor(self.coord[1]/zonesize), zpr-1)+zpr *
              floor(self.coord[0] / zonesize)].stamps.append(self)
        self.steppeddist = ssz
stamps = []  
avlstamps = []  
prestamp = True
stamps.append(stamp([round(x, ndigits=None) for x in start]))
prestamp = False
avlstamps.append(stamps[0])
PathFound = False
dispimg = copy(img)  
itern = 0
finalSt = stamps[0]
while PathFound == False:  
    itern = itern + 1
    if itern % nPrSt == 0:
        print('Step Number: '+str(itern)+'. Searching...')
    rcntstamps = []  
    best = [1000000000, 0]  
    if len(avlstamps) == 0:
        showWIP()
        print('Search unsuccessful.')
        quit()
    for x in avlstamps:  
        if x.score < best[0]:
            best = [x.score, stamps[stamps.index(x)]]
    sourcestamp = best[1]  
    angles = []  
    for mystep in [max(2, round(stepsz()/i)) for i in range(1, stepblocks+1)]:
        for i in range(mystep+2):
            if img[max(0, int(sourcestamp.coord[0]-i))][int(sourcestamp.coord[1])] != grypth:
                break  
            elif i == mystep+1:
                stamps.append(stamp(coord=[max(
                    0, sourcestamp.coord[0]-i), sourcestamp.coord[1]], stepcount=sourcestamp.stepcount, prev=id(sourcestamp), direction='Up', grandpID=sourcestamp.prevstamp, ssz=mystep))
                if stamps[len(stamps) - 1].TgtDistNumSteps < 1 and clearPath(stamps[len(stamps)-1].coord, [int(round(end[0])), int(round(end[1]))], img) == True:
                    PathFound = True
                    finalSt = stamps[len(stamps) - 1]
                avlstamps.append(stamps[len(stamps) - 1])
                rcntstamps.append(stamps[len(stamps) - 1])
                angles.append(90)
        if 90 in angles:
            break  
    for mystep in [max(2, round(stepsz()/i)) for i in range(1, stepblocks+1)]:
        for i in range(mystep+2):
            if img[min(len(img)-1, sourcestamp.coord[0]+i)][sourcestamp.coord[1]] != grypth:
                break  
            elif i == mystep+1:
                stamps.append(stamp(coord=[min(len(img)-1, sourcestamp.coord[0]+i),
                                           sourcestamp.coord[1]], stepcount=sourcestamp.stepcount, prev=id(sourcestamp), direction='Down', grandpID=sourcestamp.prevstamp, ssz=mystep))
                if stamps[len(stamps) - 1].TgtDistNumSteps < 1 and clearPath(stamps[len(stamps)-1].coord, [int(round(end[0])), int(round(end[1]))], img) == True:
                    PathFound = True
                    finalSt = stamps[len(stamps) - 1]
                avlstamps.append(stamps[len(stamps) - 1])
                rcntstamps.append(stamps[len(stamps) - 1])
                angles.append(270)
        if 270 in angles:
            break  
    for mystep in [max(2, round(stepsz()/i)) for i in range(1, stepblocks+1)]:
        for i in range(mystep+2):
            if img[sourcestamp.coord[0]][min(len(img[0])-1, sourcestamp.coord[1]+i)] != grypth:
                break  
            elif i == mystep+1:
                stamps.append(stamp(coord=[sourcestamp.coord[0],
                                           min(len(img[0])-1, sourcestamp.coord[1]+i)], stepcount=sourcestamp.stepcount, prev=id(sourcestamp), direction='Right', grandpID=sourcestamp.prevstamp, ssz=mystep))
                if stamps[len(stamps) - 1].TgtDistNumSteps < 1 and clearPath(stamps[len(stamps)-1].coord, [int(round(end[0])), int(round(end[1]))], img) == True:
                    PathFound = True
                    finalSt = stamps[len(stamps) - 1]
                avlstamps.append(stamps[len(stamps) - 1])
                rcntstamps.append(stamps[len(stamps) - 1])
                angles.append(0)
                angles.append(360)
        if 360 in angles:
            break  
    for mystep in [max(2, round(stepsz()/i)) for i in range(1, stepblocks+1)]:
        for i in range(mystep+2):
            if img[sourcestamp.coord[0]][max(0, sourcestamp.coord[1]-i)] != grypth:
                break  
            elif i == mystep+1:
                stamps.append(stamp(coord=[sourcestamp.coord[0],
                                           max(0, sourcestamp.coord[1]-i)], stepcount=sourcestamp.stepcount, prev=id(sourcestamp), direction='Left', grandpID=sourcestamp.prevstamp, ssz=mystep))
                if stamps[len(stamps) - 1].TgtDistNumSteps < 1 and clearPath(stamps[len(stamps)-1].coord, [int(round(end[0])), int(round(end[1]))], img) == True:
                    PathFound = True
                    finalSt = stamps[len(stamps) - 1]
                avlstamps.append(stamps[len(stamps) - 1])
                rcntstamps.append(stamps[len(stamps) - 1])
                angles.append(180)
        if 180 in angles:
            break  
    for mystep in [max(2, round(stepsz()/i)) for i in range(1, stepblocks+1)]:
        for i in range(0, 361, dg):
            if i not in [0, 90, 180, 270, 360] and (i not in angles):
                fit = True
                for x in angles:
                    if abs(x-i) < 54 or abs((i-360)-i) < 54:
                        fit = False  
                if fit == True:
                    xdisp = round(stepsz() * math.cos(i * math.pi / 180))
                    if math.sin(i * math.pi / 180) < 0:
                        ydisp = ceil(stepsz() * -1*math.sin(i * math.pi / 180))
                    else:
                        ydisp = floor(
                            stepsz() * -1*math.sin(i * math.pi / 180))
                    if clearPath(sourcestamp.coord, [sourcestamp.coord[0]+ydisp, sourcestamp.coord[1]+xdisp], img) == True:
                        stamps.append(stamp(coord=[min(len(img)-1, max(0, sourcestamp.coord[0]+ydisp)), min(len(img[0])-1, max(
                            0, sourcestamp.coord[1]+xdisp))], stepcount=sourcestamp.stepcount, prev=id(sourcestamp), direction='obl '+str(i), grandpID=sourcestamp.prevstamp, ssz=mystep))
                        if stamps[len(stamps) - 1].TgtDistNumSteps < 1 and clearPath(stamps[len(stamps)-1].coord, [int(round(end[0])), int(round(end[1]))], img) == True:
                            PathFound = True
                            finalSt = stamps[len(stamps) - 1]
                        avlstamps.append(stamps[len(stamps) - 1])
                        rcntstamps.append(stamps[len(stamps) - 1])
                        angles.append(i)
    nbrhd = stamps.index(sourcestamp)
    avlstamps.remove(sourcestamp)
    sourcestamp.status = 0
    chckstamps = copy(rcntstamps)
    if finalSt == stamps[0]:  
        for st in chckstamps:
            zId = min(floor(st.coord[1]/zonesize),
                      zpr-1)+zpr*floor(st.coord[0]/zonesize)
            for pr in zones[zId].stamps:
                if (pr is not st) and (pr not in chckstamps) and (pr is not finalSt) and (pr is not sourcestamp)and ((st.coord[0] - pr.coord[0]) ** 2 + (st.coord[1] - pr.coord[1]) ** 2)**(0.5) <= cutoffdist(zones[st.zoneIndex].stepsize):
                    if clearPath(st.coord, pr.coord, img) == True:
                        st.status = 0
                        rcntstamps.remove(st)
                        avlstamps.remove(st)
                        stamps.remove(st)
                        zones[zId].stamps.remove(st)
                        break
    if itern % nPrSt == 0:
        print(itern)
        showWIP()
it = stamps[stamps.index(finalSt)]  
tr = 1
while it is not stamps[0]:
    it.status = 2
    tr = tr+1  
    for x in stamps:  
        if id(x) == it.prevstamp:
            st = x  
            break  
    it = st  
    if tr % nPrTr == 0:
        print('Tracing Back Solution... ' + str(tr))
trail = []
allframes = []  
pathimg = copy(origimg)  
allimg = copy(origimg)
for st in stamps:
    if st.status == 2:
        trail.append(st)  
        for rw in range(max(0, floor(st.coord[0] - (zones[st.zoneIndex].stampsize / 2))), min(len(origimg) - 1, floor(st.coord[0] + (zones[st.zoneIndex].stampsize / 2))), 1):
            for cl in range(max(0, floor(st.coord[1] - (zones[st.zoneIndex].stampsize / 2))), min(len(origimg[1]) - 1, floor(st.coord[1] + (sz / 2))), 1):
                if img[rw][cl] != grybnd:
                    allimg[rw][cl] = np.array([28, 206, 13])
        if stamps.index(st) % gifsteps == 0:
            allframes.append(copy(allimg))
        if stamps.index(st) == len(stamps) - 1:
            allframes.extend([allimg]*endfreezeframes)
    else:
        for rw in range(max(0, floor(st.coord[0] - (zones[st.zoneIndex].stampsize / 2))), min(len(origimg) - 1, floor(st.coord[0] + (zones[st.zoneIndex].stampsize / 2))), 1):
            for cl in range(max(0, floor(st.coord[1] - (zones[st.zoneIndex].stampsize / 2))), min(len(origimg[1]) - 1, floor(st.coord[1] + (zones[st.zoneIndex].stampsize / 2))), 1):
                if img[rw][cl] != grybnd:
                    if st.status == 0:  
                        origimg[rw][cl] = np.array([255, 127, 39])
                        if allimg[rw][cl].tolist() != [28, 206, 13]:
                            allimg[rw][cl] = np.array([255, 127, 39])
                    elif st.status == 1:  
                        origimg[rw][cl] = np.array([251, 235, 43])
                        if allimg[rw][cl].tolist() != [28, 206, 13]:
                            allimg[rw][cl] = np.array([251, 235, 43])
        if stamps.index(st) % gifsteps == 0:
            allframes.append(copy(allimg))
        if stamps.index(st) == len(stamps) - 1:
            allframes.extend([allimg]*endfreezeframes)
pathframes = []
for st in trail:
    for rw in range(max(0, floor(st.coord[0] - (zones[st.zoneIndex].stampsize / 2))), min(len(origimg) - 1, floor(st.coord[0] + (zones[st.zoneIndex].stampsize / 2))), 1):
        for cl in range(max(0, floor(st.coord[1] - (zones[st.zoneIndex].stampsize / 2))), min(len(origimg[1]) - 1, floor(st.coord[1] + (zones[st.zoneIndex].stampsize / 2))), 1):
            if img[rw][cl] != grybnd:
                origimg[rw][cl] = np.array([28, 206, 13])  
                pathimg[rw][cl] = np.array([28, 206, 13])  
    if trail.index(st) % gifsteps == 0:
        pathframes.append(copy(pathimg))
    if trail.index(st) == len(trail) - 1:
        pathframes.extend([pathimg]*endfreezeframes)
imageio.mimsave('SolutionSearch.gif', allframes,
                duration=frmdrn, subrectangles=True)
imageio.mimsave('Solution Direct.gif', pathframes,
                duration=frmdrn, subrectangles=True)
print('End Reached')
imgplot = plt.imshow(origimg)
plt.show()
im = PIL.Image.fromarray(origimg)
im.save("Solved.jpg")
if tRun == True:
    print('Total Run Time: '+str(time()-sttime))
