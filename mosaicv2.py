# See Readme
# Key functions written by Taco at the SMA and 
# adapted for use for the CMZoom Survey by Nimesh Patel
# Checked and updated by C. Battersby Feb. 2020 for python 3
#!/usr/bin/env python
import pyregion,sys, random, math,numpy
import pyds9
#import pyregion,ds9,sys, random, math,numpy
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import FK5
#import Ctile

# Function to check if a point is within given polygon
def inPoly(x, y, poly):
    n = len(poly)
    inside =False
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y
    return inside

# Merit function to check matching of hexagonal grid pattern with given
# polygon
def merit(x, y, poly):
    nPoints = len(x)
    nVerts = len(poly)
    maxMin = -1.0e30
    for i in range(0, nVerts):
        minDist = 1.0e30
        for j in range(0, nPoints):
            d = (x[j]-poly[i][0])**2 + (y[j]-poly[i][1])**2
            if d > (0.5*spacing)**2:
                d *= 1.0e6
            if d < minDist:
                minDist = d
        if maxMin < minDist:
            maxMin = minDist
    return maxMin

# Monte-carlo shifting and rotating to get best merit
def tile(spacing, verticies, nTries, plot):
    invcosdec = 1.0/math.cos(center.dec.radian)
    n = len(verticies)
    xmax = -1.0e30
    xmin = 1.0e30
    ymax = -1.0e30
    ymin = 1.0e30
    for i in range(0,n):
        if verticies[i][0] > xmax:
            xmax = verticies[i][0]
        if verticies[i][0] < xmin:
            xmin = verticies[i][0]
        if verticies[i][1] > ymax:
            ymax = verticies[i][1]
        if verticies[i][1] < ymin:
            ymin = verticies[i][1]
    xpattern = []
    ypattern = []
    xp = xmin-2*spacing
    while xp < xmax+10*spacing:
        yp = ymin-10*spacing
        offset = 0.0
        while yp < ymax+2*spacing:
            xpattern.append(xp+offset)
            ypattern.append(yp)
            if offset == 0.0:
                offset = 0.5*spacing
            else:
                offset = 0.0
            yp += 0.866025403784*spacing
        xp += spacing*invcosdec #to get spacing right (walker found this may16,2014)
    xSpan = xmax-xmin
    ySpan = ymax-ymin
    bestMerit = 1.0e30
    i = 0
    while ( i < nTries) or (bestMerit > 100.0):
        shiftx = random.random()*spacing
        shifty = random.random()*spacing
        angle = random.random()*2.0*math.pi
        xTestPattern = []
        yTestPattern = []
        for j in range(0, len(xpattern)):
            r = math.sqrt(xpattern[j]**2 + ypattern[j]**2)
            theta = math.atan2(ypattern[j], xpattern[j]) + angle
            x = r*math.cos(theta) + shiftx
            y = r*math.sin(theta) + shifty
            if inPoly(x, y, verticies):
                xTestPattern.append(x)
                yTestPattern.append(y)
        m = merit(xTestPattern, yTestPattern, verticies)
        if m < bestMerit:
            bestx = xTestPattern
            besty = yTestPattern
            bestMerit = m
            print(i, shiftx, shifty, angle, m)
        i+=1
    xs = [verticies[i][0] for i in range(0, n)]
    ys = [verticies[i][1] for i in range(0, n)]
    xss = [verticies[i][0] for i in range(0, n)]
    yss = [verticies[i][1] for i in range(0, n)]
    xss.append(verticies[0][0])
    yss.append(verticies[0][1])
    if plot:
        plt.plot(xs,ys,'bo')
        plt.plot(bestx, besty, 'go')
        plt.plot(xss,yss,antialiased=True)
        box = [xmin*1.1, xmax*1.1, ymin*1.1, ymax*1.1]
        plt.axis(box)
        plt.axes().set_aspect('equal')
        plt.show()
    pointings = []
    for i in range(0, len(bestx)):
        pointings.append((bestx[i], besty[i]))
    return pointings

# above three functions are from Taco.
# Following code is for visualization using ds9
d = ds9.ds9()
d.set("file Mppix.fits")  # Input fits image
primaryBeam = 52.0 #arcsec; for 230 GHz.
spacing = primaryBeam/math.sqrt(3.)/3600. # sqrt(3) for sparse array
                       # use sqrt(2) for many antennas e.g. VLA 

# an existing ds9.reg will be used...
# First set of coordinates will be used for center position.
# Remaining coordinates for vertices of input polygonal region.
regionFile = "ds9.reg"

print("Click on the central position, followed by a series of points") 
print("to define the polygonal region to be mapped. Hit enter after") 
print("creating regions file: ds9.reg.")
ans = raw_input("Hit enter to continue...")

lineNo = 0
niters = 1
verticies = []
r = pyregion.open(regionFile)
for item in r:
     vertex=FK5(r[lineNo].coord_list[0],r[lineNo].coord_list[1],\
                       unit=(u.degree,u.degree))
     x = vertex.ra.degree
     y = vertex.dec.degree
     print(lineNo+1,x,y)
     if (lineNo == 0):
         center = FK5(x,y,unit=(u.degree,u.degree))
     else:
         verticies.append((x,y))
     lineNo += 1

pointings = tile(spacing,verticies, niters, False)
#pointings = Ctile.tile(spacing,verticies, niters)

# ds9out.reg file has the final pointings for later plotting
# mosaicOffsets.txt file contains pointingNo, delta-ra (") , delta-dec (")
# for observing scripts

fileRegionsOut = open("ds9out.reg","w")
fileOffsetsOut = open("mosaicOffsets.txt","w")
zerozero = 0
centerpos = "%d %.5f %.5f\n" % (zerozero, center.ra.degree, center.dec.degree)
fileOffsetsOut.write(centerpos)

i = 1
for pointing in pointings:
     print(i, pointings[i-1])
     (ra,dec)=pointings[i-1]
     regioncmd="fk5; circle(%f,%f,%.2f\")" % (ra,dec,primaryBeam/2.0)
     xoff = (ra*3600.-center.ra.arcsec)#*math.cos(center.dec.radian)
                   #positions not distances, not required.
     yoff = (dec*3600.-center.dec.arcsec)
     offsets = "%d %.2f %.2f\n" % (i,xoff,yoff)
     regionsOut = regioncmd+"\n"
     fileRegionsOut.write(regionsOut)
     fileOffsetsOut.write(offsets)
     d.set('regions',regioncmd)
     i += 1
fileRegionsOut.close()
fileOffsetsOut.close()
