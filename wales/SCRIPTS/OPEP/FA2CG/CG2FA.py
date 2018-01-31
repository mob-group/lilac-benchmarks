#! /usr/bin/env python
# -*- coding: UTF8 -*-

import gzip
import sys
from string import split
from math import *

# --------------------- Functions ------------------------#


def rotateB(x, y, z, sxy, a, b):
    r = sqrt(x*x + y*y)
    R = sqrt(x*x + y*y + z*z)
    st = y/r
    ct = x/r
    sf = z/R
    cf = r/R
    for i in range(a, b):
        x = (sxy[i][6])
        y = (sxy[i][7])
        z = (sxy[i][8])
        X = cf*ct*x - sf*ct*z - st*y
        Y = st*cf*x - st*sf*z + ct*y
        Z = sf*x + cf*z
        sxy[i][6] = X
        sxy[i][7] = Y
        sxy[i][8] = Z


def rotateF(x, y, z, sxy, a, b):
    r = sqrt(x*x + y*y)
    R = sqrt(x*x + y*y + z*z)
    st = y/r
    ct = x/r
    sf = z/R
    cf = r/R
    for i in range(a, b):
        x = sxy[i][6]
        y = sxy[i][7]
        z = sxy[i][8]
        X = cf*ct*x + cf*st*y + sf*z
        Y = -st*x + ct*y
        Z = -sf*ct*x - st*sf*y + cf*z
        sxy[i][6] = X
        sxy[i][7] = Y
        sxy[i][8] = Z


def translate(x, y, z, sxy, a, b):
    for i in range(a, b):
        sxy[i][6] = sxy[i][6]-x
        sxy[i][7] = sxy[i][7]-y
        sxy[i][8] = sxy[i][8]-z


def prstr(sxy):
    print "*"
    for i in range(L):
        print '%3.3f\t %3.3f\t %3.3f' % (sxy[i][6], sxy[i][7], sxy[i][8])


def toPDB(s):
    for i in range(len(s)):
        print '%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f %2s%2s' % (s[i][0], int(s[i][1]), s[i][2], ' ', s[i][3], s[i][4], int(s[i][5]), ' ', float(s[i][6]), float(s[i][7]), float(s[i][8]), float(s[i][9]), float(s[i][10]), s[i][11], ' ')


def GG_sugar(sxy, qxy, PL, i, Na):
    O5x = sxy[3][6]
    O5y = sxy[3][7]
    O5z = sxy[3][8]
    C5x = sxy[4][6]
    C5y = sxy[4][7]
    C5z = sxy[4][8]
    CAx = sxy[5][6]
    CAy = sxy[5][7]
    CAz = sxy[5][8]
    CYx = sxy[11][6]
    CYy = sxy[11][7]
    CYz = sxy[11][8]

    O5gx = qxy[1][5]
    O5gy = qxy[1][6]
    O5gz = qxy[1][7]
    C5gx = qxy[2][5]
    C5gy = qxy[2][6]
    C5gz = qxy[2][7]
    CAgx = qxy[3][5]
    CAgy = qxy[3][6]
    CAgz = qxy[3][7]
    CYgx = qxy[4][5]
    CYgy = qxy[4][6]
    CYgz = qxy[4][7]

    DO5 = sqrt((O5x-O5gx)*(O5x-O5gx)+(O5y-O5gy)*(O5y-O5gy)+(O5z-O5gz)*(O5z-O5gz))
    DC5 = sqrt((C5x-C5gx)*(C5x-C5gx)+(C5y-C5gy)*(C5y-C5gy)+(C5z-C5gz)*(C5z-C5gz))
    DCA = sqrt((CAx-CAgx)*(CAx-CAgx)+(CAy-CAgy)*(CAy-CAgy)+(CAz-CAgz)*(CAz-CAgz))
    DCY = sqrt((CYx-CYgx)*(CYx-CYgx)+(CYy-CYgy)*(CYy-CYgy)+(CYz-CYgz)*(CYz-CYgz))

    dist = sqrt(DO5*DO5+DC5*DC5+DCA*DCA+DCY*DCY)
    return dist


def GG_A(sxy, qxy, PL, i, Na):
    B1x = (sxy[12][6]+sxy[13][6]+sxy[14][6]+sxy[15][6]+sxy[21][6])/5
    B1y = (sxy[12][7]+sxy[13][7]+sxy[14][7]+sxy[15][7]+sxy[21][7])/5
    B1z = (sxy[12][8]+sxy[13][8]+sxy[14][8]+sxy[15][8]+sxy[21][8])/5
    B2x = (sxy[15][6]+sxy[16][6]+sxy[18][6]+sxy[19][6]+sxy[20][6]+sxy[21][6])/6
    B2y = (sxy[15][7]+sxy[16][7]+sxy[18][7]+sxy[19][7]+sxy[20][7]+sxy[21][7])/6
    B2z = (sxy[15][8]+sxy[16][8]+sxy[18][8]+sxy[19][8]+sxy[20][8]+sxy[21][8])/6
    G1x = qxy[5][5]
    G1y = qxy[5][6]
    G1z = qxy[5][7]
    G2x = qxy[6][5]
    G2y = qxy[6][6]
    G2z = qxy[6][7]

    DB1 = sqrt((B1x-G1x)**2+(B1y-G1y)**2+(B1z-G1z)**2)
    DB2 = sqrt((B2x-G2x)**2+(B2y-G2y)**2+(B2z-G2z)**2)

    dist = sqrt(DB1**2+DB2**2)
    return dist


def GG_G(sxy, qxy, PL, i, Na):
    B1x = (sxy[12][6]+sxy[13][6]+sxy[14][6]+sxy[15][6]+sxy[22][6])/5
    B1y = (sxy[12][7]+sxy[13][7]+sxy[14][7]+sxy[15][7]+sxy[22][7])/5
    B1z = (sxy[12][8]+sxy[13][8]+sxy[14][8]+sxy[15][8]+sxy[22][8])/5
    B2x = (sxy[15][6]+sxy[16][6]+sxy[18][6]+sxy[19][6]+sxy[21][6]+sxy[22][6])/6
    B2y = (sxy[15][7]+sxy[16][7]+sxy[18][7]+sxy[19][7]+sxy[21][7]+sxy[22][7])/6
    B2z = (sxy[15][8]+sxy[16][8]+sxy[18][8]+sxy[19][8]+sxy[21][8]+sxy[22][8])/6

    G1x = qxy[5][5]
    G1y = qxy[5][6]
    G1z = qxy[5][7]
    G2x = qxy[6][5]
    G2y = qxy[6][6]
    G2z = qxy[6][7]

    DB1 = sqrt((B1x-G1x)*(B1x-G1x)+(B1y-G1y)*(B1y-G1y)+(B1z-G1z)*(B1z-G1z))
    DB2 = sqrt((B2x-G2x)*(B2x-G2x)+(B2y-G2y)*(B2y-G2y)+(B2z-G2z)*(B2z-G2z))

    dist = sqrt(DB1*DB1+DB2*DB2)
    return dist


def GG_CU(sxy, qxy, PL, i, Na):
    B1x = (sxy[12][6]+sxy[13][6]+sxy[15][6]+sxy[16][6]+sxy[18][6]+sxy[19][6])/6
    B1y = (sxy[12][7]+sxy[13][7]+sxy[15][7]+sxy[16][7]+sxy[18][7]+sxy[19][7])/6
    B1z = (sxy[12][8]+sxy[13][8]+sxy[15][8]+sxy[16][8]+sxy[18][8]+sxy[19][8])/6
    G1x = qxy[5][5]
    G1y = qxy[5][6]
    G1z = qxy[5][7]
    DB1 = sqrt((B1x-G1x)*(B1x-G1x)+(B1y-G1y)*(B1y-G1y)+(B1z-G1z)*(B1z-G1z))

    dist = sqrt(DB1*DB1)
    return dist


def rebuildA(qxy, Sxy, PL, i, Na):
    L = len(Sxy)
    LG = len(qxy)
    wxy = []
    N = 33
    for i in range(0, L):
        sxy = []
        if(Sxy[i][2] == 'P'):
            for j in range(N):
                sxy.append(Sxy[i+j])

            x = sxy[0][6]
            y = sxy[0][7]
            z = sxy[0][8]
            translate(x, y, z, sxy, 0, N)

            # Align CA
            x = sxy[5][6]
            y = sxy[5][7]
            z = sxy[5][8]
            rotateF(x, y, z, sxy, 0, N)
            x = qxy[3][5]-qxy[0][5]
            y = qxy[3][6]-qxy[0][6]
            z = qxy[3][7]-qxy[0][7]
            rotateB(x, y, z, sxy, 0, N)

            # Center Cy and B2
            x = sxy[11][6]
            y = sxy[11][7]
            z = sxy[11][8]
            translate(x, y, z, sxy, 0, N)
            x = (sxy[15][6]+sxy[16][6]+sxy[18][6]+sxy[19][6]+sxy[20][6]+sxy[21][6])/6
            y = (sxy[15][7]+sxy[16][7]+sxy[18][7]+sxy[19][7]+sxy[20][7]+sxy[21][7])/6
            z = (sxy[15][8]+sxy[16][8]+sxy[18][8]+sxy[19][8]+sxy[20][8]+sxy[21][8])/6
            # rotateF(x, y, z, sxy, 12, 22)
            # rotateF(x, y, z, sxy, 29, N)

            # Align B2
            x = qxy[6][5] - qxy[4][5]
            y = qxy[6][6] - qxy[4][6]
            z = qxy[6][7] - qxy[4][7]
            # rotateB(x, y, z, sxy, 12, 22)
            # rotateB(x, y, z, sxy, 29, N)

            # translation ref base to base coordinates: P atoms coincide
            x = sxy[0][6]-qxy[0][5]
            y = sxy[0][7]-qxy[0][6]
            z = sxy[0][8]-qxy[0][7]
            translate(x, y, z, sxy, 0, N)

            if(i == 0):
                wxy = sxy
                dist = 100

            ndist = GG_sugar(sxy, qxy, PL, i, Na)
            if ndist < dist:
                wxy = sxy
                dist = ndist

    for i in range(0, L):
        # print 'A', i
        hxy = []
        if Sxy[i][2] == 'P':
            for j in range(N):
                hxy.append(Sxy[i+j])

            # Center Cy and B2
            x = hxy[11][6]
            y = hxy[11][7]
            z = hxy[11][8]
            translate(x, y, z, hxy, 0, N)
            x = (hxy[15][6]+hxy[16][6]+hxy[18][6]+hxy[19][6]+hxy[20][6]+hxy[21][6])/6
            y = (hxy[15][7]+hxy[16][7]+hxy[18][7]+hxy[19][7]+hxy[20][7]+hxy[21][7])/6
            z = (hxy[15][8]+hxy[16][8]+hxy[18][8]+hxy[19][8]+hxy[20][8]+hxy[21][8])/6
            rotateF(x, y, z, hxy, 12, 22)
            rotateF(x, y, z, hxy, 29, N)

            # Align B2
            x = qxy[6][5] - qxy[4][5]
            y = qxy[6][6] - qxy[4][6]
            z = qxy[6][7] - qxy[4][7]
            rotateB(x, y, z, hxy, 12, 22)
            rotateB(x, y, z, hxy, 29, N)

            # print 'hxy-B'
            # for k in hxy:
            #     print k

            # translation ref base to base coordinates: CY atoms coincide
            x = hxy[11][6]-qxy[4][5]
            y = hxy[11][7]-qxy[4][6]
            z = hxy[11][8]-qxy[4][7]
            translate(x, y, z, hxy, 0, N)

            if(i == 0):
                dist = 100

            # print 'hxy'
            # for k in hxy:
            #     print k

            ndist = GG_A(hxy, qxy, PL, i, Na)
            # print dist
            if ndist < dist:
                for i in range(12, 22):
                    wxy[i] = hxy[i]
                for i in range(29, N):
                    wxy[i] = hxy[i]
                dist = ndist

            # print 'wxy'
            # for k in wxy:
            #     print k

    return wxy


def rebuildG(qxy, Sxy, PL, i, Na):
    L = len(Sxy)
    LG = len(qxy)
    wxy = []
    N = 34
    for i in range(0, L):
        sxy = []
        if(Sxy[i][2] == 'P'):
            for j in range(N):
                sxy.append(Sxy[i+j])

            # print i
            x = sxy[0][6]
            y = sxy[0][7]
            z = sxy[0][8]
            translate(x, y, z, sxy, 0, N)

            # Align CA
            x = sxy[5][6]
            y = sxy[5][7]
            z = sxy[5][8]

            rotateF(x, y, z, sxy, 0, N)
            x = qxy[3][5]-qxy[0][5]
            y = qxy[3][6]-qxy[0][6]
            z = qxy[3][7]-qxy[0][7]
            rotateB(x, y, z, sxy, 0, N)

            # Center Cy and B2
            x = sxy[11][6]
            y = sxy[11][7]
            z = sxy[11][8]
            translate(x, y, z, sxy, 0, N)
            x = (sxy[14][6]+sxy[16][6]+sxy[18][6]+sxy[19][6]+sxy[21][6]+sxy[22][6])/6
            y = (sxy[14][7]+sxy[16][7]+sxy[18][7]+sxy[19][7]+sxy[21][7]+sxy[22][7])/6
            z = (sxy[14][8]+sxy[16][8]+sxy[18][8]+sxy[19][8]+sxy[21][8]+sxy[22][8])/6
            # rotateF(x, y, z, sxy, 12, 23)
            # rotateF(x, y, z, sxy, 29, N)

            # Align B2
            x = qxy[6][5] - qxy[4][5]
            y = qxy[6][6] - qxy[4][6]
            z = qxy[6][7] - qxy[4][7]
            # rotateB(x, y, z, sxy, 12, 23)
            # rotateB(x, y, z, sxy, 29, N)

            # translation ref base to base coordinates: P atoms coincide
            x = sxy[0][6]-qxy[0][5]
            y = sxy[0][7]-qxy[0][6]
            z = sxy[0][8]-qxy[0][7]
            translate(x, y, z, sxy, 0, N)

            if(i == 0):
                wxy = sxy
                dist = 100

            ndist = GG_sugar(sxy, qxy, PL, i, Na)
            if ndist < dist:
                wxy = sxy
                dist = ndist

            # print 'wxy - sugar', ndist
            # for k in wxy:
            #     print k

    for i in range(0, L):
        # print 'G', i
        hxy = []
        if Sxy[i][2] == 'P':
            for j in range(N):
                hxy.append(Sxy[i+j])

            # Center Cy and B2
            x = hxy[11][6]
            y = hxy[11][7]
            z = hxy[11][8]
            translate(x, y, z, hxy, 0, N)
            x = (hxy[14][6]+hxy[16][6]+hxy[18][6]+hxy[19][6]+hxy[21][6]+hxy[22][6])/6
            y = (hxy[14][7]+hxy[16][7]+hxy[18][7]+hxy[19][7]+hxy[21][7]+hxy[22][7])/6
            z = (hxy[14][8]+hxy[16][8]+hxy[18][8]+hxy[19][8]+hxy[21][8]+hxy[22][8])/6
            rotateF(x, y, z, hxy, 12, 23)
            rotateF(x, y, z, hxy, 29, N)

            # Align B2
            x = qxy[6][5] - qxy[4][5]
            y = qxy[6][6] - qxy[4][6]
            z = qxy[6][7] - qxy[4][7]
            rotateB(x, y, z, hxy, 12, 23)
            rotateB(x, y, z, hxy, 29, N)

            # print 'hxy-prima'
            # for k in hxy:
            #     print k

            # translation ref base to base coordinates: P atoms coincide
            x = hxy[11][6]-qxy[4][5]
            y = hxy[11][7]-qxy[4][6]
            z = hxy[11][8]-qxy[4][7]
            translate(x, y, z, hxy, 0, N)

            # print 'hxy'
            # for k in hxy:
            #     print k

            if i == 0:
                dist = 100

            ndist = GG_G(hxy, qxy, PL, i, Na)
            # print ndist
            if ndist < dist:
                for i in range(12, 23):
                    wxy[i] = hxy[i]
                for i in range(29, N):
                    wxy[i] = hxy[i]
                dist = ndist

            # print 'wxy - final', ndist
            # for k in wxy:
            #     print k

    return wxy


def rebuildC(qxy, Sxy, PL, i, Na):
    L = len(Sxy)
    LG = len(qxy)
    wxy = []
    N = 31
    for i in range(0, L):
        sxy = []
        if Sxy[i][2] == 'P':
            for j in range(N):
                sxy.append(Sxy[i+j])

            # print i
            x = sxy[0][6]
            y = sxy[0][7]
            z = sxy[0][8]
            translate(x, y, z, sxy, 0, N)

            # Align CA
            x = sxy[5][6]
            y = sxy[5][7]
            z = sxy[5][8]

            rotateF(x, y, z, sxy, 0, N)
            x = qxy[3][5]-qxy[0][5]
            y = qxy[3][6]-qxy[0][6]
            z = qxy[3][7]-qxy[0][7]
            rotateB(x, y, z, sxy, 0, N)

            # Center Cy and B2
            x = sxy[11][6]
            y = sxy[11][7]
            z = sxy[11][8]
            translate(x, y, z, sxy, 0, N)
            x = (sxy[12][6]+sxy[13][6]+sxy[15][6]+sxy[16][6]+sxy[18][6]+sxy[19][6])/6
            y = (sxy[12][7]+sxy[13][7]+sxy[15][7]+sxy[16][7]+sxy[18][7]+sxy[19][7])/6
            z = (sxy[12][8]+sxy[13][8]+sxy[15][8]+sxy[16][8]+sxy[18][8]+sxy[19][8])/6
            # rotateF(x, y, z, sxy, 12, 20)
            # rotateF(x, y, z, sxy, 27, N)

            # Align B2
            x = qxy[5][5] - qxy[4][5]
            y = qxy[5][6] - qxy[4][6]
            z = qxy[5][7] - qxy[4][7]
            # rotateB(x, y, z, sxy, 12, 20)
            # rotateB(x, y, z, sxy, 27, N)

            # translation ref base to base coordinates: P atoms coincide
            x = sxy[0][6]-qxy[0][5]
            y = sxy[0][7]-qxy[0][6]
            z = sxy[0][8]-qxy[0][7]
            translate(x, y, z, sxy, 0, N)

            if(i == 0):
                wxy = sxy
                dist = 100

            ndist = GG_sugar(sxy, qxy, PL, i, Na)
            if ndist < dist:
                wxy = sxy
                dist = ndist

    for i in range(0, L):
        # print 'C', i
        hxy = []
        if Sxy[i][2] == 'P':
            for j in range(N):
                hxy.append(Sxy[i+j])

            # Center Cy and B2
            x = hxy[11][6]
            y = hxy[11][7]
            z = hxy[11][8]
            translate(x, y, z, hxy, 0, N)
            x = (hxy[12][6]+hxy[13][6]+hxy[15][6]+hxy[16][6]+hxy[18][6]+hxy[19][6])/6
            y = (hxy[12][7]+hxy[13][7]+hxy[15][7]+hxy[16][7]+hxy[18][7]+hxy[19][7])/6
            z = (hxy[12][8]+hxy[13][8]+hxy[15][8]+hxy[16][8]+hxy[18][8]+hxy[19][8])/6
            rotateF(x, y, z, hxy, 12, 20)
            rotateF(x, y, z, hxy, 27, N)

            # Align B2
            x = qxy[5][5] - qxy[4][5]
            y = qxy[5][6] - qxy[4][6]
            z = qxy[5][7] - qxy[4][7]
            rotateB(x, y, z, hxy, 12, 20)
            rotateB(x, y, z, hxy, 27, N)

            # translation ref base to base coordinates: P atoms coincide
            x = hxy[11][6]-qxy[4][5]
            y = hxy[11][7]-qxy[4][6]
            z = hxy[11][8]-qxy[4][7]
            translate(x, y, z, hxy, 0, N)

            if(i == 0):
                dist = 100

            ndist = GG_CU(hxy, qxy, PL, i, Na)
            if ndist < dist:
                for i in range(12, 20):
                    wxy[i] = hxy[i]
                for i in range(27, N):
                    wxy[i] = hxy[i]
                dist = ndist

    return wxy


def rebuildU(qxy, Sxy, PL, i, Na):
    L = len(Sxy)
    LG = len(qxy)
    wxy = []
    N = 30
    for i in range(0, L):
        sxy = []
        if(Sxy[i][2] == 'P'):
            for j in range(N):
                sxy.append(Sxy[i+j])

            x = sxy[0][6]
            y = sxy[0][7]
            z = sxy[0][8]
            translate(x, y, z, sxy, 0, N)

            # Align CA
            x = sxy[5][6]
            y = sxy[5][7]
            z = sxy[5][8]

            rotateF(x, y, z, sxy, 0, N)
            x = qxy[3][5]-qxy[0][5]
            y = qxy[3][6]-qxy[0][6]
            z = qxy[3][7]-qxy[0][7]
            rotateB(x, y, z, sxy, 0, N)

            # Center Cy and B2
            x = sxy[11][6]
            y = sxy[11][7]
            z = sxy[11][8]
            translate(x, y, z, sxy, 0, N)
            x = (sxy[12][6]+sxy[13][6]+sxy[15][6]+sxy[16][6]+sxy[18][6]+sxy[19][6])/6
            y = (sxy[12][7]+sxy[13][7]+sxy[15][7]+sxy[16][7]+sxy[18][7]+sxy[19][7])/6
            z = (sxy[12][8]+sxy[13][8]+sxy[15][8]+sxy[16][8]+sxy[18][8]+sxy[19][8])/6
            # rotateF(x, y, z, sxy, 12, 20)
            # rotateF(x, y, z, sxy, 27, N)

            # Align B2
            x = qxy[5][5] - qxy[4][5]
            y = qxy[5][6] - qxy[4][6]
            z = qxy[5][7] - qxy[4][7]
            # rotateB(x, y, z, sxy, 12, 20)
            # rotateB(x, y, z, sxy, 27, N)

            # translation ref base to base coordinates: P atoms coincide
            x = sxy[0][6]-qxy[0][5]
            y = sxy[0][7]-qxy[0][6]
            z = sxy[0][8]-qxy[0][7]
            translate(x, y, z, sxy, 0, N)

            if(i == 0):
                wxy = sxy
                dist = 100

            ndist = GG_sugar(sxy, qxy, PL, i, Na)
            if ndist < dist:
                wxy = sxy
                dist = ndist

    for i in range(0, L):
        # print 'U', i
        hxy = []
        if(Sxy[i][2] == 'P'):
            for j in range(N):
                hxy.append(Sxy[i+j])

            # Center Cy and B2
            x = hxy[11][6]
            y = hxy[11][7]
            z = hxy[11][8]
            translate(x, y, z, hxy, 0, N)
            x = (hxy[12][6]+hxy[13][6]+hxy[15][6]+hxy[16][6]+hxy[18][6]+hxy[19][6])/6
            y = (hxy[12][7]+hxy[13][7]+hxy[15][7]+hxy[16][7]+hxy[18][7]+hxy[19][7])/6
            z = (hxy[12][8]+hxy[13][8]+hxy[15][8]+hxy[16][8]+hxy[18][8]+hxy[19][8])/6
            rotateF(x, y, z, hxy, 12, 20)
            rotateF(x, y, z, hxy, 27, N)

            # Align B2
            x = qxy[5][5] - qxy[4][5]
            y = qxy[5][6] - qxy[4][6]
            z = qxy[5][7] - qxy[4][7]
            rotateB(x, y, z, hxy, 12, 20)
            rotateB(x, y, z, hxy, 27, N)

            # translation ref base to base coordinates: P atoms coincide
            x = hxy[11][6]-qxy[4][5]
            y = hxy[11][7]-qxy[4][6]
            z = hxy[11][8]-qxy[4][7]
            translate(x, y, z, hxy, 0, N)

            if(i == 0):
                dist = 100

            ndist = GG_CU(hxy, qxy, PL, i, Na)
            if ndist < dist:
                for i in range(12, 20):
                    wxy[i] = hxy[i]
                for i in range(27, N):
                    wxy[i] = hxy[i]
                dist = ndist

    return wxy


def main(filename):
    # --------------------- Main Routine ------------------------ #

    Na = 0

    Gxy = []
    q = []
    for l in open(filename, "r"):
        if l[0:4] == "ATOM":
            q.append(l)
            Gxy.append(split(l))
    LG = len(Gxy)

    PL = []
    for i in range(len(Gxy)):
        # print i
        Gxy[i][5] = float(q[i][30:38])
        Gxy[i][6] = float(q[i][38:46])
        Gxy[i][7] = float(q[i][46:54])
        Na = max(int(q[i][22:26]), Na)
        if(Gxy[i][2] == 'P'):
            PL.append([Gxy[i][5], Gxy[i][6], Gxy[i][7]])

    SxyA = []
    SxyG = []
    SxyC = []
    SxyU = []
    fname_list = ["baseA.lib.gz", "baseG.lib.gz", "baseC.lib.gz", "baseU.lib.gz"]
    Sxy_list = [SxyA, SxyG, SxyC, SxyU]
    for Sxyi, fi in zip(Sxy_list, fname_list):
        for l in gzip.open(fi, 'r'):
            Sxyi.append(split(l))
            Sxyi[-1][6] = float(Sxyi[-1][6])
            Sxyi[-1][7] = float(Sxyi[-1][7])
            Sxyi[-1][8] = float(Sxyi[-1][8])

    Fxy = []
    ATOM = 1
    for i in range(2, Na+1):
        qxy = []
        for j in range(LG):
            if int(q[j][22:26]) == i:
                qxy.append(Gxy[j])

        if 'G' in qxy[0][3]:
            wxy = rebuildG(qxy, SxyG, PL, i, Na)
        elif 'A' in qxy[0][3]:
            wxy = rebuildA(qxy, SxyA, PL, i, Na)
        elif 'C' in qxy[0][3]:
            wxy = rebuildC(qxy, SxyC, PL, i, Na)
        elif 'U' in qxy[0][3]:
            wxy = rebuildU(qxy, SxyU, PL, i, Na)

    #    for k in range(len(wxy)):
    #        wxy[k][5] = i
    #        Fxy.append(wxy[k])

        for k in range(len(wxy)):
            wxy[k][1] = ATOM+k
            wxy[k][5] = i

    #    for k in wxy:
    #        print k

        ATOM = ATOM+len(wxy)
        toPDB(wxy)

    # for i in range(len(Fxy)):
    #      Fxy[i][1] = i

    # toPDB(Fxy)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        fname = "6TNA_fit.pdb"
    else:
        fname = sys.argv[1]

    main(fname)
