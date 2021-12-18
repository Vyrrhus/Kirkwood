import numpy as np
import csv
import sys, os

"""
    Clean le code
    Plot a(e), a(w), e(w), etc.
    Faire un système de sauvegarde d'images qui soit correct (classé par rev, nb ast puis type : histo, plots, etc)
    Avec et sans Mars influenza btw
"""

DEFAULT_FILE = 'test.dat'

if '--file' in sys.argv:
    i = sys.argv.index('--file')
    filearg = sys.argv[i+1].split(sep='-')
    filename = filearg[0] + 'MB_' + filearg[1] + 'rev_' + filearg[2] + 'ast'
else:
    filename = DEFAULT_FILE

# Read file
datafile    = open('output/'+filename+'.dat', 'r')
reader      = csv.reader(datafile)
lines       = [line for line in reader]

# Get data
""" header au format:
    [nb Objet], [x1], [x2], [x3]...
"""
header  = lines.pop(0)
nbLines = len(lines)
nbObjs  = int(header.pop(0))
nbArgs  = len(header[1::])

# Fill arrays
idList = []
for i in range(nbObjs):
    idList.append(str(lines[0][1+i*(1+nbArgs)]))
nbIniObjs = len(idList)
data = {nb:{el: [] for el in header[1::]} for nb in idList}
time = np.zeros(nbLines)

for i in range(nbLines):
    line = lines[i]
    time[i] = line.pop(0)
    for j in range(1,nbArgs+1):
        el = header[j]
        for k in range(0,nbObjs):
            try:
                nb = line[k*(nbArgs+1)]
            except:
                nbObjs -= 1
                continue
            data[nb][el].append(line[j+k*(nbArgs+1)])

for i in range(nbIniObjs):
    for el in header[1::]:
        data[idList[i]][el] = np.array(data[idList[i]][el]).astype(float)

print("Data loaded")

# Figures
"""
    - position
    - velocity
    - elements : a(t), e(t), w(t)
"""
from matplotlib import pyplot as plt

position    = '--pos'       in sys.argv
velocity    = '--vel'       in sys.argv
elements    = '--elements'  in sys.argv

# Figures and axes
if position:
    figPos, axPos = plt.subplots()
    plt.suptitle("Trajectory on the ecliptic frame")
    axPos.axis('equal')
    axPos.set_xlabel('x [AU]')
    axPos.set_ylabel('y [AU]')

if velocity:
    figVel, axVel = plt.subplots()
    plt.suptitle("Velocity on the ecliptic frame")
    axVel.axis('equal')
    axVel.set_xlabel('Vx [AU/year]')
    axVel.set_ylabel('Vy [AU/year]')

if elements:
    figSMA, axSMA = plt.subplots()
    plt.suptitle("Semi-major axis")
    axSMA.set_ylabel('semi-major axis [AU]')

    figEcc, axEcc = plt.subplots()
    plt.suptitle("Eccentricity")
    axEcc.set_ylabel("eccentricity []")

    figArg, axArg = plt.subplots()
    plt.suptitle("True longitude of periapsis")
    axArg.set_ylabel("w [°]")

    for ax in [axSMA, axEcc, axArg]:
        ax.set_xlabel("time [year]")

# Plot
for id in idList:
    body = data[id]
    if position:
        axPos.plot(body['x'], body['y'], '.')
    if velocity:
        axVel.plot(body['vx'], body['vy'], '.')
    if elements:
        if (id == '-1'):
            name = 'Sun'
        elif (id == '-2'):
            name = 'Jupiter'
        elif (id == '-3'):
            name = 'Mars'
        elif (id == '-4'):
            name = 'Saturn'
        else:
            name = ''
        size = len(body[list(body.keys())[0]])
        t = time[:size]
        axSMA.plot(t, body['semi-major axis'],label=name)
        axEcc.plot(t, body['eccentricity'],label=name)
        axArg.plot(t, body['true longitude'],label=name)

if elements:
    for ax in [axSMA,axEcc,axArg]:
        ax.legend()

plt.show()