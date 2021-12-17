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
datafile    = open('output3d/'+filename+'.dat', 'r')
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
data = {el:np.zeros((nbLines, nbObjs)) for el in header[1::]}
data['time'] = np.zeros(nbLines)

for i in range(nbLines):
    line = lines[i]
    data['time'][i] = line.pop(0)
    for j in range(nbObjs):
        k = 1
        # print('Obj n°', j)
        for el in header[1::]:
            data[el][i,j] = line[j*(nbArgs+1)+k]
            k+=1

def kepler(x,y,vx,vy):
    # Constantes
    AU_CONST = 149597870700
    YR_CONST = 86400 * 365.2524
    GM_CONST = 1.32712440018e+020
    GM = GM_CONST * YR_CONST**2 / AU_CONST**3

    # Elements
    r2 =  x**2 +  y**2
    v2 = vx**2 + vy**2
    if not np.all(r2):
        return (np.zeros_like(x), np.zeros_like(x))
    h  = x * vy - y * vx # orbital momentum vector
    eV = [(vy * h) / GM - x / np.sqrt(r2),
        - (vx * h) / GM - y / np.sqrt(r2)]
    e  = np.sqrt(eV[0]**2 + eV[1]**2)
    a  = 1 / (2/np.sqrt(r2) - v2/GM)
    return (a,e)

# Figures
"""
    - position + norme (distance)
    - vitesse  + norme (vitesse relative)
    - elements : a(t), e(t), w(t), E(t) 
"""
from matplotlib import pyplot as plt

state       = '--state'     in sys.argv
norm        = '--norm'      in sys.argv
elements    = '--elements'  in sys.argv
recompute   = '--recompute' in sys.argv
phase       = '--phase'     in sys.argv
hist        = '--hist'      in sys.argv
hide        = '--hide'      in sys.argv
if hide:
    START = 2
else:
    START = 0

# Figures and axes
if state:
    figPos, axPos = plt.subplots()
    plt.suptitle("Position y(x) [AU]")
    axPos.axis('equal')
    figVel, axVel = plt.subplots()
    plt.suptitle("Velocity Vy(Vx) [AU/year]")
    axVel.axis('equal')

if phase and norm:
    figPhase, axPhase = plt.subplots()
    plt.suptitle("Espace des phases (v(r))")
    axPhase.axis('equal')

if norm:
    figDist, axDist = plt.subplots()
    plt.suptitle("Distance to Sun [AU]")
    figRelV, axRelV = plt.subplots()
    plt.suptitle("Relative velocity to Sun [AU/year]")

if elements:
    k = 1
    if recompute:
        k = 2
    figA, axA = plt.subplots(k,1)
    plt.suptitle("Semi-major axis [AU]")
    figE, axE = plt.subplots(k,1)
    plt.suptitle("Eccentricity")
    # figW, axW = plt.subplots(k,1, subplot_kw={'projection': 'polar'})
    # plt.suptitle("True longitude to perihelion [rad]")

if hist:
    figHist, axHist = plt.subplots()
    # plt.suptitle('True semi-major axis histogram, {} asteroids (t=0, t={})'.format(filename.split('rev')[0], filename.split('_')[1][:-3]))
    figAE, axAE = plt.subplots()


# Plot
for obj in range(START,nbObjs):
    if state:
        axPos.plot(data ['x'][:,obj], data ['y'][:,obj], '-')
        axVel.plot(data['vx'][:,obj], data['vy'][:,obj], '-')
    if norm:
        r2 = (data ['x'][:,obj] -  data['x'][:,0])**2 + (data ['y'][:,obj] - data ['y'][:,0])**2
        v2 = (data['vx'][:,obj] - data['vx'][:,0])**2 + (data['vy'][:,obj] - data['vy'][:,0])**2
        axDist.plot(data['time'], np.sqrt(r2))
        axRelV.plot(data['time'], np.sqrt(v2))
        if phase:
            axPhase.plot(np.sqrt(r2), np.sqrt(v2))
    if elements and not recompute:
        axA.plot(data['time'], data['semi-major axis'][:,obj])
        axE.plot(data['time'], data['eccentricity'][:,obj])
        # axW.plot(data['time'], data['true longitude'][:,obj])
    if elements and recompute:
        axA[0].plot(data['time'], data['semi-major axis'][:,obj])
        axE[0].plot(data['time'], data['eccentricity'][:,obj])
        # axW[0].plot(data['time'], data['true longitude'][:,obj])

        # Recompute
        a,e = kepler( data['x'][:,obj] -  data['x'][:,0],
                      data['y'][:,obj] -  data['y'][:,0],
                     data['vx'][:,obj] - data['vx'][:,0],
                     data['vy'][:,obj] - data['vy'][:,0])
        axA[1].plot(data['time'], a)
        axE[1].plot(data['time'], e)
    
if hist:
    # Read file
    datahist   = open('output/'+filename+'.dat', 'r')
    reader      = csv.reader(datafile)
    lines       = [line for line in reader] 
    i = sys.argv.index('--hist')
    try:
        bins = int(sys.argv[i+1])
    except:
        bins = 50
    axHist[0].hist(180/np.pi*data['semi-major axis'][0,2:], bins)
    axHist[1].hist(180/np.pi*data['semi-major axis'][-1,2:], bins)
    axHist[0].set_xlim(0, 360.)
    axHist[1].set_xlim(0, 360.)
    axHist[1].annotate('{} bins'.format(bins), xy=(0.9,0.9), xycoords='figure fraction')
    plt.savefig("img/histo_SEMIMAJORAXIS_{}.png".format(bins, filename[:-4]))

plt.show()