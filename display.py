import numpy as np
import sys, os
from matplotlib import pyplot as plt

"""
    Display all the plots from output sources.
    Options:
    --file {x-y-z}: x major bodies
                    y number of Jupiter revolutions
                    z amount of asteroids
    --state:    show state vectors (x,y,vx,vy)
    --elements: show kepler elements (semi-major axis, eccentricity, true longitude)
    --hist :    show histograms
    --save :    save figs in ./img
"""

FORMAT_FILE = "output3d/{}MB_{}rev_{}ast"
HIST_LIM    = [1.2,4.8]

# File name
if '--file' in sys.argv:
    i = sys.argv.index('--file')
    argFiles = tuple(sys.argv[i+1].split(sep='-'))
else:
    argFiles = (1,0,1)

tmax = round(int(argFiles[1]) * 11.87)

# Store data
if '--hist' in sys.argv:
    hist_name = FORMAT_FILE.format(*argFiles)+'-hist.dat'
    dataUnfiltered = np.genfromtxt(hist_name, delimiter=',')

    # Filter the asteroids on a hyperbolic trajectory or too far from the Sun
    data = dataUnfiltered[(dataUnfiltered[:,4] < 1.0) & (dataUnfiltered[:,3] <= 9)]
    dataErased = dataUnfiltered[(dataUnfiltered[:,4] >= 1.0) & (dataUnfiltered[:,3] < 9)]
    # data = dataUnfiltered[np.where(dataUnfiltered[:,4] >= 1.0, False, True)]
    # data = data[np.where(data[:,3] > 9,  False, True)]
    ini = {"SMA": data[:,0], "ecc": data[:,1], "i": data[:,2], "W": data[:,3], "w": data[:,4]}
    end = {"SMA": data[:,5], "ecc": data[:,6], "i": data[:,7], "W": data[:,8], "w": data[:,9]}
    """
        Format of -hist.dat files should be:
        {a(t=0), e(t=0), w(t=0), a(t=tmax), e(t=tmax), w(t=tmax)}
    
    """
    length = len(data)    
    figHist, axHist = plt.subplots(figsize=(12,5))
    figSMA_ECC, axSMA_ECC = plt.subplots()
    figSMA_LON, axSMA_LON = plt.subplots()
    figSMA_INC, axSMA_INC = plt.subplots()

    # HISTOGRAM
    # Bins width or amount
    try:
        binsWidth = float(sys.argv[sys.argv.index('--hist')+1])
        bins = np.arange(HIST_LIM[0], HIST_LIM[1]+binsWidth, binsWidth)
    except:
        bins = int(length/20)
        binsWidth = (min(data[:,3]) + max(data[:,3])) / bins
    
    # Histograms
    axHist.hist(end["SMA"],bins, color='steelblue', edgecolor='k', lw=0.3) # Final distribution
    # axHist.hist(data[:,0],bins, color='green', alpha=0.2) # Initial distribution

    # Resonances
    axHist.axvline(x=2.8253, color='purple',  ls='-.', lw=1)
    axHist.annotate('5:2', xy=(2.8253,1), color='purple', horizontalalignment='left', verticalalignment='top', xycoords=axHist.get_xaxis_transform(), xytext=(2,-2), textcoords='offset points')
    axHist.axvline(x=3.2785, color='orange',  ls='-.', lw=1)
    axHist.annotate('2:1', xy=(3.2785,1), color='orange', horizontalalignment='left', verticalalignment='top', xycoords=axHist.get_xaxis_transform(), xytext=(2,-2), textcoords='offset points')
    axHist.axvline(x=2.5020, color='darkgreen',  ls='-.', lw=1)
    axHist.annotate('3:1', xy=(2.5020,1), color='darkgreen', horizontalalignment='left', verticalalignment='top', xycoords=axHist.get_xaxis_transform(), xytext=(2,-2), textcoords='offset points')
    axHist.set_xlim(HIST_LIM[0],HIST_LIM[1])

    # Titles, labels
    axHist.set_title('Kirkwood gaps : asteroid belt semi-major axis histogram ({} year)'.format(tmax))
    axHist.set_xlabel('Semi-major axis [UA]')
    axHist.set_ylabel('Number of asteroids (per {} UA bins)'.format(binsWidth))
    # axHist.legend(loc='upper right')

    # e(a)
    marsApoapsis     = 1.666
    jupiterPeriapsis = 4.951

    axSMA_ECC.plot(end["SMA"], end["ecc"],'k.',ms=1)

    if '--div' in sys.argv:
        figDiv, axDiv = plt.subplots()
        axDiv.plot(ini["SMA"], ini["ecc"],'g.',ms=1)
        axDiv.plot(dataErased[:,0], dataErased[:,1],'r.',ms=1, alpha=0.5)
    else:
        axSMA_ECC.plot(ini["SMA"], ini["ecc"],'g.',ms=1, alpha=0.5)
        axSMA_ECC.plot(dataErased[:,0], dataErased[:,1],'r.',ms=1, alpha=0.5)



    # Asteroid periapsis in Mars range isoclines
    """
        r = a (1-e)
        1 - r/a = e 
    """
    SMA  = np.linspace(marsApoapsis, jupiterPeriapsis, 1000)
    eRpMars = 1 - marsApoapsis/SMA     # Isoclines Rp = 1.666
    eRaJup  = jupiterPeriapsis/SMA - 1
    if '--div' in sys.argv:
        axes = [axSMA_ECC, axDiv]
    else:
        axes = [axSMA_ECC]
    for axis in axes:
        axis.plot(SMA, eRpMars, 'r-', ms=0.1, label='Rp = {} UA'.format(marsApoapsis))
        axis.plot(SMA, eRaJup,  'b-', ms=0.1, label='Ra = {} UA'.format(jupiterPeriapsis))

        axis.axvline(x=2.8253, color='purple',  ls='-.', lw=1)
        axis.axvline(x=3.2785, color='orange',  ls='-.', lw=1)
        axis.axvline(x=2.5020, color='darkgreen',  ls='-.', lw=1)

        axis.set_title('e(a)')
        axis.set_xlabel('Semi-major axis [UA]')
        axis.set_ylabel('Eccentricity')
        axis.legend(loc='upper right')

    if '--r' in sys.argv:
        figR, axR = plt.subplots()
        rp = end["SMA"] * (1 - end["ecc"])
        ra = end["SMA"] * (1 + end["ecc"])
        axR.plot(rp, ra,'k.',ms=1)
        axR.axvline(x=marsApoapsis, color='red', ls='-', lw=1, label='Mars apoapsis')
        axR.axhline(y=jupiterPeriapsis, color='blue', ls='-', lw=1, label='Jupiter periapsis')

        axR.plot(SMA, 2*2.8253 - SMA, color='purple', ls='-.', ms=0.1)
        axR.plot(SMA, 2*3.2785 - SMA, color='orange', ls='-.', ms=0.1)
        axR.plot(SMA, 2*2.5020 - SMA, color='darkgreen', ls='-.', ms=0.1)

        axR.set_title('ra(rp)')
        axR.set_xlabel('Periapsis [UA]')
        axR.set_ylabel('Apoapsis [UA]')
        axR.legend(loc='upper right')
        axR.axis('equal')


    # a(w)
    axSMA_LON.plot(end["w"], end["SMA"],'.',ms=1)

    axSMA_LON.set_title('a(w)')
    axSMA_LON.set_xlabel('True longitude of perihelion [rad]')
    axSMA_LON.set_ylabel('Semi-major axis [UA]')

    # a(i)
    axSMA_INC.plot(end["i"], end["SMA"], '.', ms=1)

    axSMA_INC.set_title('a(i)')
    axSMA_INC.set_xlabel('Inclination [rad]')
    axSMA_INC.set_ylabel('Semi-major axis [UA]')

plt.show()