#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
import os

def eoc( a, b ):
    if b == 0:
        return -1
    return np.log( a / b )

rundirs = [ f for f in os.listdir('.') if 'run_' in f ]
rundirs.sort()

for run in rundirs:
    hs = []
    L2s = []
    H1s = []

    mainfile = open( run + '/main.out' )
    for line in mainfile.readlines():
        if '#' in line:
            ll = line.split()
            for l in ll:
                data = l.split(':')
                if data[0] == 'mcesfem.M':
                    M = int( data[1] )
                    print 'M = ', M

    mainfile.close()

    l2GammaOld = -1
    h1GammaOld = -1

    hOld = -1
    tauOld = -1

    files = [f for f in os.listdir( run ) if '.txt' in f]

    for i,fn in enumerate(files):
        filename = run + '/' + fn
        # try:
        my_data = np.genfromtxt( filename , delimiter='  ')

        h = -1
        tau = -1

        file = open( filename )
        lines = file.readlines()
        for line in lines:
            if '# h:' in line:
                data = line.split()
                h = float( data[-1] )
            if '# tau:' in line:
                data = line.split()
                tau = float( data[-1] )
        file.close()
        # except:
        #     print 'unable to open', filename
        #     continue

        # extract data
        time    = [ d[0] for d in my_data ]
        l2Gamma = [ d[1] for d in my_data ]
        h1Gamma = [ d[2] for d in my_data ]

        if len(time) == 0:
            continue

        # # plot
        # plt.plot( time, l2Omega, 'b.-' )
        # plt.semilogy( time, l2Gamma, 'r.-' )

        # print results
        if i == 0:
            print '    h         tau      L2 Gamma  (eoc h)'

        if i > 0 and h >= 0:
            eocGammah = eoc( l2Gamma[-1], l2GammaOld ) / eoc( h, hOld )
            eocGammaTau = eoc( l2Gamma[-1], l2GammaOld ) / eoc( tau, tauOld )

            print '{0:7.4e} {1:7.4e} {2:7.4e} {3:7.4f} {4:7.4f}'.format( h, tau,
                                                                         l2Gamma[-1], eocGammah, eocGammaTau )
        else:
            print '{0:7.4e} {1:7.4e} {2:7.4e} {3} {4}'.format( h, tau,
                                                               l2Gamma[-1], '  ---  ', '  ---  ' )

        if h < 0:
            print 'read to time ', time[-1]

        hs.append( h )
        L2s.append( l2Gamma[-1] )
        H1s.append( h1Gamma[-1] )

        l2GammaOld = l2Gamma[-1]
        h1GammaOld = h1Gamma[-1]
        hOld = h
        tauOld = tau

    plt.figure(1)
    plt.loglog( hs, L2s, label=run )
    plt.figure(2)
    plt.loglog( hs, H1s, label=run )

plt.figure(1)
plt.legend(loc='best')
plt.figure(2)
plt.legend(loc='best')

plt.show()
