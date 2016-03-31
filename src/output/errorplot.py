#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
import os
from collections import OrderedDict

def eoc( a, b ):
    if b == 0:
        return -1
    return np.log( a / b )

colors = 'bgrcmyk'

rundirs = [ f for f in os.listdir('.') if 'run_' in f ]
rundirs.sort()

L2ss = {}
H1ss = {}

y1, y2 = [], []

for run in rundirs:
    rn = int(run.split('_')[-1])
    M = -1

    hs = []
    L2s = []
    H1s = []

    try:
        mainfile = open( run + '/main.out' )
    except:
        continue
    for line in mainfile.readlines():
        if '#' in line:
            ll = line.split()
            for l in ll:
                data = l.split(':')
                if data[0] == 'mcesfem.M':
                    M = int( data[1] )
                    print 'M = ', M

        if 'Y1:' and 'Y2:' in line:
            data = line.split()
            y1.append( float( data[1] ) )
            y2.append( float( data[3] ) )

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

        # # print results
        # if i == 0:
        #     print '    h         tau      L2 Gamma  (eoc h)'

        # if i > 0 and h >= 0:
        #     eocGammah = eoc( l2Gamma[-1], l2GammaOld ) / eoc( h, hOld )
        #     eocGammaTau = eoc( l2Gamma[-1], l2GammaOld ) / eoc( tau, tauOld )

        #     print '{0:7.4e} {1:7.4e} {2:7.4e} {3:7.4f} {4:7.4f}'.format( h, tau,
        #                                                                  l2Gamma[-1], eocGammah, eocGammaTau )
        # else:
        #     print '{0:7.4e} {1:7.4e} {2:7.4e} {3} {4}'.format( h, tau,
        #                                                        l2Gamma[-1], '  ---  ', '  ---  ' )

        if h < 0:
            print 'read to time ', time[-1], ' at level ', i

        l2Gmax = max( l2Gamma )

        hs.append( h )
        L2s.append( l2Gmax )
        H1s.append( h1Gamma[-1] )

        if M >= 0 or 1:
            if M in L2ss:
                L2ss[ M ].append( (l2Gmax,h) )
            else:
                L2ss[ M ] = [ (l2Gmax,h) ]

            if M in H1ss:
                H1ss[ M ].append( (h1Gamma[-1],h) )
            else:
                H1ss[ M ] = [ (h1Gamma[-1],h) ]

        l2GammaOld = l2Gamma[-1]
        h1GammaOld = h1Gamma[-1]
        hOld = h
        tauOld = tau

        plt.figure(5)
        plt.semilogy( time, l2Gamma, colors[i], label='Level {0}'.format(i) )

    plt.figure(1)
    plt.loglog( hs, L2s, '.-', label=run )
    plt.figure(2)
    plt.loglog( hs, H1s, '.-', label=run )

plt.figure(1)
plt.xlabel('$h$')
plt.ylabel('$L^2$ error')
plt.legend(loc='best')
plt.figure(2)
plt.xlabel('$h$')
plt.ylabel('$H^1$ error')
plt.legend(loc='best')

plt.figure(5)
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc='best')

plt.figure(3)
for i, (M, errorPairs) in enumerate(L2ss.iteritems()):
    print 'M', M

    byhdict = {}

    for e,h in errorPairs:
        plt.loglog( h, e, colors[i] + 'x' )

        if h in byhdict:
            byhdict[ h ].append( e )
        else:
            byhdict[ h ] = [e]


    hs = []
    means = []

    for h, errors in byhdict.iteritems():
        mean = sum( errors ) / float( len(errors) )

        hs.append(h)
        means.append(mean)

    hsSorted = [x for (x,y) in sorted(zip(hs,means), key=lambda pair: -pair[0])]
    meansSorted = [y for (x,y) in sorted(zip(hs,means), key=lambda pair: -pair[0])]

    meanOld = -1
    hOld = -1

    print '    h         tau      L2 Gamma  (eoc h)'
    for h, mean in zip( hsSorted, meansSorted ):
        if meanOld >= 0:
            eocGammah = eoc( mean, meanOld ) / eoc( h, hOld )
            print '{0:7.4e} {1:7.4e} {2:7.4e} {3:7.4f}'.format( h, tau,
                                                                mean, eocGammah )
        else:
            print '{0:7.4e} {1:7.4e} {2:7.4e} {3} {4}'.format( h, tau,
                                                               mean, '  ---  ', '  ---  ' )

        meanOld = mean
        hOld = h


    plt.loglog( hsSorted, meansSorted, colors[i] + '-o', label='$M = {0}$'.format(M) )

plt.xlabel('$h$')
plt.ylabel('$L^2$ error')

plt.legend( loc='best' )


plt.figure(4)
for i, (M, errorPairs) in enumerate(H1ss.iteritems()):
    print 'M', M

    byhdict = {}

    for e,h in errorPairs:
        plt.loglog( h, e, colors[i] + 'x' )

        if h in byhdict:
            byhdict[ h ].append( e )
        else:
            byhdict[ h ] = [e]


    hs = []
    means = []

    for h, errors in byhdict.iteritems():
        mean = sum( errors ) / float( len(errors) )

        hs.append(h)
        means.append(mean)

    hsSorted = [x for (x,y) in sorted(zip(hs,means), key=lambda pair: -pair[0])]
    meansSorted = [y for (x,y) in sorted(zip(hs,means), key=lambda pair: -pair[0])]

    meanOld = -1
    hOld = -1

    print '    h         tau      H1 Gamma  (eoc h)'
    for h, mean in zip( hsSorted, meansSorted ):
        if meanOld >= 0:
            eocGammah = eoc( mean, meanOld ) / eoc( h, hOld )
            print '{0:7.4e} {1:7.4e} {2:7.4e} {3:7.4f} {4}'.format( h, tau,
                                                                mean, eocGammah, len(byhdict[h]) )
        else:
            print '{0:7.4e} {1:7.4e} {2:7.4e} {3} {4}'.format( h, tau,
                                                               mean, '  ---  ', len(byhdict[h]) )

        meanOld = mean
        hOld = h


    plt.loglog( hsSorted, meansSorted, colors[i] + '-o', label='$M = {0}$'.format(M) )

plt.xlabel('$h$')
plt.ylabel('$L^2$ error')

plt.legend( loc='best' )


plt.figure(6)

# the histogram of the data
plt.title('Distribution of $Y_1$ and $Y_2$')
plt.hist( y1 + y2, 50, facecolor='green', alpha=0.5)

plt.show()

