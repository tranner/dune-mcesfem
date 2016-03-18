#!/usr/bin/env python

import subprocess
from multiprocessing import Pool

name='sin'

subprocess.call(["rm", "-rf", "../output/${name}_run-*"])

def run( idx ):

    outdir = '../output/run_{0}'.format( idx )
    subprocess.call( ["mkdir", "-p", outdir] )
    print outdir

    M = 2**idx
    seed = 1729 + M
    flags = 'mcesfem.M:{0} mcesfem.seed:{1}'.format( M, seed )
    print 'job', idx, 'flags', flags

    outfile = open( outdir + '/main.out', 'w' )
    outfile.write( '# {0}\n'.format( flags) )

    errfile = open( outdir + '/main.err', 'w' )

    R = subprocess.Popen( ['./main' ] + flags.split(), stdout=outfile, stderr=errfile )
    R.wait()

    outfile.close()
    errfile.close()

    print 'job ', idx, ' complete'

MM = 6
idx = range(MM)
p = Pool( 20 )
p.map( run, idx )

