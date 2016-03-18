#!/usr/bin/env python

import subprocess
from multiprocessing import Pool

subprocess.call(["rm", "-rf", "../output/run_*"])

def run( idx ):

    outdir = '../output/run_{0}'.format( idx )
    subprocess.call( ["mkdir", "-p", outdir] )

    M = 2**idx
    seed = 1729 + M
    flags = 'mcesfem.M:{0} mcesfem.seed:{1} fem.prefix:{2}'.format( M, seed, outdir )
    print 'job', idx, 'flags', flags

    outfile = open( outdir + '/main.out', 'w' )
    outfile.write( '# {0}\n'.format( flags) )

    errfile = open( outdir + '/main.err', 'w' )

    R = subprocess.Popen( ['./main' ] + flags.split(), stdout=outfile, stderr=errfile )
    R.wait()

    outfile.close()
    errfile.close()

    print 'job', idx, 'complete'

MM = 12
idx = range(MM)
p = Pool( 20 )
p.map( run, idx )
