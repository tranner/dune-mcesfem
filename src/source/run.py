#!/usr/bin/env python

import subprocess
from multiprocessing import Pool

subprocess.call(["rm", "-rf", "../output/run_*"])

def run( idx ):

    outdir = '../output/run_{0}'.format( idx )
    subprocess.call( ["mkdir", "-p", outdir] )

    m = idx / 10
    r = idx % 10

    M = 2**m
    seed = 1729 + idx
    flags = 'mcesfem.M:{0} mcesfem.rng.seed:{1} fem.prefix:{2}'.format( M, seed, outdir )
    print 'job', idx, 'flags', flags

    outfile = open( outdir + '/main.out', 'w' )
    outfile.write( '# {0}\n'.format( flags) )

    errfile = open( outdir + '/main.err', 'w' )

    R = subprocess.Popen( ['./main' ] + flags.split(), stdout=outfile, stderr=errfile )
    R.wait()

    outfile.close()
    errfile.close()

    print 'job', idx, 'complete'

MM = 10
idx = range(MM*10)
p = Pool( 10 )
p.map( run, idx )
