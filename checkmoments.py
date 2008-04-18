#!/usr/bin/env python

import subprocess, sys

modelfile = sys.argv[1]
forcefile = None
node = None

if len(sys.argv)>2:
	forcefile = sys.argv[2]

if len(sys.argv)>3:
	node = sys.argv[3]

sys.stderr.write("READING FORCE BALANCE...\n")

cmd = ['build/forcebalance','-m',modelfile,'-o','-FM']

if forcefile:
	cmd += ['-f',forcefile]
	
if node:
	cmd += ['-n',node]

subprocess.Popen(cmd).communicate()

sys.stderr.write("CONVERTING MODEL TO GRAPHICAL FORMAT...\n")

subprocess.Popen(['build/arc2iv',modelfile,'-o']).communicate()

sys.stderr.write("RENDERING...\n")

subprocess.Popen(['../optx/viewer','-ht','forcebalance.iv','microstranmodel.iv']).communicate()


