#!/usr/bin/env python

import subprocess, sys

modelfile = sys.argv[1]

sys.stderr.write("READING FORCE BALANCE...\n")

subprocess.Popen(['build/forcebalance','-m',modelfile,'-o','-M']).communicate()

sys.stderr.write("CONVERTING MODEL TO GRAPHICAL FORMAT...\n")

subprocess.Popen(['build/arc2iv',modelfile,'-to']).communicate()

sys.stderr.write("RENDERING...\n")

subprocess.Popen(['../optx/viewer','-ht','forcebalance.iv','microstranmodel.iv']).communicate()


