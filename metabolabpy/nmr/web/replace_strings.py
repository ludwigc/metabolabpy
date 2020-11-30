#!/usr/bin/env python

import os
import re

files = []
rx = re.compile(r'index\.html')
for path, dnames, fnames in os.walk('.'):
    files.extend([os.path.join(path, x) for x in fnames if rx.search(x)])


for k in range(len(files)):
    print(files[k])
    f = open(files[k],'r')
    p = f.read().replace('http://./','./')
    f.close()
    p = p.replace('http://../','../')
    p = p.replace('http://beregond.bham.ac.uk/~ludwigc/tutorials/2d-data-processing-tutorial/','./2d-data-processing-tutorial/index.html')
    p = p.replace('http://beregond.bham.ac.uk/~ludwigc/tutorials/','../tutorials/index.html')
    p = p.replace('"2d-data-processing-tutorial/"','"./2d-data-processing-tutorial/index.html"')
    p = p.replace('"../"','"../index.html"')
    f = open(files[k],'w')
    f.write(p)
    f.close()

