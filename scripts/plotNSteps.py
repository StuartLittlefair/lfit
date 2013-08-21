import tempfile
import pylab as p
import numpy as n
import sys
import os
import re

f = open(sys.argv[1])
lines = n.array(f.readlines())

olines = lines[n.arange(0,len(lines),2)]
x=[]
pattern = 'accepted after (\d*)'
for line in olines:
	m = re.search(pattern,line)
	if m:
		x.append(m.group(1))

x = n.array(x).astype(int)

success = float(x.size)/olines.size
print 'Success rate = ', success

bins = n.linspace(0.0,100.0,50)
p.hist(x,bins=bins)
p.show()
