#!/usr/bin/python

import datetime
import math

fcfile="counts_salt_mock.txt"

fclist=[]
with open(fcfile) as inputfile:
    for line in inputfile:
        fclist.append(line.strip().split('\t'))

fc1=[item[0] for item in fclist]
fc2=[float(item[1]) for item in fclist]
fc3=[float(item[2]) for item in fclist]

readnum = 0
sum = 0
rpkm = 0
rpkm_sum = 0
for i in xrange(0, len(fclist)):
    readnum = fc3[i]
    sum = sum + readnum
print sum


outfile="counts_salt_mock_RPKM_vals.txt"
myfile= open(outfile, "w")

for j in xrange(0, len(fclist)):
    rpkm = (1000000000*fc3[j])/(sum*fc2[j])
    rpkm_sum= rpkm_sum + rpkm
    if rpkm > 0:
        myfile.write("%s\t%d\t%d\t%.3f\n" % (fc1[j], fc2[j], fc3[j], rpkm))
    else: continue

    print fc1[j], fc2[j], fc3[j], rpkm
print rpkm_sum
myfile.close()
