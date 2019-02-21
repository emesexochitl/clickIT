#!/usr/bin/python

import datetime
import math

fcfile="counts_transcript_click.txt"

fclist=[]
with open(fcfile) as inputfile:
    for line in inputfile:
        fclist.append(line.strip().split(' '))

fc1=[item[0] for item in fclist]
fc2=[float(item[1]) for item in fclist]
fc3=[float(item[2]) for item in fclist]


sum = 0
tpm = 0
fpkm_sum = 0
for i in xrange(0, len(fclist)):
    rpk = float(fc3[i]/fc2[i])
    #print rpk
    sum = sum + rpk
print sum

scaling_factor = float(sum/1000000)

print scaling_factor

outfile="counts_transcript_click_TPM_vals.txt"
myfile= open(outfile, "w")

for j in xrange(0, len(fclist)):
    tpm = float(fc3[j]/fc2[j])/scaling_factor
    fpkm_sum= fpkm_sum +  float(fc3[j]/fc2[j])

    if tpm > 0:
        tpmlog2 = math.log(tpm,2)
    else:
        continue

    #print fc1[j], fc2[j], fc3[j], tpm, tpmlog2
    myfile.write("%s\t%d\t%.3f\t%.3f\t%.3f\n" % (fc1[j], fc2[j], fc3[j], tpm, tpmlog2))    
print fpkm_sum
myfile.close()
