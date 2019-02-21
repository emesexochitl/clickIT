#!/usr/bin/python

import datetime
import operator
import math

now = datetime.datetime.now()

print now

click="iso1_feature_click_TPM_vals_0.001.txt"
clicklist=[]
with open(click) as inputfile:
    for line in inputfile:
       clicklist.append(line.strip().split('\t'))

print clicklist[:10]

c1=[item[0] for item in clicklist] #isoform
c2=[float(item[3]) for item in clicklist] #TPM, or 0 (no represented reads).

print c1[:5], c2[:5]

total="iso1_feature_total_TPM_vals_0.001.txt"
totallist=[]
with open(total) as inputfile:
    for line in inputfile:
       totallist.append(line.strip().split('\t'))

print totallist[:10]

t1=[item[0] for item in totallist] #isoform
t2=[float(item[3]) for item in totallist] #TPM, or 0 (no represented reads).

print t1[:5], t2[:5]

output1="feature_corr_logTPM_vals_fc.txt"
myfile1= open(output1, "w")

output2="intronless_feature_corr_TPM_vals.txt"
myfile2= open(output2, "w")

now = datetime.datetime.now()


gene_iso = None
c_tpm = None
t_tpm = None
fc =None
note = None
label = None
hit = False

print "Start merging..."

for i in xrange(0, len(clicklist)):
    gene_iso = None
    c_tpm = None
    t_tpm = None
    fc = None
    note = None
    label = "average fc"
    hit = False

    for j in xrange(0, len(totallist)):

        if c1[i] == t1[j] and c2[i] != 0 and t2[j] != 0:
            hit = True
            gene_iso = c1[i]
            c_tpm = c2[i]
            t_tpm = t2[j]
            fc = c_tpm/t_tpm
            note = "normal"
            break


    if fc >= 2:
        label = "more than 2-fold increase"
    if fc <= 0.5 and fc != 0:
        label = "more than 2-fold decrease"

    if hit == True:
        print gene_iso, c_tpm, t_tpm, fc, label, note
        myfile1.write("%s\t%.10f\t%.10f\t%.5f\t%s\t%s\n" % (gene_iso, c_tpm, t_tpm, fc, label, note))

    if hit == False:
        #print e1[i], e2[i]
        myfile2.write("%s\t%.3f\n" % (c1[i], c2[i]))

myfile1.close()
myfile2.close()

print now
