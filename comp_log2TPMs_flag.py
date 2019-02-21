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
c2=[float(item[4]) for item in clicklist] #log2TPM, or 0 (no represented reads).
c3=[float(item[2]) for item in clicklist] # average read counts.

print c1[:5], c2[:5]

total="iso1_feature_total_TPM_vals_0.001.txt"
totallist=[]
with open(total) as inputfile:
    for line in inputfile:
       totallist.append(line.strip().split('\t'))

print totallist[:10]

t1=[item[0] for item in totallist] #isoform
t2=[float(item[4]) for item in totallist] #log2TPM, or 0 (no represented reads).
t3=[float(item[2]) for item in totallist] # average read counts.

print t1[:5], t2[:5]

output1="comp_feature_log2TPM_fc_0.001_flag_corr.txt"
myfile1= open(output1, "w")

output2="2intronless_iso1_feature_TPM_vals.txt"
myfile2= open(output2, "w")

now = datetime.datetime.now()


gene_iso = None
c_tpm = None
t_tpm = None
log2fc = None
fc =None
note = None
label = None
hit = False

print "Start merging..."

for i in xrange(0, len(clicklist)):
    gene_iso = None
    c_tpm = None
    t_tpm = None
    log2fc = None
    fc = None
    note = None
    label = "average fc"
    hit = False

    for j in xrange(0, len(totallist)):

        if c1[i] == t1[j] and c3[i] > 0.001 and t3[j] > 0.001:
            hit = True
            gene_iso = c1[i]
            c_tpm = c2[i]
            t_tpm = t2[j]
            log2fc = c_tpm-t_tpm
            fc = math.pow(2, log2fc)
            note = "normal"
            break

        elif c1[i] == t1[j] and c3[i] > 0.001 and t3[j] == 0.001:
            hit = True
            gene_iso = c1[i]
            c_tpm = c2[i]
            t_tpm = t2[j]
            log2fc = c_tpm-t_tpm
            fc = math.pow(2, log2fc)
            note = "no total reads"
            break

        elif c1[i] == t1[j] and c3[i] == 0.001 and t3[j] > 0.001:
            hit = True
            gene_iso = c1[i]
            c_tpm = c2[i]
            t_tpm = t2[j]
            log2fc = c_tpm-t_tpm
            fc = math.pow(2, log2fc)
            note = "no click reads"
            break

        elif c1[i] == t1[j] and c3[i] == 0.001 and t3[j] == 0.001:
            hit = True
            gene_iso = c1[i]
            c_tpm = c2[i]
            t_tpm = t2[j]
            log2fc = c_tpm-t_tpm
            fc = math.pow(2, log2fc)
            note = "no reads"
            break

    if fc >= 2:
        label = "more than 2-fold increase"
    if fc <= 0.5 and log2fc != 0:
        label = "more than 2-fold decrease"

    if hit == True:
        print gene_iso, c3[i], t3[j], c_tpm, t_tpm, log2fc, fc, label, note
        myfile1.write("%s\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%s\n" % (gene_iso, c_tpm, t_tpm, log2fc, fc, label, note))

    if hit == False:
        #print e1[i], e2[i]
        myfile2.write("%s\t%.3f\n" % (c1[i], c2[i]))

myfile1.close()
myfile2.close()

print now
