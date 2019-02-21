#!/usr/bin/python

import datetime
import operator

now = datetime.datetime.now()

print now

exon="iso1_feature_click_TPM_vals.txt"
exonlist=[]
with open(exon) as inputfile:
    for line in inputfile:
       exonlist.append(line.strip().split('\t'))

print exonlist[:10]

e1=[item[0] for item in exonlist] #isoform
e2=[float(item[3]) for item in exonlist] #TPM, or 0 (no represented reads).

print e1[:5], e2[:5]

intron="iso1_feature_total_TPM_vals.txt"
intronlist=[]
with open(intron) as inputfile:
    for line in inputfile:
       intronlist.append(line.strip().split('\t'))

print intronlist[:10]

i1=[item[0] for item in intronlist] #isoform
i2=[float(item[3]) for item in intronlist] #TPM, or 0 (no represented reads).

print i1[:5], i2[:5]



output1="comp_iso1_feature_TPM_vals.txt"
myfile1= open(output1, "w")

output2="intronless_iso1_feature_TPM_vals.txt"
myfile2= open(output2, "w")

now = datetime.datetime.now()


gene_iso = None
e_tpm = None
i_tpm = None
ratio = None
Note = None
hit = False

print "Start merging..."

for i in xrange(0, len(exonlist)):
    gene_iso = None
    e_tpm = None
    i_tpm = None
    ratio = None
    note = None
    hit = False

    for j in xrange(0, len(intronlist)):

        if e1[i] == i1[j] and e2[i] > 0:
            hit = True
            gene_iso = e1[i]
            e_tpm = e2[i]
            i_tpm = i2[j]
            ratio = i_tpm/e_tpm
            percent = (i_tpm/e_tpm)*100
            note = "normal"
            break

        elif e1[i] == i1[j] and e2[i] == 0 and i2[j] > 0:
            hit = True
            gene_iso = e1[i]
            e_tpm = e2[i]
            i_tpm = i2[j]
            ratio = 0
            percent = 0
            note = "no exonic reads"
            break

        elif e1[i] == i1[j] and e2[i] == 0 and i2[j] == 0:
            hit = True
            gene_iso = e1[i]
            e_tpm = e2[i]
            i_tpm = i2[j]
            ratio = 0
            percent = 0
            note = "no reads"
            break

    if hit == True:
        #print gene_iso, e_tpm, i_tpm, ratio, percent, note
        myfile1.write("%s\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n" % (gene_iso, e_tpm, i_tpm, ratio, percent, note))


    if hit == False:
        #print e1[i], e2[i]
        myfile2.write("%s\t%.3f\n" % (e1[i], e2[i]))

myfile1.close()
myfile2.close()

print now
