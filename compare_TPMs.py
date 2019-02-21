#!/usr/bin/python

import datetime
import operator

now = datetime.datetime.now()

print now

transcript="counts_ol_transcript_total_TPM_vals.txt"
#transcript="test.txt"
transcriptlist=[]
with open(transcript) as inputfile:
    for line in inputfile:
       transcriptlist.append(line.strip().split('\t'))

#print transcriptlist[:10]

tr1=[item[0][0:9] for item in transcriptlist] #gene
tr2=[item[0] for item in transcriptlist] #isoform
tr3=[float(item[3]) for item in transcriptlist] #TPM, or 0 (no represented reads).

print tr1[:5], tr2[:5], tr3[:5]

gene="counts_ol_gene_total_TPM_vals.txt"
genelist=[]
with open(gene) as inputfile:
    for line in inputfile:
       genelist.append(line.strip().split('\t'))

#print genelist[:10]

g1=[item[0] for item in genelist] #gene
g2=[float(item[3]) for item in genelist] #TPM, or 0 (no represented reads).

print g1[:5], g2[:5]


output="compared_TPMs_total_ol_tr-g.txt"

myfile= open(output, "w")

now = datetime.datetime.now()


gene = None
gene_iso = None
tr_tpm = None
g_tpm = None
ratio = None 
hit = False

for i in xrange(0, len(transcriptlist)):
    gene = None
    gene_iso = None
    tr_tpm = None
    g_tpm = None
    ratio = None 
    hit = False

    for j in xrange(0, len(genelist)):

        if tr1[i] == g1[j]:
            hit = True
            gene = tr1[i]
            gene_iso = tr2[i]
            tr_tpm = tr3[i]
            g_tpm = g2[j] + 0.0000000001
            ratio = tr_tpm/g_tpm
    if hit == True:
        print gene, gene_iso, tr_tpm, g_tpm, ratio
        myfile.write("%s\t%s\t%.3f\t%.3f\t%.3f\n" % (gene, gene_iso, tr_tpm, g_tpm, ratio))        

