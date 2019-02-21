#!/usr/bin/python

import datetime
import operator

now = datetime.datetime.now()

print now

transcript="counts_all_totalclick_transcript_mod.txt"
#transcript="test.txt"
transcriptlist=[]
with open(transcript) as inputfile:
    for line in inputfile:
       transcriptlist.append(line.strip().split('\t'))

print transcriptlist[:10]

tcr1=[item[0] for item in transcriptlist] #gene
tcr2=[item[1] for item in transcriptlist] #transcript isoform
tcr3=[item[2] for item in transcriptlist] #trlength
tcr4=[int(item[3]) for item in transcriptlist] #C rep1
tcr5=[int(item[4]) for item in transcriptlist] #C rep2
tcr6=[int(item[5]) for item in transcriptlist] #C rep3
tcr7=[int(item[6]) for item in transcriptlist] #T rep1
tcr8=[int(item[7]) for item in transcriptlist] #T rep2
tcr9=[int(item[8]) for item in transcriptlist] #T rep1


exon="counts_all_overlap_totalclick_mod.txt"
#exon="test.txt"
exonlist=[]
with open(exon) as inputfile:
    for line in inputfile:
       exonlist.append(line.strip().split('\t'))

print exonlist[:10]

exon1=[item[0] for item in exonlist] #gene
exon2=[item[1] for item in exonlist] #exonlength
exon3=[int(item[2]) for item in exonlist] #C rep1
exon4=[int(item[3]) for item in exonlist] #C rep2
exon5=[int(item[4]) for item in exonlist] #C rep3
exon6=[int(item[5]) for item in exonlist] #T rep1
exon7=[int(item[6]) for item in exonlist] #T rep2
exon8=[int(item[7]) for item in exonlist] #T rep1


output="Arabidopsis_thaliana.TAIR10.36_clickIT_transcript_exon_ratios.txt"
myfile = open(output, "w")

gene = None
isoform = None
trlength = None
exonlength = None
tr_c_rep1 = None
tr_c_rep2 = None
tr_c_rep3 = None
tr_t_rep1 = None
tr_t_rep2 = None
tr_t_rep3 = None
e_c_rep1 = None
e_c_rep2 = None
e_c_rep3 = None
e_t_rep1 = None
e_t_rep2 = None
e_t_rep3 = None
c1_percent = None
c2_percent = None
c3_percent = None
cav_percent = None
t1_percent = None
t2_percent = None
t3_percent = None
tav_percent = None
hit = False

print "Merging has been started..."


for i in xrange(0, len(transcriptlist)):
    gene = None
    isoform = None
    trlength = None
    exonlength = None
    tr_c_rep1 = None
    tr_c_rep2 = None
    tr_c_rep3 = None 
    tr_t_rep1 = None
    tr_t_rep2 = None
    tr_t_rep3 = None 
    e_c_rep1 = None
    e_c_rep2 = None
    e_c_rep3 = None 
    e_t_rep1 = None
    e_t_rep2 = None
    e_t_rep3 = None
    c1_percent = None
    c2_percent = None
    c3_percent = None
    cav_percent = None
    t1_percent = None
    t2_percent = None
    t3_percent = None
    tav_percent = None

    hit = False

    for j in xrange(0, len(exonlist)):
    
        if tcr1[i] == exon1[j]:
            gene = tcr1[i]
            isoform = tcr2[i]
            trlength = tcr3[i]
            exonlength = exon2[j]
            tr_c_rep1 = tcr4[i]
            tr_c_rep2 = tcr5[i]
            tr_c_rep3 = tcr6[i] 
            tr_t_rep1 = tcr7[i]
            tr_t_rep2 = tcr8[i]
            tr_t_rep3 = tcr9[i] 
            e_c_rep1 = exon3[j]+0.000001
            e_c_rep2 = exon4[j]+0.000001
            e_c_rep3 = exon5[j]+0.000001
            e_t_rep1 = exon6[j]+0.000001
            e_t_rep2 = exon7[j]+0.000001
            e_t_rep3 = exon8[j]+0.000001
            c1_percent = ((tr_c_rep1/float(e_c_rep1))*100)
            c2_percent = ((tr_c_rep2/float(e_c_rep2))*100)
            c3_percent = ((tr_c_rep3/float(e_c_rep3))*100)
            cav_percent = (((tr_c_rep1+tr_c_rep2+tr_c_rep3)/(float(e_c_rep1+e_c_rep2+e_c_rep3)))*100)
            t1_percent = ((tr_t_rep1/float(e_t_rep1))*100)
            t2_percent = ((tr_t_rep2/float(e_t_rep2))*100)
            t3_percent = ((tr_t_rep3/float(e_t_rep3))*100)
            tav_percent = (((tr_t_rep1+tr_t_rep2+tr_t_rep3)/(float(e_t_rep1+e_t_rep2+e_t_rep3)))*100)
            hit = True

            print isoform, cav_percent, tav_percent
            #print gene, isoform, trlength, exonlength, tr_c_rep1, tr_c_rep2, tr_c_rep3, \
            #tr_t_rep1, tr_t_rep2, tr_t_rep3, e_c_rep1-1, e_c_rep2-0.000001, e_c_rep3-0.000001, e_t_rep1-0.000001, e_t_rep2-0.000001, e_t_rep3-0.000001, \
            #c1_percent, c2_percent, c3_percent, cav_percent, t1_percent, t2_percent, t3_percent, tav_percent

            myfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t\
            %.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n" \
            % (gene, isoform, trlength, exonlength, tr_c_rep1, tr_c_rep2, tr_c_rep3, \
            tr_t_rep1, tr_t_rep2, tr_t_rep3, e_c_rep1-0.000001, e_c_rep2-0.000001, e_c_rep3-0.000001, \
            e_t_rep1-0.000001, e_t_rep2-0.000001, e_t_rep3-0.000001, c1_percent, c2_percent, c3_percent, \
            t1_percent, t2_percent, t3_percent, cav_percent, tav_percent))
            break

        elif tcr1[i] == exon1[j].lower():
            gene = tcr1[i]
            isoform = tcr2[i]
            trlength = tcr3[i]
            exonlength = exon2[j]
            tr_c_rep1 = tcr4[i]
            tr_c_rep2 = tcr5[i]
            tr_c_rep3 = tcr6[i]
            tr_t_rep1 = tcr7[i]
            tr_t_rep2 = tcr8[i]
            tr_t_rep3 = tcr9[i]
            e_c_rep1 = exon3[j]+0.000001
            e_c_rep2 = exon4[j]+0.000001
            e_c_rep3 = exon5[j]+0.000001
            e_t_rep1 = exon6[j]+0.000001
            e_t_rep2 = exon7[j]+0.000001
            e_t_rep3 = exon8[j]+0.000001
            c1_percent = ((tr_c_rep1/float(e_c_rep1))*100)
            c2_percent = ((tr_c_rep2/float(e_c_rep2))*100)
            c3_percent = ((tr_c_rep3/float(e_c_rep3))*100)
            cav_percent = (((tr_c_rep1+tr_c_rep2+tr_c_rep3)/(float(e_c_rep1+e_c_rep2+e_c_rep3)))*100)
            t1_percent = ((tr_t_rep1/float(e_t_rep1))*100)
            t2_percent = ((tr_t_rep2/float(e_t_rep2))*100)
            t3_percent = ((tr_t_rep3/float(e_t_rep3))*100)
            tav_percent = (((tr_t_rep1+tr_t_rep2+tr_t_rep3)/(float(e_t_rep1+e_t_rep2+e_t_rep3)))*100)

            hit = True

      	    print isoform, cav_percent,	tav_percent
            #print gene, isoform, trlength, exonlength, tr_c_rep1, tr_c_rep2, tr_c_rep3, \
            #tr_t_rep1, tr_t_rep2, tr_t_rep3, e_c_rep1-0.000001, e_c_rep2-0.000001, e_c_rep3-0.000001, e_t_rep1-0.000001, e_t_rep2-0.000001, e_t_rep3-0.000001, \
            #c1_percent, c2_percent, c3_percent, cav_percent, t1_percent, t2_percent, t3_percent, tav_percent

            myfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t\
            %.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n" \
            % (gene, isoform, trlength, exonlength, tr_c_rep1, tr_c_rep2, tr_c_rep3, \
            tr_t_rep1, tr_t_rep2, tr_t_rep3, e_c_rep1-0.000001, e_c_rep2-0.000001, e_c_rep3-0.000001, \
            e_t_rep1-0.000001, e_t_rep2-0.000001, e_t_rep3-0.000001, c1_percent, c2_percent, c3_percent, \
            t1_percent, t2_percent, t3_percent, cav_percent, tav_percent))
            break

        else:
            continue

myfile.close()

print now
