#!/usr/bin/python

import datetime
import operator
from operator import itemgetter

print "Loading inputs."

filenamefile ="filename_list.txt"
filenames=[]

with open(filenamefile) as inputfile:
    for line in inputfile:
        filenames.append(line.strip().split('\t'))

print filenames

fnames=[item[0] for item in filenames]
ftypes=[item[1] for item in filenames]

output1="read_split_types2.txt"
output2="read_split_types_norm2.txt"

myfile1= open(output1, "w")
myfile2= open(output2, "w")

myfile1.write("file\texonic\tintronic\texonic-intronic\tintergenic\n")
myfile2.write("file\texonic\tintronic\texonic-intronic\tintergenic\n")


for a in xrange(0,len(filenames)):

    inputs = "%s_%s_sorted.bam.featureCounts" % (fnames[a], ftypes[a])
    now = datetime.datetime.now()
    output1="read_split_types.txt"
    output2="read_split_types_norm.txt"

    inputlist=[]
    with open(inputs) as inputfile:
        for line in inputfile:
            inputlist.append(line.strip().split('\t'))

    in1=[item[1] for item in inputlist] # Category
    in2=[item[2] for item in inputlist] # Feature ID
    
    total = 0
    counts = {'exonic':0, 'intronic':0, 'exonic-intronic':0, 'intergenic':0}

    for i in xrange(0, len(inputlist)):
        if in1[i] == "Assigned":
            if "e" in in2[i] and "i" not in in2[i]:
                counts['exonic'] += 1
                total = total + 1
       	    if "e" not in in2[i] and "i" in in2[i]:
       	        counts['intronic'] += 1
       	       	total =	total + 1
       	    if "e" in in2[i] and "i" in in2[i]:
       	        counts['exonic-intronic'] += 1
       	       	total =	total + 1

        elif in1[i] == "Unassigned_NoFeatures":
            counts['intergenic'] += 1
       	    total = total + 1

    factor = 1000000/float(total)
    print inputs, total, factor

    count_list =[[inputs], ['exonic'], ['intronic'], ['exonic-intronic'], ['intergenic']]
    count_list_n = [[inputs], ['exonic'], ['intronic'], ['exonic-intronic'], ['intergenic']]

    counts_mod = sorted(counts.items())
    sums = 0 
    for k,v in counts_mod:
        print k, v
        sums = sums + (float(v)*factor)
        for j in xrange(0, len(count_list)):
            if k in count_list[j]:
                count_list[j].append(v)
                count_list_n[j].append(float(v)*factor)                
                
    print count_list
    print count_list_n, sums       
    
    myfile1.write("%s\t%.3f\t%.3f\t%.3f\t%.3f\n" % (count_list[0][0], count_list[1][1], count_list[2][1], count_list[3][1], count_list[4][1]))
    myfile2.write("%s\t%.3f\t%.3f\t%.3f\t%.3f\n" % (count_list_n[0][0], count_list_n[1][1], count_list_n[2][1], count_list_n[3][1], count_list_n[4][1]))

myfile1.close()
myfile2.close()
