#!/usr/bin/python

import datetime
import operator
import math
import scipy
import numpy
from scipy import stats
from numpy import *
import os
import matplotlib.pyplot as plt


now = datetime.datetime.now()

print now

class gene_set():

    def __init__(self,gene_id):
        
        self.gene_id = gene_id
        self.exps = [float(y)+0.001 for y in [x[1:] for x in reference if x[0] == gene_id][0]] #TPM.

        self.av_exp =[(sum(self.exps[i:i+3])/3) for i in range(0, len(self.exps), 3)] #average transcript TPM.

        self.norm_time()
        self.limit_togo()
        self.average_togo()
        #self.printing()
        self.saving()

    def norm_time(self):

        self.base = self.av_exp[0]
        self.norm_exp = [a/self.base for a in self.av_exp[1:]]
        self.norm_vals = [1] + self.norm_exp
        return self.base, self.norm_exp, self.norm_vals

    def limit_togo(self):
    
        self.greenlight1 = 0

        for i in self.norm_vals:
            if i >= 1.25:
                self.greenlight1 = 1
        return self.greenlight1

    def average_togo(self):

        self.greenlight2 = 0

        for i in self.av_exp:
            if i == 0.001:
                self.greenlight2 = 1
        return self.greenlight2


    def printing(self):

            print  self.greenlight1, self.greenlight2, self.gene_id, self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5]
            print  self.greenlight2, self.gene_id, self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5]

        
    def saving(self):
            with  open('expression_norm1_vals_1.25_average_filtering.txt', 'a') as f:
                if self.greenlight1 == 0 and self.greenlight2 == 0:
                    f.write("%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (self.gene_id, self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5]))





try:
    os.remove("expression_norm1_vals_1.25_average_filtering.txt")
except:
    pass


### input looped in!

####I will need one original referencelist!

referencefile = "time_course_normcouns.txt"
reference = []
with open(referencefile) as inputfile:
    for line in inputfile:
        reference.append(line.strip().split('\t'))

print reference[:25]

del reference[0]

#ids=[item[0] for item in reference if "0" not in item[1:]]
ids=[item[0] for item in reference]

#idsfile = "test_curve2.txt"
#idslist = []
#with open(idsfile) as inputfile:
#    for line in inputfile:
#        idslist.append(line.strip().split('\t'))


#ids=[item[0] for item in idslist]

print len(ids)
print "-----------"

x = numpy.array([0,1,2,6,12,24])

for i in ids:
   gene_id = gene_set(i)
   y = numpy.array(gene_id.norm_vals)
   #plt.plot(x,y)

#plt.show()   
