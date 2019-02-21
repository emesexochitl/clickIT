#!/usr/bin/python

import datetime
import operator
import math
import scipy
import numpy
from scipy import stats
from numpy import *
import os

now = datetime.datetime.now()

print now

class gene_set():
    
    def __init__(self,gene_id):
   
        self.features=[x[0] for x in tpmlist if x[0][0:9] == gene_id] #feature ids
        self.tpms =[float(x[4]) for x in tpmlist if x[0][0:9] == gene_id] #log2TPM.

        self.st_dev()
        self.z_score()
        self.printing()
        self.save_files()

    def st_dev(self):

       self.mean = sum(self.tpms)/len(self.tpms)
       self.length = len(self.tpms)
       self.diff_mean = [i - self.mean for i in self.tpms]
       self.squares = [d**2 for d in self.diff_mean]
       self.ssd = sum(self.squares)/self.length # this is the SD
       self.sd = sqrt(self.ssd) #this is the variance
       return self.sd, self.mean, self.length

   
    def z_score(self):

        self.z = [((j - self.mean)/self.sd) for j in self.tpms]
        self.z_array = array(self.z)
        self.p_values = scipy.stats.norm.sf(abs(self.z_array))*2
        self.sigs = ['significant' if k < 0.05 else 'not significant' for k in numpy.nditer(self.p_values) ]
        
        return self.z, self.p_values, self.sigs

    def printing(self):
    
        for l in xrange(0, len(self.features)):
            print self.features[l], self.tpms[l], self.z[l],  self.p_values[l], self.sigs[l], self.length, self.mean 

    def save_files(self):
        with  open('z_score_feature_click_corr_log2TPM.txt', 'a') as f:
            for l in xrange(0, len(self.features)):
                f.write("%s\t%.3f\t%.3f\t%.3f\t%s\t%d\t%.3f\n" % (self.features[l], self.tpms[l], self.z[l],  self.p_values[l], self.sigs[l], self.length, self.mean))
try:
    os.remove("z_score_feature_click_corr_log2TPM.txt")
except:
    pass

tpm="feature_click_corr_TPM_vals_0.001_sorted.txt"
tpmlist=[]
with open(tpm) as inputfile:
    for line in inputfile:
        tpmlist.append(line.strip().split('\t'))

l1=[item[0][0:9] for item in tpmlist] #isoform edited to gene id.

now = datetime.datetime.now()
counts = dict()

transformtuple = tuple(l1)
for p in transformtuple:
   counts[p] = counts.get(p, 0) + 1

for k,v in counts.items():
   gene_id = gene_set(k)
