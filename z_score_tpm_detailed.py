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

    def __init__(self,gene_id, listname):
         
        self.filetype = listname 
        self.features=[x[0] for x in self.filetype if x[0][0:9] == gene_id] #feature ids
        self.tpms =[float(x[3]) for x in self.filetype if x[0][0:9] == gene_id] #TPM.
        
        self.st_dev()
        self.z_score()
        #self.printing()
        #self.save_files()

    def st_dev(self):

       self.mean = sum(self.tpms)/len(self.tpms)
       self.length = len(self.tpms)
       self.diff_mean = [i - self.mean for i in self.tpms]
       self.squares = [d**2 for d in self.diff_mean]
       self.ssd = sum(self.squares)/self.length # this is the SD
       self.sd = sqrt(self.ssd) #this is the variance
       return self.sd, self.mean, self.length
   
    def z_score(self): #difference from the mean

        self.z = [((j - self.mean)/self.sd) for j in self.tpms]
        self.z_array = array(self.z)
        self.p_values = scipy.stats.norm.sf(abs(self.z_array))*2
        self.sigs = ['significant' if k < 0.05 else 'not significant' for k in numpy.nditer(self.p_values) ]
        
        return self.z, self.p_values, self.sigs

    def printing(self):
    
        for l in xrange(0, len(self.features)):
            print self.features[l], self.tpms[l], self.z[l],  self.p_values[l], self.sigs[l], self.length, self.mean#, self.levene 

    def save_files(self):
        with  open('test_z_score_feature_click_corr_TPM.txt', 'a') as f:
            for l in xrange(0, len(self.features)):
                f.write("%s\t%.3f\t%.3f\t%.3f\t%s\t%d\t%.3f\n" % (self.features[l], self.tpms[l], self.z[l],  self.p_values[l], self.sigs[l], self.length, self.mean))


class variance_comp():

    def __init__(self, gene_id):    

       self.gene_id = gene_id

       self.click_tpms = clicks.tpms
       self.click_features = clicks.features         
       self.click_z = clicks.z
       self.click_z_pval = clicks.p_values
       self.click_z_sigs = clicks.sigs
       self.gene_len = clicks.length
       self.click_mean =clicks.mean

       self.total_tpms = totals.tpms
       self.total_features = totals.features
       self.total_z = totals.z       
       self.total_z_pval = totals.p_values
       self.total_z_sigs = totals.sigs
       self.total_mean = totals.mean

       self.levene_stat()
       self.covar_stat()
       self.pearson_corr_stat()

       self.save_files1()
       self.save_files2()
    def levene_stat(self): # test for equal variances

        self.levene = scipy.stats.levene(self.click_tpms, self.total_tpms)[0]
        self.levene_p = scipy.stats.levene(self.click_tpms, self.total_tpms)[1]
        self.levene_sig = 'significant' if self.levene_p < 0.05 else 'not significant'
        return self.levene, self.levene_p, self.levene_sig

    def covar_stat(self): #test for direction

       self.covar = numpy.cov(self.click_tpms, self.total_tpms, bias=True)[0][1]

       if self.covar >= 1:
           self.covar_note = 'positive'
       if self.covar <= -1:
           self.covar_note = 'negative'
       if self.covar < 1 and  self.covar > -1:
           self.covar_note = 'near zero'

       return self.covar, self.covar_note

    def pearson_corr_stat(self): #test for correlation strength

       self.pearson_corr = scipy.stats.pearsonr(self.click_tpms, self.total_tpms)[0]
       self.pearson_p = scipy.stats.pearsonr(self.click_tpms, self.total_tpms)[1]
       self.pearson_sig = 'significant' if self.pearson_p < 0.05 else 'not significant'

       return self.pearson_corr, self.pearson_p, self.pearson_sig

    def save_files1(self):
        with  open('correlation_feature.txt', 'a') as f:
            f.write("%s\t%.3f\t%.3f\t%s\t%.3f\t%s\t%.3f\t%.3f\t%s\n" % \
            (self.gene_id, self.levene, self.levene_p, self.levene_sig, \
            self.covar, self.covar_note, \
            self.pearson_corr, self.pearson_p, self.pearson_sig))    



    def save_files2(self):
        with  open('z_score_mean_correlation_feature.txt', 'a') as f:
            for l in xrange(0, len(self.click_features)):
                f.write("%s\t%.3f\t%.3f\t%.3f\t%s\t%d\t%.3f\t%s\t%.3f\t%.3f\t%.3f\t%s\t%d\t%.3f\t%s\t%.3f\t%.3f\t%s\t%.3f\t%s\t%.3f\t%.3f\t%s\n" % \
                (self.click_features[l], self.click_tpms[l], self.click_z[l],  self.click_z_pval[l], self.click_z_sigs[l], self.gene_len, self.click_mean, \
                self.total_features[l], self.total_tpms[l], self.total_z[l],  self.total_z_pval[l], self.total_z_sigs[l], self.gene_len, self.total_mean, \
                self.gene_id, self.levene, self.levene_p, self.levene_sig, \
                self.covar, self.covar_note, \
                self.pearson_corr, self.pearson_p, self.pearson_sig))

try:
    os.remove("correlation_feature.txt")
except:
    pass

try:
    os.remove("z_score_mean_correlation_feature.txt")
except:
    pass


click_tpm="feature_click_corr_TPM_vals_0.001_sorted.txt"
click_tpmlist=[]
with open(click_tpm) as inputfile:
    for line in inputfile:
        click_tpmlist.append(line.strip().split('\t'))


total_tpm="feature_total_corr_TPM_vals_0.001_sorted.txt"
total_tpmlist=[]
with open(total_tpm) as inputfile:
    for line in inputfile:
        total_tpmlist.append(line.strip().split('\t'))

l1=[item[0][0:9] for item in click_tpmlist] #isoform edited to gene id.

now = datetime.datetime.now()
counts = dict()

transformtuple = tuple(l1)
for p in transformtuple:
   counts[p] = counts.get(p, 0) + 1

for k,v in counts.items():
    print k, v
    clicks = gene_set(k, click_tpmlist)
    totals = gene_set(k, total_tpmlist)
    variances = variance_comp(k)

    #print clicks.tpms
    #print totals.tpms
    #print variances.levene_stat()

    #print variances.covar_stat()
    #print variances.pearson_corr_stat()

