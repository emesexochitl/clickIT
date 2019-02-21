#!/usr/bin/python

import datetime
import operator
import math
import scipy
import numpy as np
import os
import matplotlib.pyplot as plt
import lmfit
from lmfit.models import GaussianModel, ExponentialModel, VoigtModel

now = datetime.datetime.now()

print now

class gene_set():

    def __init__(self,gene_id, cluster):

        self.cluster= cluster

        self.gene_id = gene_id
        self.norm_vals = [float(y) for y in [x[1:] for x in reference if x[0] == gene_id][0]] #TPM.

        self.get_model()
        #self.half_life()
        #self.printing()
        #self.saving()

    def get_model(self):

        self.x = np.array([0,1,2,6,12,24])
        self.y = np.array(self.norm_vals)
        self.model_flag = None
        self.model_list = []

        # First model with Gaussian curve.
        self.background1 = ExponentialModel(prefix='e1_')
        self.pars1 = self.background1.guess(self.y, x=self.x)
        self.peak1 = GaussianModel(prefix='p1_')
        self.pars1 += self.peak1.guess(self.y, x=self.x)
        self.comp_mod1 = self.peak1 + self.background1
        self.init1 = self.comp_mod1.eval(self.pars1, x=self.x)
        self.comp_out1 = self.comp_mod1.fit(self.y, x=self.x, fit_kws={'nan_policy': 'omit'})
        self.comp_list1 = self.comp_out1.fit_report().split('\n')
        self.comp_chisq1 = float(self.comp_list1[6][-5:])

        # Second model with Voigt curve.
        self.background2 = ExponentialModel(prefix='e2_')
        self.pars2 = self.background2.guess(self.y, x=self.x)
        self.peak2 = VoigtModel(prefix='p2_')
        self.pars2 += self.peak2.guess(self.y, x=self.x)
        self.comp_mod2 = self.peak2 + self.background2
        self.init2 = self.comp_mod2.eval(self.pars2, x=self.x)
        self.comp_out2 = self.comp_mod2.fit(self.y, x=self.x, fit_kws={'nan_policy': 'omit'})
        self.comp_list2 = self.comp_out2.fit_report().split('\n')
        self.comp_chisq2 = float(self.comp_list2[6][-5:])

        # Exponential model for reference
        self.exp_mod = ExponentialModel(prefix='onlye_')
        self.pars = self.exp_mod.guess(self.y, x=self.x)
        self.init = self.exp_mod.eval(self.pars, x=self.x)

        self.exp_out = self.exp_mod.fit(self.y, x=self.x, missing='drop')
        self.exp_list = self.exp_out.fit_report().split('\n')
        self.exp_chisq = float(self.exp_list[6][-5:])

        self.model_list = [self.comp_chisq1, self.comp_chisq2, self.exp_chisq]

        if np.count_nonzero(np.isinf(self.comp_out1.best_fit)) == 5 and np.count_nonzero(np.isinf(self.comp_out2.best_fit)):
             model_flag = "exponential"
             self.out = self.exp_out

        elif len(self.model_list) == len(set(self.model_list)):

             if min(self.model_list) == self.comp_chisq1:
                 self.model_flag = "Gaussian compound"
                 self.out = self.comp_out1

             elif min(self.model_list) == self.comp_chisq2:
                 self.model_flag = "Voigt compound"
                 self.out = self.comp_out2

             elif min(self.model_list) == self.exp_chisq:
                 self.model_flag = "exponential"
                 self.out = self.exp_out

        elif len(self.model_list) != len(set(self.model_list)):

             if min(self.model_list) == self.comp_chisq1:
                 self.model_flag = "Gaussian compound"
                 self.out = self.comp_out1

             elif min(self.model_list) == self.comp_chisq2:
                 self.model_flag = "Voigt compound"
                 self.out = self.comp_out2

             elif min(self.model_list) == self.exp_chisq:
                 self.model_flag = "exponential"
                 self.out = self.exp_out


             if min(self.model_list) == self.comp_chisq1 and self.comp_chisq1 == self.comp_chisq2:
                 self.model_flag = "Both compounds"
                 self.out = self.comp_out2

             if min(self.model_list) == self.comp_chisq2 and self.comp_chisq2 == self.exp_chisq:
                 self.model_flag = "Voigt compound and exponential"
                 self.out = self.comp_out2

             if min(self.model_list) == self.exp_chisq and self.exp_chisq == self.comp_chisq1:
                 self.model_flag = "Gaussian compound and exponential"
                 self.out = self.comp_out1


        return self.comp_out1, self.comp_chisq1, self.comp_out2, self.comp_chisq2, self.exp_out, self.exp_chisq, self.model_flag

#    def half_life(self):
        
#        self.new_x = np.array([0])
#        self.hl_eval = np.array([0]) 
#        self.hl_array = np.array([0])
#        self.hl_coord = np.array([0])
#        self.bestfit = self.out.best_fit
#        self.idx = np.argmin(np.abs(self.bestfit - 0.5))

#        if self.idx == 5 and self.bestfit[self.idx-1] < 0.5:
#            self.bestfit = self.out.best_fit[:-1]
#            self.idx = np.argmin(np.abs(self.bestfit - 0.5))

#        if self.bestfit[self.idx] == 0.5:
#            self.half_life_y = self.bestfit[self.idx]
#            self.half_life_x = self.idx
        # New
#        elif 0.5 < self.bestfit[self.idx] and self.bestfit[self.idx] == self.bestfit[5]:
#            self.min = 0
#            self.max = 0

#        elif 0.5 > self.bestfit[self.idx] and self.bestfit[self.idx-1] > 0.5:
#            self.max = self.x[self.idx]
#            self.min = self.x[self.idx-1]
#            self.flag= 'first'
#        elif 0.5 < self.bestfit[self.idx] and self.bestfit[self.idx+1] < 0.5:
#            self.min = self.x[self.idx]
#            self.max = self.x[self.idx+1]
#            self.flag= 'second'

        #New
#        elif 0.5 < self.bestfit[self.idx] and self.bestfit[self.idx+1] > 0.5 and self.bestfit[self.idx-1] > 0.5 and  self.bestfit[self.idx-3] < 0.5:
#            self.min = self.bestfit[self.idx-2]
#            self.max = self.bestfit[self.idx-3]

#        elif 0.5 < self.bestfit[self.idx] and self.bestfit[self.idx+1] > 0.5 and  self.bestfit[self.idx+2] < 0.5:
#            self.min = self.x[self.idx+1]
#            self.max = self.x[self.idx+2]

        #New

#        elif 0.5 > self.bestfit[self.idx] and self.bestfit[self.idx-1] < 0.5 and self.bestfit[self.idx-2] > 0.5:
#            self.min = self.x[self.idx-2]
#            self.max = self.x[self.idx-1]


#        elif 0.5 > self.bestfit[self.idx] and self.bestfit[self.idx+1] < 0.5 and self.bestfit[self.idx-2] > 0.5:
#            self.min = self.x[self.idx-2]
#            self.max = self.x[self.idx]

#        if self.min == 0 and self.max == 0:
#            self.half_life_y = 0
#            self.half_life_x = 0
#        else:
#            self.ranging = np.arange(self.min, self.max, 0.001)
#            for j in np.nditer(self.ranging):

#                self.new_x = np.array([j])
#                self.hl_eval = self.out.eval(x=self.new_x)

#                if  self.hl_eval >= 0.50 and self.hl_eval <= 0.51:

#                    self.hl_array = np.append(self.hl_array,self.hl_eval)
#                    self.hl_coord = np.append(self.hl_coord,self.new_x)

#            self.half_life_id = np.argmin(np.abs(self.hl_array - 0.5))
#            self.half_life_y = self.hl_array[self.half_life_id]
#            self.half_life_x = self.hl_coord[self.half_life_id]
#            self.bestfit = self.out.best_fit  
#        return self.half_life_y, self.half_life_x

#    def printing(self):


#        if self.comp_chisq < self.exp_chisq:
#            print  (self.gene_id, self.cluster, "composite", self.comp_chisq, self.half_life_y, self.half_life_x,
#                 self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
#                 self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5])

#        elif self.comp_chisq == self.exp_chisq:
#            print  (self.gene_id,  self.cluster, "equal", self.comp_chisq, self.half_life_y, self.half_life_x,
#                 self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
#                 self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5])

#        elif self.comp_chisq > self.exp_chisq:
#            print  (self.gene_id,  self.cluster, "exponential", self.exp_chisq, self.half_life_y, self.half_life_x,
#                self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
#                self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5])
 
#    def saving(self):
#            with  open('model_fit_c3_average_filtering_compound.txt', 'a') as f:

#                if self.comp_chisq < self.exp_chisq:
#                    (f.write("%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" 
#                        % (self.gene_id, self.cluster, "composite", self.comp_chisq, self.half_life_y, self.half_life_x,
#                        self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
#                        self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5])))
 
#                elif self.comp_chisq == self.exp_chisq:
#                    (f.write("%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"
#                        % (self.gene_id, self.cluster, "equal", self.comp_chisq, self.half_life_y, self.half_life_x,
#                        self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
#                        self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5])))

#                elif self.comp_chisq > self.exp_chisq:
#                    (f.write("%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"
#                        % (self.gene_id, self.cluster, "exponential", self.exp_chisq, self.half_life_y, self.half_life_x,
#                        self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
#                        self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5])))



try:
    os.remove("model_fit_c3_average_filtering_compound.txt")
except:
    pass

cluster = 1
referencefile = "expression_norm1_vals_1.25_average_filtering.txt"
reference = []
with open(referencefile) as inputfile:
    for line in inputfile:
        reference.append(line.strip().split('\t'))


idsfile = "average_chron_cluster1_clean.txt"
#idsfile = "test3.txt"
idslist = []
with open(idsfile) as inputfile:
    for line in inputfile:
        idslist.append(line.strip().split('\t'))



ids=[item[0] for item in idslist]

print len(ids)
print "-----------"


for i in ids:
   gene_id = gene_set(i, cluster)
   #x = gene_id.x
   #y = np.array(gene_id.norm_vals)
   #y = np.array(gene_id.bestfit)
   #plt.plot(x,y)
   print gene_id.model_flag #, "Gaussian", gene_id.comp_chisq1, "Voigt", gene_id.comp_chisq2, "exponential", gene_id.exp_chisq

#plt.show()   


