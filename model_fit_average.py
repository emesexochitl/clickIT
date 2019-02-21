#!/usr/bin/python

import datetime
import operator
import math
import scipy
import numpy as np
import os
import matplotlib.pyplot as plt
import lmfit
from lmfit.models import GaussianModel, ExponentialModel

now = datetime.datetime.now()

print now

class gene_set():

    def __init__(self,gene_id, cluster):

        self.cluster= cluster

        self.gene_id = gene_id
        self.norm_vals = [float(y) for y in [x[1:] for x in reference if x[0] == gene_id][0]] #TPM.

        self.get_model()
        self.half_life()
        self.printing()
        self.saving()

    def get_model(self):

        self.x = np.array([0,1,2,6,12,24])
        self.y = np.array([self.norm_vals])
        self.peak = GaussianModel(prefix='g_')
        self.background  = ExponentialModel(prefix='e_')
        self.comp_mod = self.peak + self.background
        #guess?pars?
        #pars = exp_mod.guess(y, x=x)
        #pars = gauss_mod.guess(y, x=x)
        #pars += step_mod.guess()
        #pars.update(gauss1.make_params())
        #init = mod.eval(pars, x=x)
        self.comp_out = self.comp_mod.fit(self.y, x=self.x, fit_kws={'nan_policy': 'omit'})
        self.comp_list = self.comp_out.fit_report().split('\n') 
        self.comp_chisq = float(self.comp_list[6][-5:]) 

        self.exp_mod = ExponentialModel(prefix='onlye_')
        self.exp_out = self.exp_mod.fit(self.y, x=self.x, missing='drop')        
        self.exp_list = self.exp_out.fit_report().split('\n')
        self.exp_chisq = float(self.exp_list[6][-5:])

        if np.count_nonzero(np.isinf(self.comp_out.best_fit)) == 5:
             self.out = self.exp_out

        elif self.comp_chisq < self.exp_chisq:
            self.out = self.comp_out

        elif self.comp_chisq == self.exp_chisq:
            self.out = self.comp_out

        elif self.comp_chisq > self.exp_chisq:
            self.out = self.exp_out

        return self.comp_out, self.comp_chisq, self.exp_out, self.exp_chisq, self.out

    def half_life(self):
        
        self.new_x = np.array([0])
        self.hl_eval = np.array([0]) 
        self.hl_array = np.array([0])
        self.hl_coord = np.array([0])
        self.bestfit = self.out.best_fit
        self.idx = np.argmin(np.abs(self.bestfit - 0.5))

        if self.idx == 5 and self.bestfit[self.idx-1] < 0.5:
            self.bestfit = self.out.best_fit[:-1]
            self.idx = np.argmin(np.abs(self.bestfit - 0.5))

        if self.bestfit[self.idx] == 0.5:
            self.half_life_y = self.bestfit[self.idx]
            self.half_life_x = self.idx
        # New
        elif 0.5 < self.bestfit[self.idx] and self.bestfit[self.idx] == self.bestfit[5]:
            self.min = 0
            self.max = 0

        elif 0.5 > self.bestfit[self.idx] and self.bestfit[self.idx-1] > 0.5:
            self.max = self.x[self.idx]
            self.min = self.x[self.idx-1]
            self.flag= 'first'
        elif 0.5 < self.bestfit[self.idx] and self.bestfit[self.idx+1] < 0.5:
            self.min = self.x[self.idx]
            self.max = self.x[self.idx+1]
            self.flag= 'second'

        #New
        elif 0.5 < self.bestfit[self.idx] and self.bestfit[self.idx+1] > 0.5 and self.bestfit[self.idx-1] > 0.5 and  self.bestfit[self.idx-3] < 0.5:
            self.min = self.bestfit[self.idx-2]
            self.max = self.bestfit[self.idx-3]

        elif 0.5 < self.bestfit[self.idx] and self.bestfit[self.idx+1] > 0.5 and  self.bestfit[self.idx+2] < 0.5:
            self.min = self.x[self.idx+1]
            self.max = self.x[self.idx+2]

        #New

        elif 0.5 > self.bestfit[self.idx] and self.bestfit[self.idx-1] < 0.5 and self.bestfit[self.idx-2] > 0.5:
            self.min = self.x[self.idx-2]
            self.max = self.x[self.idx-1]


        elif 0.5 > self.bestfit[self.idx] and self.bestfit[self.idx+1] < 0.5 and self.bestfit[self.idx-2] > 0.5:
            self.min = self.x[self.idx-2]
            self.max = self.x[self.idx]

        if self.min == 0 and self.max == 0:
            self.half_life_y = 0
            self.half_life_x = 0
        else:
            self.ranging = np.arange(self.min, self.max, 0.001)
            for j in np.nditer(self.ranging):

                self.new_x = np.array([j])
                self.hl_eval = self.out.eval(x=self.new_x)

                if  self.hl_eval >= 0.50 and self.hl_eval <= 0.51:

                    self.hl_array = np.append(self.hl_array,self.hl_eval)
                    self.hl_coord = np.append(self.hl_coord,self.new_x)

            self.half_life_id = np.argmin(np.abs(self.hl_array - 0.5))
            self.half_life_y = self.hl_array[self.half_life_id]
            self.half_life_x = self.hl_coord[self.half_life_id]
            self.bestfit = self.out.best_fit  
        return self.half_life_y, self.half_life_x

    def printing(self):


        if self.comp_chisq < self.exp_chisq:
            print  (self.gene_id, self.cluster, "composite", self.comp_chisq, self.half_life_y, self.half_life_x,
                 self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
                 self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5])

        elif self.comp_chisq == self.exp_chisq:
            print  (self.gene_id,  self.cluster, "equal", self.comp_chisq, self.half_life_y, self.half_life_x,
                 self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
                 self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5])

        elif self.comp_chisq > self.exp_chisq:
            print  (self.gene_id,  self.cluster, "exponential", self.exp_chisq, self.half_life_y, self.half_life_x,
                self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
                self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5])
 
    def saving(self):
            with  open('model_fit_c3_average_filtering_compound.txt', 'a') as f:

                if self.comp_chisq < self.exp_chisq:
                    (f.write("%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" 
                        % (self.gene_id, self.cluster, "composite", self.comp_chisq, self.half_life_y, self.half_life_x,
                        self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
                        self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5])))
 
                elif self.comp_chisq == self.exp_chisq:
                    (f.write("%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"
                        % (self.gene_id, self.cluster, "equal", self.comp_chisq, self.half_life_y, self.half_life_x,
                        self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
                        self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5])))

                elif self.comp_chisq > self.exp_chisq:
                    (f.write("%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"
                        % (self.gene_id, self.cluster, "exponential", self.exp_chisq, self.half_life_y, self.half_life_x,
                        self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
                        self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5])))



try:
    os.remove("model_fit_c3_average_filtering_compound.txt")
except:
    pass

cluster = 3
referencefile = "expression_norm1_vals_1.25_average_filtering.txt"
reference = []
with open(referencefile) as inputfile:
    for line in inputfile:
        reference.append(line.strip().split('\t'))


idsfile = "average_chron_cluster3_clean.txt"
#idsfile = "test.txt"
idslist = []
with open(idsfile) as inputfile:
    for line in inputfile:
        idslist.append(line.strip().split('\t'))



ids=[item[0] for item in idslist]

print len(ids)
print "-----------"


for i in ids:
   gene_id = gene_set(i, cluster)
   x = gene_id.x
   #y = np.array(gene_id.norm_vals)
   y = np.array(gene_id.bestfit)
   plt.plot(x,y)

plt.show()   


