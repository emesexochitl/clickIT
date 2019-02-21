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
        self.y = np.array(self.norm_vals)
        #x = np.array([0,1,2,6,12,24])

        # Exponential model for reference
        self.exp_mod = ExponentialModel(prefix='onlye_')
        self.pars = self.exp_mod.guess(self.y, x=self.x)
        self.init = self.exp_mod.eval(self.pars, x=self.x)

        self.exp_out = self.exp_mod.fit(self.y, x=self.x, missing='drop')
        self.exp_list = self.exp_out.fit_report().split('\n')
        self.exp_chisq = float(self.exp_list[6][-5:])

        return self.exp_out, self.exp_chisq



    def half_life(self):

        self.new_x = np.array([0])
        self.hl_eval = np.array([0])
        self.hl_array = np.array([0])
        self.hl_coord = np.array([0])
        self.bestfit = self.exp_out.best_fit
        self.idx = np.argmin(np.abs(self.bestfit - 0.5))

        if self.idx == 5 and self.bestfit[self.idx-1] < 0.5:
            self.bestfit = self.exp_out.best_fit[:-1]
            self.idx = np.argmin(np.abs(self.bestfit - 0.5))

        if self.bestfit[self.idx] == 0.5:
            self.half_life_y = self.bestfit[self.idx]
            self.half_life_x = self.idx

        elif 0.5 > self.bestfit[self.idx] and self.bestfit[self.idx-1] > 0.5:
            self.max = self.x[self.idx]
            self.min = self.x[self.idx-1]

        elif 0.5 < self.bestfit[self.idx] and self.bestfit[self.idx] == self.bestfit[5]:
            self.min = 0
            self.max = 0

        elif 0.5 < self.bestfit[self.idx] and self.bestfit[self.idx+1] < 0.5:
            self.min = self.x[self.idx]
            self.max = self.x[self.idx+1]

        elif 0.5 < self.bestfit[self.idx] and self.bestfit[self.idx+1] > 0.5 and  self.bestfit[self.idx+2] < 0.5:
            self.min = self.x[self.idx+1]
            self.max = self.x[self.idx+2]

        elif 0.5 > self.bestfit[self.idx] and self.bestfit[self.idx+1] < 0.5 and self.bestfit[self.idx-2] > 0.5:
            self.min = self.x[self.idx-2]
            self.max = self.x[self.idx]

        self.ranging = np.arange(self.min, self.max, 0.001)

        if self.max > 0:

#        if self.min > 0 and self.max > 0:
            for j in np.nditer(self.ranging):

                self.new_x = np.array([j])
                self.hl_eval = self.exp_out.eval(self.exp_out.params,x=self.new_x)

                if self.hl_eval >= 0.50 and self.hl_eval <= 0.51:

                    self.hl_array = np.append(self.hl_array,self.hl_eval)
                    self.hl_coord = np.append(self.hl_coord,self.new_x)

            self.half_life_id = np.argmin(np.abs(self.hl_array - 0.5))
            self.half_life_y = self.hl_array[self.half_life_id]
            self.half_life_x = self.hl_coord[self.half_life_id]
            self.bestfit = self.exp_out.best_fit

        else:
            self.half_life_y = 0
            self.half_life_x = 0
            self.bestfit = self.exp_out.best_fit
        return self.half_life_y, self.half_life_x

    def printing(self):


        print  (self.gene_id,  self.cluster, "exponential", self.exp_chisq, self.half_life_y, self.half_life_x,
            self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
            self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5])

    def saving(self):
        with  open('model_fit_c5_average_filtering_newexp.txt', 'a') as f:

            (f.write("%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"
                % (self.gene_id, self.cluster, "exponential", self.exp_chisq, self.half_life_y, self.half_life_x,
                self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
                self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5])))


try:
    os.remove("model_fit_c5_average_filtering_newexp.txt")
except:
    pass


cluster = 5

with  open('model_fit_c5_average_filtering_newexp.txt', 'a') as f:

    (f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
        % ("Gene ID","cluster","model","chisq","half_life_y","self.half_life_x",
        "norm_val_0","norm_val_1","norm_val_2","norm_val_6","norm_val_12","norm_val_24",
        "bestfit_0","bestfit_1","bestfit_2","bestfit_6","bestfit_12","bestfit_24")))

referencefile = "expression_norm1_vals_1.25_average_filtering.txt"
reference = []
with open(referencefile) as inputfile:
    for line in inputfile:
        reference.append(line.strip().split('\t'))


idsfile = "average_chron_cluster5_clean.txt"
idslist = []
with open(idsfile) as inputfile:
    for line in inputfile:
        idslist.append(line.strip().split('\t'))



ids=[item[0] for item in idslist]

print len(ids)
print "-----------"


for i in ids:
   gene_id = gene_set(i, cluster)
   x = np.sort(np.append(gene_id.x,gene_id.half_life_x))
   #y = np.array(gene_id.norm_vals)
   #y = np.array(gene_id.bestfit)
   #plt.plot(x, gene_id.exp_out.init_fit, 'b--')
   #plt.plot(x, gene_id.exp_out.best_fit, 'b-')
   #plt.plot(x,y,'ko-')
   y1 = np.sort(np.append(gene_id.bestfit, gene_id.half_life_y))[::-1]
   #y2 = np.array(gene_id.y)
   plt.plot(x,y1,label= i)
   #plt.plot(x,y2,label=i)
   plt.legend(loc='best')

   plt.vlines(gene_id.half_life_x, 0, 1, colors='g', linestyles='dashed')

   print x, y1

plt.hlines(0.5, 0, 25, colors='r', linestyles='dashed')

plt.show()   


