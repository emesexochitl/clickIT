#!/usr/bin/python

import datetime
import operator
import math
import scipy
import numpy as np
import os
import matplotlib.pyplot as plt
import lmfit
from lmfit.models import GaussianModel,ExponentialModel,VoigtModel

now = datetime.datetime.now()

print now

class gene_set():

    def __init__(self,gene_id, cluster):

        self.cluster= cluster

        self.gene_id = gene_id
        self.norm_vals = [float(y) for y in [x[1:] for x in reference if x[0] == gene_id][0]] #TPM.

        self.get_model()
        self.def_peaks()
        self.model_resetting()
        self.half_life()
        #self.printing()
        self.saving()

    def get_model(self):

        self.x = np.array([0,1,2,6,12,24])
        self.y = np.array(self.norm_vals)

        # Compound model with Voigt curve.
        self.background = ExponentialModel(prefix='b_')
        self.pars = self.background.guess(self.y, x=self.x)
        self.peak = VoigtModel(prefix='p_')
        self.pars += self.peak.guess(self.y, x=self.x)
        self.comp_mod = self.peak + self.background
        self.init = self.comp_mod.eval(self.pars, x=self.x)
        self.comp_out = self.comp_mod.fit(self.y, x=self.x, fit_kws={'nan_policy': 'propagate'}) # instead of 'omit', it keeps up the zero vals.
        self.comp_list = self.comp_out.fit_report().split('\n')
        self.comp_chisq = float(self.comp_list[6][-5:])

        self.out = self.comp_out
        self.chisq = float(self.comp_list[6][-5:])
        self.usedmod = self.comp_mod
        self.model_flag = "composite (exponential+Voigt)"

        return self.comp_out, self.comp_chisq, self.out, self.chisq, self.usedmod, self.model_flag

    def def_peaks(self):

        self.bestfit = self.comp_out.best_fit
        self.idx = np.argmin(np.abs(self.bestfit - 0.5))
        self.mysort = np.argsort(np.abs(self.bestfit - 0.5))
        self.peak_flag = None

        if all(i > 0.5 for i in self.bestfit):

        # Meaning that it is never reaching the half-life, and we don't do extrapolation (not enough data points). 
            
            self.min = 0
            self.max = 0 
            self.peak_flag = "No predictable half-life"

        else:

            if self.bestfit[self.idx] == 0.5:

            # If by accident one time point hits the half-life.

                self.half_life_y = self.bestfit[self.idx]
                self.half_life_x = self.idx
                self.peak_flag = "Exact compound half-life"

            elif self.bestfit[0] > 0.5 and self.bestfit[1] < 0.5:

                self.min = self.x[0]
                self.max = self.x[1]
                self.peak_flag = "Compound"


            elif self.idx == 5 and self.bestfit[self.idx-1] < 0.5:

            # Last value crosses only

                self.max = self.x[self.idx]
                self.min = self.x[self.idx-1]
                self.peak_flag = "Compound"

            elif np.abs(self.idx-self.mysort[1]) == 1:

                if self.bestfit[self.idx] < 0.5:

                    self.min = self.x[self.idx-1]
                    self.max = self.x[self.idx]
                    self.peak_flag = "Compound"
                
                elif self.bestfit[self.idx] > 0.5:

                    self.min = self.x[self.idx]
                    self.max = self.x[self.idx+1]
                    self.peak_flag = "Compound"


            elif np.abs(self.idx-self.mysort[1]) > 1:
            
           # Meaning that the steps are not linear, there's a bump.

                if self.bestfit[self.idx] < 0.5:

                    self.min = self.x[self.idx-1]
                    self.max = self.x[self.idx]
                    self.peak_flag = "Compound"

                elif self.bestfit[self.idx] > 0.5 and self.bestfit[self.mysort[1]] < 0.5:

                    if self.bestfit[self.idx+1] < 0.5:

                        self.min = self.x[self.idx]
                        self.max = self.x[self.idx+1]
                        self.peak_flag = "Compound"

                    #resetting!!
                    else:

                        self.min = self.x[self.mysort[1]-1]
                        self.max = self.x[self.mysort[1]]
                        self.peak_flag = "Resetting"

                elif self.bestfit[self.idx] > 0.5 and self.bestfit[self.mysort[1]] > 0.5:

                    if self.bestfit[self.idx+1] < 0.5:

                        self.min = self.x[self.idx]
                        self.max = self.x[self.idx+1]
                        self.peak_flag = "Compound"

                    #resetting!!
                    elif self.bestfit[self.idx+1] > 0.5 and self.bestfit[self.mysort[1]+1] < 0.5:

                        self.min = self.x[self.mysort[1]-1]
                        self.max = self.x[self.mysort[1]]
                        self.peak_flag = "Resetting"

        return self.min, self.max, self.peak_flag, self.bestfit

    def model_resetting(self):

        if self.peak_flag != "Resetting":
            #go for the previous  method
            pass

        elif self.peak_flag == "Resetting":

# mostly for plotting, half-life needs new zeros

            self.scnd_peak = np.sort(self.bestfit)[-2]
            self.scnd_idx = np.argsort(self.bestfit)[-2]
            self.newzero = self.x[self.scnd_idx]
   
            # Cutting the new time scale, reset to 0.
            self.x2 = np.array([i-self.newzero for i in self.x[self.scnd_idx:]])
            #x2 = np.array([i for i in x[scnd_idx:]])

            # Re-normalized and cutted array
            self.y2 = np.array([i/self.y[self.scnd_idx] for i in self.y[self.scnd_idx:]])
            #newarray = myarray[scnd_idx:]

            self.exp_mod = ExponentialModel(prefix='e_')
            self.pars = self.exp_mod.guess(self.y2, x=self.x2)
            self.init = self.exp_mod.eval(self.pars, x=self.x2)

            self.exp_out = self.exp_mod.fit(self.y2, x=self.x2, missing='drop')
            self.exp_list = self.exp_out.fit_report().split('\n')
            self.exp_chisq = float(self.exp_list[6][-5:])

            self.out = self.exp_out
            self.chisq = float(self.exp_list[6][-5:])
            self.usedmod = self.exp_mod

            self.bestfit = self.exp_out.best_fit
            self.idx = np.argmin(np.abs(self.bestfit - 0.5))
            self.mysort = np.argsort(np.abs(self.bestfit - 0.5))
            self.peak_flag = None

            if self.bestfit[self.idx] < 0.5:

                self.min = self.x2[self.idx-1]
                self.max = self.x2[self.idx]
                self.peak_flag = "Resetted exponential"

            elif self.bestfit[self.idx] > 0.5:

                self.min = self.x2[self.idx]
                self.max = self.x2[self.idx+1]
                self.peak_flag = "Resetted exponential"

            # For printing.
            if len(self.bestfit) < 6:

                l = [self.bestfit[-1]]*(6-len(self.bestfit))
                self.bestfit = np.append(self.bestfit,l)

            self.x = self.x2
            self.y = self.y2
            self.model_flag = "exponential"


        return self.min, self.max, self.peak_flag, self.bestfit, self.x, self.out, self.chisq, self.usedmod, self.model_flag

    def half_life(self):

        self.new_x = np.array([0])
        self.hl_eval = np.array([0])
        self.hl_array = np.array([0])
        self.hl_coord = np.array([0])
        self.step = None

        if  self.max == 0:
            self.half_life_y = 0
            self.half_life_x = 0
            self.peak_flag = "No predictable half-life"

        else:
            self.half_life_y = 0
            self.half_life_x = 0
            self.step = 0.1
            self.max_allowed = 3 
            self.attempt = 0

            #while self.attempt < 3 or self.half_life_y == 0:
            while self.half_life_y == 0 and self.attempt < 3:
                self.attempt += 1 
                self.step = self.step/100
 
                self.ranging = np.arange(self.min, self.max, self.step) # normally it 0.001, but the slope is so radical, can't catxh half-life.
                for j in np.nditer(self.ranging):

                    self.new_x = np.array([j])

                    #self.h = self.out.eval_components(self.out.params,x=self.new_x)
                    #self.hl_eval = list(self.h.values())[-1]

                    self.hl_eval = self.out.eval(self.out.params,x=self.new_x)

                    if  self.hl_eval >= 0.50 and self.hl_eval <= 0.51:

                        self.hl_array = np.append(self.hl_array,self.hl_eval)
                        self.hl_coord = np.append(self.hl_coord,self.new_x)

                self.half_life_id = np.argmin(np.abs(self.hl_array - 0.5))
                self.half_life_y = self.hl_array[self.half_life_id]
                self.half_life_x = self.hl_coord[self.half_life_id]
                self.peak_flag = self.peak_flag

            if self.half_life_y == 0:
                self.peak_flag = "Above permitted interpolation iterations"

        return self.half_life_y, self.half_life_x, self.peak_flag

    def saving(self):
            with  open('model_fit_c4_average_filtering_compound.txt', 'a') as f:

                    f.write("%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n" 
                        % (self.gene_id, self.cluster, self.model_flag, self.chisq, self.half_life_y, self.half_life_x,
                        self.norm_vals[0], self.norm_vals[1], self.norm_vals[2], self.norm_vals[3], self.norm_vals[4], self.norm_vals[5],
                        self.bestfit[0], self.bestfit[1], self.bestfit[2], self.bestfit[3], self.bestfit[4], self.bestfit[5], self.peak_flag))
 
try:
    os.remove("model_fit_c4_average_filtering_compound.txt")
except:
    pass

with  open('model_fit_c4_average_filtering_compound.txt', 'a') as f:

    (f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
        % ("Gene ID","cluster","model","chisq","half_life_y","self.half_life_x",
        "norm_val_0","norm_val_1","norm_val_2","norm_val_6","norm_val_12","norm_val_24",
        "bestfit_0","bestfit_1","bestfit_2","bestfit_6","bestfit_12","bestfit_24","Notes")))

cluster = 4
referencefile = "expression_norm1_vals_1.25_average_filtering.txt"
reference = []
with open(referencefile) as inputfile:
    for line in inputfile:
        reference.append(line.strip().split('\t'))


idsfile = "average_chron_cluster4_clean.txt"
#idsfile = "test2.txt"
idslist = []
with open(idsfile) as inputfile:
    for line in inputfile:
        idslist.append(line.strip().split('\t'))



ids=[item[0] for item in idslist]

print len(ids)
print "-----------"
print cluster

for i in ids:
   gene_id = gene_set(i, cluster)
   #x = np.sort(np.append(gene_id.x,gene_id.half_life_x))
   #x = gene_id.x
   #y = np.array(gene_id.bestfit)
   #y = np.array(gene_id.bestfit)
   #plt.plot(x, gene_id.exp_out.init_fit, 'b--')
   #plt.plot(x, gene_id.exp_out.best_fit, 'b-')
   #plt.plot(x,y,'ko-')
   #y = np.sort(np.append(gene_id.bestfit, gene_id.half_life_y))[::-1]
   #plt.plot(x,y,label= i)
   #plt.legend(loc='best')
   #print gene_id.__dict__
   #print i, gene_id.idx, gene_id.min, gene_id.max, gene_id.peak_flag, gene_id.bestfit, gene_id.half_life_y, gene_id.half_life_x
   print i, gene_id.half_life_y, gene_id.half_life_x, gene_id.peak_flag, gene_id.step
   #plt.vlines(gene_id.half_life_x, 0, 1, colors='g', linestyles='dashed')
   #print gene_id.bestfit


#plt.hlines(0.5, 0, 25, colors='r', linestyles='dashed')
#plt.show()   

