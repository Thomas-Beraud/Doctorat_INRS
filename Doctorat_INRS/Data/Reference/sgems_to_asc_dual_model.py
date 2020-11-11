# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:22:28 2020

@author: tbera

To import GSLIB grid properties exported from SGEMS, add the suffix "back" after the back_transform
"""


import numpy as np
import pandas as pd

f = open("data_permea_ensemble_low")
f1 = open("data_permea_ensemble_high")


f.readline()
f.readline()
f1.readline()
f1.readline()

index = np.zeros(90,dtype=int)
index1 = np.zeros(90,dtype=int)

j=0
k = 0
for i in range(180): # 180 is the number of simulation in the GSLIB file
	if 'gauss' in f.readline():
		index[j] = i
		j+=1
	if 'gauss' in f1.readline():
		index1[k] = i
		k+=1


final_data = np.loadtxt("data_permea_ensemble_low",skiprows=182).T  #delimiter=',',
final_data_1 = np.loadtxt("data_permea_ensemble_high",skiprows=182).T  #delimiter=',',

path = "k_fields/"

i = 0
for col in index :
	file = open(path+"fgen92_sim"+str(i)+".asc",'w')
	file.write("          realization #    "+str(i)+"       \n")
	file.write("  41    41   21     2.5     2.5     1    \n")

	for lines in range(len(final_data[col])):
		if lines < 18491 or lines > 25215 :
			file.write(str(final_data_1[col][lines])+"\n")
		else :
			file.write(str(final_data[col][lines])+"\n")

	file.close()
	i+=1