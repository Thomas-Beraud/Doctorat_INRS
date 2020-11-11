# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 14:32:03 2020

@author: tbera
"""

import numpy as np
import matplotlib.pyplot as plt

global_permea = np.loadtxt("FGEN_global/fgen92.asc",skiprows=2)
layered_permea = np.loadtxt("FGEN_layer/fgen92.asc",skiprows=2)

nx = 41
ny = 41

section_2d = nx*ny

# Depth of layer and depth of layer :

index_begin = section_2d*11
index_end = section_2d*15

global_permea[index_begin:index_end] = layered_permea[index_begin:index_end]

section = global_permea.reshape(21,41,41)


plt.imshow(section[12,:,:],origin='lower')
plt.colorbar()
plt.show()

plt.imshow(section[:,10,:],origin='lower')
plt.colorbar()
plt.show()

plt.imshow(section[:,:,10],origin='lower')
plt.colorbar()
plt.show()


file = open("fgen92.asc",'w')

file.write('          realization #     1\n   41   41   21     2.500     2.500     1.000\n')

for line in global_permea :
	file.write(str(line)+'\n')
file.close()



file = open("wells_permea_high.csv",'w')
file1  = open("wells_permea_low.csv",'w')

file.write('x,y,z,ln_permea\n')
file1.write('x,y,z,ln_permea\n')

wells = [[22,18],[27,37],[25,32],[38,6],[10,40]]

pumping = [[25,25],[15,10],[40,35]]


for well in wells:
	for z in range(21):
		if z < 11 or z >14 :
			file.write(str(well[0]*2)+','+str(well[1]*2)+','+str(z)+','+str(global_permea[z*nx*ny+well[1]*nx+well[0]])+'\n')
		else :
			file1.write(str(well[0]*2)+','+str(well[1]*2)+','+str(z)+','+str(global_permea[z*nx*ny+well[1]*nx+well[0]])+'\n')

file.close()
file1.close()

plt.imshow(section[12,:,:],origin='lower')
for i in range(len(wells)):
	plt.scatter(wells[i][0],wells[i][1],color='red')
for i in range(len(pumping)):
	plt.scatter(pumping[i][0],pumping[i][1],color='blue')
plt.show()

plt.imshow(section[:,10,:],origin='lower')
plt.show()

plt.imshow(section[:,:,10],origin='lower')
plt.show()
