import numpy as np
import re
import Functions.Import_Data as ID
import matplotlib.pyplot as plt


def load_results(assim_number, path, bool_plot):

    model = 90
    number_of_well = 5
    wells = ["Well O"+str(i)+"" for i in range(number_of_well)]
    wells[3] = "Well C3"
    wells[4] = "Well C4"
    coordinates = [[55,75,8],
                   [75,25,8],
                   [25,22.5,13],
                   [20,80,7],
                   [55,15,11]]


    if bool_plot :
        plt.figure(figsize=(6,6))
        for i in range(number_of_well):
            if i < 3 :
                plt.scatter(coordinates[i][0],coordinates[i][1],c=['lightsteelblue'])
            else :
                plt.scatter(coordinates[i][0],coordinates[i][1],c=['navy'])

            plt.text(coordinates[i][0]*1.1,coordinates[i][1]*1.1,wells[i],fontsize=14)
            plt.xlim(0,100)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.ylim(0,100)

        plt.scatter(50,50,c=['red'],label='Pumping Well')
        plt.scatter(30,20,c=['red'])
        plt.scatter(80,70,c=['red'])
        plt.legend()
        plt.show()



    brk_ref  = open(path+'Data/Reference/break.out', 'r')

    time = []
    cpk = []
    # skipping the 2 header lines
    brk_ref.readline()
    brk_ref.readline()

    for line in brk_ref:

        words = re.split(',|  | ',line)
        while '' in words:
            del words[words.index('')]

        if words[0]=='zone':
            time_step = int(words[2])

            time_temp = np.zeros(time_step)
            cpk_temp = np.zeros(time_step)
            iter_ = 0
            time.append(time_temp)
            cpk.append(cpk_temp)

        if words[0]!='zone' and iter_ < time_step:
            time_temp[iter_] = float(words[0])
            cpk_temp[iter_] = float(words[2])
            iter_+= 1


    time = np.array(time)
    cpk = np.array(cpk)


    # Writing the output file

    f = open("Data/measures.txt",'w')

    f.write(str(number_of_well)+" ; number of observation point\n")

    for i in range(number_of_well):
        f.write(str(wells[i])+" ; name of observation point\n")
        f.write(str(len(cpk[i]))+" ; number of measures\n")
        f.write(str(coordinates[i][0])+" ; x coordinate of observation point\n")
        f.write(str(coordinates[i][1])+" ; y coordinate of observation point\n")
        f.write(str(coordinates[i][2])+" ; z coordinate of observation point\n")
        for j in range(len(time[i])):
            f.write(str(time[i][j])+","+str(cpk[i][j])+"\n")
    f.close()


    brk_obs  = open(path+'Data/break_update_'+str(assim_number)+'/break_'+str(0)+'.out', 'r')

    iter_ = 0
    time = []
    cpk_global = []
    cpk = []
    # skipping the 2 header lines
    brk_obs.readline()
    brk_obs.readline()

    for line in brk_obs:

        words = re.split(',|  | ',line)
        while '' in words:
            del words[words.index('')]

        if words[0]=='zone':
            time_step = int(words[2])

            time_temp = np.zeros(time_step)
            cpk_temp = np.zeros(time_step)
            iter_ = 0
            time.append(time_temp)
            cpk.append(cpk_temp)

        if words[0]!='zone' and iter_ < time_step:
            time_temp[iter_] = float(words[0])
            cpk_temp[iter_] = float(words[2])
            iter_+= 1


    cpk_global.append(cpk)

    brk_obs.close()

    for k in range(model-1):

        brk_obs  = open(path+'Data/break_update_'+str(assim_number)+'/break_'+str(k+1)+'.out', 'r')
        iter_ = 0
        cpk = []
        # skipping the 2 header lines
        brk_obs.readline()
        brk_obs.readline()

        for line in brk_obs:

            words = re.split(',|  | ',line)
            while '' in words:
                del words[words.index('')]

            if words[0]=='zone':
                time_step = int(words[2])
                cpk_temp = np.zeros(time_step)
                iter_ = 0
                cpk.append(cpk_temp)

            if words[0]!='zone' and iter_ < time_step:
                cpk_temp[iter_] = float(words[2])
                iter_+= 1

        cpk_global.append(cpk)
        brk_obs.close()

    cpk_global = np.array(cpk_global)


    # Writing the output file

    f = open("Data/observed_"+str(assim_number)+".txt",'w')

    f.write(str(number_of_well)+" ; number of observation point\n")
    f.write(str(model)+" ; number of ensemble\n")

    for i in range(number_of_well):
        f.write(str(wells[i])+" ; name of observation point\n")
        f.write(str(len(time[i]))+" ; number of measures\n")
        f.write(str(coordinates[i][0])+" ; x coordinate of observation point\n")
        f.write(str(coordinates[i][1])+" ; y coordinate of observation point\n")
        f.write(str(coordinates[i][2])+" ; z coordinate of observation point\n")
        for j in range(len(time[i])):
            f.write(str(time[i][j])+",")
            for n in range(model):
                f.write(str(cpk_global[n][i][j])+",")
            f.write(str(cpk_global[n][i][j])+"\n")

    f.close()
