import numpy as np
import matplotlib.pyplot as plt
import re

def import_permeability_field(filename,ensemble_number, cells_3d_model) :
    """
    filename : folder path and name to the permeability field ensemble

    k_ensemble : shape ((ensemble,cells_3d_model))

    """

    # Loop to read *.asc file coming from Saltflow, final ensemble_matrix must be of
    # dimension parameters as row, ensemble as column
    #
    #
    # For 100 ensemble and 1000 parameters :
    #
    # ensemble_matrix :
    #
    #          100
    #      |   ...   |
    # 1000 |   ...   |
    #      |   ...   |
    #

    k_ensemble = np.zeros((ensemble_number,cells_3d_model))

    for simu in range(ensemble_number) :

        k  = open(filename+str(simu)+'.asc', 'r')
        iter_ = 0

        for line in k:
            if iter_ == 1 :
                words = re.split(',|  | ',line)
                while '' in words:
                    del words[words.index('')]

                i = int(words[0])
                j = int(words[2])
                break
            iter_ += 1

        iter_ = 0
        for line in k:
            k_ensemble[simu][iter_] = float(line)
            iter_ += 1

        k.close()


    return np.matrix(k_ensemble.T)


class Well():
    def __init__(self, name, size, x, y, z):

        self.name = name
        self.size = size
        self.time = np.zeros(size)
        self.data = np.zeros(size)
        self.x = x
        self.y = y
        self.z = z


    def set_data(self,data,time):
        if len(data) == self.size :
            self.data = data
        if len(time) == self.size :
            self.time = time
        else :
            print("Your value array size is different than the one defining this well \nYour data array should be of size : "+str(self.size)+" ")

    def set_data_index(self,data,time,index):
        if self.size >= index :
            self.data[index] = data
            self.time[index] = time
        else :
            print("Your index is greater than the array size. Max index : "+str(self.size)+" ")

    def show_well_description(self) :
        print("Name : "+str(self.name)+"")
        print("Size : "+str(self.size)+"")
        print("X : "+str(self.x)+"")
        print("Y : "+str(self.y)+"")
        print("Z : "+str(self.z)+"")

    def show_well_data(self) :
        print(self.time,self.data)

    def get_well_description(self) :
        return(self.name,self.size,self.x,self.y,self.z)

    def get_well_data(self) :
        return(self.time,self.data)

    def plot(self):
        plt.figure(figsize=(5,3))

        plt.plot(self.time,self.data, 'r', label='Measured',linewidth=3)

        plt.xlabel('Time (days)',fontsize=16)
        plt.ylabel('Concentration (g/l)',fontsize=16)


        plt.xlim(0,max(self.time))
        plt.ylim(0,1.1*max(self.data))
        plt.legend(fontsize=16)



        plt.title(self.name+" \n x:"+str(self.x)+" y:"+str(self.y)+" z:"+str(self.z),fontsize=20)
        #plt.axvline(x=time[250])
        #plt.axvline(x=time[375])

        plt.show()

class Measures():
    """
     Take the path to your measured files
     It will be only a long 1D vector of measured made on your terrain or model
     The simulated data will be loaded with another function
    """


    def __init__(self, path):
        f = open(path,'r')

        self.nbr_obs_points = int(f.readline().partition(';')[0])

        self.wells = [Well for i in range(self.nbr_obs_points)]
        for i in range(self.nbr_obs_points):
            name = f.readline().partition(';')[0]
            size = int(f.readline().partition(';')[0])
            x = float(f.readline().partition(';')[0])
            y = float(f.readline().partition(';')[0])
            z = float(f.readline().partition(';')[0])

            self.wells[i] = Well(name, size, x, y, z)

            temp_data = np.zeros(size)
            temp_time = np.zeros(size)

            for j in range(size):
                temp = np.fromstring(f.readline(), dtype=float, sep=',')
                temp_data[j] = temp[1]
                temp_time[j] = temp[0]
                #temp_data[j] = np.(f.readline())

            self.wells[i].set_data(temp_data,temp_time)

        f.close()


class Observation():
    """
     Take the path to your observation files
     It will be only a long 2D vector of observations made in your model
    """


    def __init__(self, path):
        f = open(path,'r')

        self.nbr_obs_points = int(f.readline().partition(';')[0])
        self.ensemble_size = int(f.readline().partition(';')[0])

        self.wells = [[Well for j in range(self.ensemble_size)] for i in range(self.nbr_obs_points)]

        for i in range(self.nbr_obs_points):
            name = f.readline().partition(';')[0]
            size = int(f.readline().partition(';')[0])
            x = float(f.readline().partition(';')[0])
            y = float(f.readline().partition(';')[0])
            z = float(f.readline().partition(';')[0])

            for k in range(size):
                temp_data = np.fromstring(f.readline(), float,self.ensemble_size+1, sep=',')
                # Index 0 contain the time stamp

                for j in range(self.ensemble_size):
                    if k ==0 :

                        self.wells[i][j] = Well(name, size, x, y, z)
                        self.wells[i][j].set_data_index(temp_data[j+1], temp_data[0],  k)
                    else :
                        self.wells[i][j].set_data_index(temp_data[j+1], temp_data[0],  k)

        f.close()
