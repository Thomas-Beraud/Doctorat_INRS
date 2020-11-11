import shutil
import os
import numpy as np
from subprocess import call
import multiprocessing as mp
    
def call_saltflow(sim):


    tmp_saltflow = 'temp_saltflow_'+str(sim)

    newpath = path+'/temp_saltflow/'+tmp_saltflow
    if not os.path.exists(newpath):
        os.makedirs(newpath)


    shutil.copy(path+folder+'/fgen92_sim'+str(sim)+'.asc',path+'temp_saltflow/'+tmp_saltflow+'/fgen92.asc')

    shutil.copy(path+'temp_saltflow/salt.data',path+'temp_saltflow/'+tmp_saltflow+'/salt.data')
    shutil.copy(path+'temp_saltflow/saltflow.exe',path+'temp_saltflow/'+tmp_saltflow+'/saltflow.exe')

    os.chdir(path+'temp_saltflow/'+tmp_saltflow)

    filelist = [ f for f in os.listdir(path+'temp_saltflow/'+tmp_saltflow)]


    while len(filelist) < 3:

        shutil.copy(path+folder+'/fgen92_sim'+str(sim)+'.asc',path+'temp_saltflow/'+tmp_saltflow+'/fgen92.asc')

        shutil.copy(path+'temp_saltflow/salt.data',path+'temp_saltflow/'+tmp_saltflow+'/salt.data')
        shutil.copy(path+'temp_saltflow/saltflow.exe',path+'temp_saltflow/'+tmp_saltflow+'/saltflow.exe')

        filelist = [ f for f in os.listdir(path+'temp_saltflow/'+tmp_saltflow)]

    cmd = path+'temp_saltflow/'+tmp_saltflow+'/saltflow.exe'

    print(cmd)
    call(cmd)

    shutil.move(path+'temp_saltflow/'+tmp_saltflow+'/break.out',path+break_+'/break_'+str(sim)+'.out')
    shutil.rmtree(path+'temp_saltflow/'+tmp_saltflow)
    os.rmdir(path+'temp_saltflow/'+tmp_saltflow)



path = "D:/Dropbox/Doctorat_INRS-master/Doctorat_INRS/Doctorat_INRS/"


f =  open(path+"temp_saltflow/folder.txt", "r")

folder = f.readline()[:-1]
break_ = f.readline()
f.close()


if __name__ == '__main__':
    

       
    
    pool = mp.Pool(90)
    
    
    sim = np.arange(0,15,1)
    # Make the Pool of workers
    
    result = pool.imap_unordered(call_saltflow, sim)
    
    # Close the pool and wait for the work to finish
    pool.close()
    pool.join()

