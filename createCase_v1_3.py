## -*- coding: utf-8 -*-
#"""
#Created on Fri Oct  4 17:36:22 2019
#
#@author: LavrinenkoAV
#"""
#CreateCase version 1.2

import os
import numpy as np
import shutil
import subprocess
import datetime
import time
from openFoam_func_lib import *

start = time.time()
now = datetime.datetime.now()
print('Running python script')
print( "Start time and date: "+ now.strftime("%H:%M %d-%m-%Y "))

data1 = np.genfromtxt('inletData',skip_header=1,invalid_raise = False)
mesh = np.genfromtxt('meshData',skip_header=1,invalid_raise = False)

n = (sum(1 for line in data1))    # -1 пропускаем первую строку
m = (sum(1 for line in mesh))
print(n, m)
#create folder
for j in range (0, m):
    mesh_path = 'mesh1_'+ str(mesh[j][0])+'_c=' + str(mesh[j][1])
    meshpath = os.path.join(mesh_path)
    meshpath1 = os.path.join(mesh_path+'/constant/polyMesh')
    if os.path.exists(meshpath):
        print('old mesh are!')
        if os.path.exists(meshpath1):
            print('mesh are!'+ mesh_path)
        else:
            shutil.rmtree(meshpath)
            print('Old mesh deleted')
            createMesh(j, mesh_path)
            subprocess.call(case_location + '/'+mesh_path +'/mesh')# 
    else:
        createMesh(j, mesh_path)
        subprocess.call(case_location + '/'+mesh_path +'/mesh')

    for i in range (0, n):
        start_cycle = time.time()
        print('*************************')
        case_path = 'case1_U=' + str(data1[i][0])
        dirpath = os.path.join(case_path)
        dirpath1 = os.path.join(case_path+'/postProcessing')
        if os.path.exists(dirpath):
            print('old Case are!')
            if os.path.exists(dirpath1):
                print('post Case are!'+ case_path)
            else:
                shutil.rmtree(dirpath)
                print('Old case_U='+str(data1[i][0])+'_n='+str(data1[i][1]) + ' deleted')

                createCase(i, case_path, mesh_path)

                print('Running bash script')
                print('*************************')

                subprocess.call(case_location + '/'+case_path +'/run.sh')

                print('Bash script done')
                print('*************************')
                cycle_time = (time.time()-start_cycle)/3600
                print('It took', cycle_time, 'hours.')
                print('*************************')
        else:
            createCase(i, case_path, mesh_path)
            print('Running bash script')
            print('*************************')

            subprocess.call(case_location + '/'+case_path +'/run.sh')

            print('Bash script done')
            print('*************************')
            cycle_time = (time.time()-start_cycle)/3600
            print('It took', cycle_time, 'hours.')
            print('*************************')
