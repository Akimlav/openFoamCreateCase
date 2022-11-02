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
import fileinput
import sys
import initialConditions as ic

data1 = np.genfromtxt('inletData',skip_header=1,invalid_raise = False)
mesh = np.genfromtxt('meshData',skip_header=1,invalid_raise = False)

case_location = os.getcwd()
air_data = np.genfromtxt(case_location + '/meanProp/air',skip_header=1,invalid_raise = False)
fresh_data = np.genfromtxt(case_location + '/meanProp/fresh_water',skip_header=1,invalid_raise = False)
sea_data = np.genfromtxt(case_location + '/meanProp/sea_water',skip_header=1,invalid_raise = False)

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

def air (T):
    air_rho = np.interp(T,air_data[:,][:,0],air_data[:,][:,1])
    air_nu = np.interp(T,air_data[:,][:,0],air_data[:,][:,2])
    return (air_rho, air_nu)

def sea (T):
    sea_rho = np.interp(T, sea_data[:,][:,0], sea_data[:,][:,1])
    sea_nu = np.interp(T, sea_data[:,][:,0], sea_data[:,][:,2])
    return (sea_rho, sea_nu)

def fresh (T):
    fresh_rho = np.interp(T, fresh_data[:,][:,0], fresh_data[:,][:,1])
    fresh_nu = np.interp(T,fresh_data[:,][:,0], fresh_data[:,][:,2])
    return (fresh_rho, fresh_nu)

T = ic.T
nutnu = ic.nutnu
Tu = ic.Tu
endT = ic.endTime
writeInt = ic.writeInterval
n_of_proc = ic.n_of_proc

water_prop = fresh(T)
air_prop = air(T)

#turbulence calculator function
def turb_calc( U, nu, nutnu, Tu ):
    Cnu = 0.09
    u1 = (Tu/100)*U
    k = 3/2*u1**2
    eps = (k**2*Cnu)/(nutnu*nu)
    o = eps/(k*Cnu)
    return (k, eps, o)

def createCase (i, case_path, mesh_path):
    src = 'headCase'
    dst = case_path
    shutil.copytree(src, dst)
    src = 'headMesh/constant/polyMesh'
    dst = case_path + '/constant/polyMesh'
    shutil.copytree(src, dst)
    src = case_path +'/0.orig'
    dst = case_path +'/0'
    shutil.copytree(src, dst)
    turb_parameters = (turb_calc(abs(data1[i][0]), water_prop[1], nutnu, Tu))
    path1 = case_path +'/0/include/initialConditions'
    path2 = case_path +'/constant/transportProperties'
    path3 = case_path +'/constant/dynamicMeshDict'
    path4 = case_path +'/system/controlDict'
    path7 = case_path +'/system/decomposeParDict'
    f1 = open(path1, 'r')
    f2 = open(path1+'1', 'w')
    checkWords = ("U*","k*","omega*")
    repWords = (str(data1[i][0]),str(turb_parameters[0]),str(turb_parameters[2]))
    for line in f1:
        for check, rep in zip(checkWords, repWords):
            line = line.replace(check, rep)
        f2.write(line)
    f1.close()
    f2.close()
    os.remove(path1)
    os.rename(path1+'1',path1)
    f3 = open(path2, 'r')
    f4 = open(path2+'1', 'w')
    checkWords2 = ("rho_air*","nu_air*", "rho_water*","nu_water*")
    repWords2 = (str(air_prop[0]),str(air_prop[1]), str(water_prop[0]), str(water_prop[1]))
    for line in f3:
        for check2, rep2 in zip(checkWords2, repWords2):
            line = line.replace(check2, rep2)
        f4.write(line)
    f3.close()
    f4.close()
    os.remove(path2)
    os.rename(path2+'1',path2)
   # replaceAll(path3, 'n*', str(data1[i][1]))
    replaceAll(path4, 'endTime*', str(endT))
    replaceAll(path4, 'writeInterval*', str(writeInt))
    replaceAll(path4, 'rhoInf*', str(water_prop[0]))
    replaceAll(path7, 'n_proc*', str(n_of_proc))
    path5 = case_location + '/'+case_path
    path6 = case_location + '/'+case_path + '/run.sh'
    replaceAll(path6, "loc1", case_path)
    replaceAll(path6, "n_of_proc*", str(n_of_proc))
    with fileinput.FileInput(path5 + '/run.sh' , inplace=True) as file:
        for line in file:
            print(line.replace('loc', path5 ), end='')
    return()

def createMesh (i, mesh_path):
    src = 'headMesh'
    dst = mesh_path
    shutil.copytree(src, dst)
    path1 = mesh_path +'/system/blockMeshDict'
    replaceAll(path1, 'cells*', str(int(mesh[i][1])))
    
def prop_forces(path):
    #replace data    
    f1 = open(path + '/postProcessing/forces/0/force.dat', 'r') 
    f2 = open(path + '/postProcessing/forces/0/force_cleared.dat', 'w') 
    f3 = open(path + '/postProcessing/forces/0/moment.dat', 'r') 
    f4 = open(path + '/postProcessing/forces/0/moment_cleared.dat', 'w') 
    for line in f1: 
        f2.write(line.replace('(', ' ').replace(')', ' '))    
    f1.close() 
    f2.close() 
    for line in f3: 
        f4.write(line.replace('(', ' ').replace(')', ' '))    
    f3.close() 
    f4.close() 
    path = path
    force = np.genfromtxt(path + '/postProcessing/forces/0/force_cleared.dat',invalid_raise = False)
    moment = np.genfromtxt(path + '/postProcessing/forces/0/moment_cleared.dat',invalid_raise = False) 
    split = np.array_split(force,4) 
    force = split[3] 
    split = np.array_split(moment,4) 
    moment= split[3]
    time = force[:,][:,0] 
    time2 = moment[:,][:,0] 
    F_x = force[:,][:,1] 
    F_p = force[:,][:,4] 
    F_v = force[:,][:,7] 
    M_x = moment[:,][:,1] 
    M_p = moment[:,][:,4] 
    M_v = moment[:,][:,7]
    avg_F = np.mean(F_x) 
    avg_Fp = np.mean(F_p) 
    avg_Fv = np.mean(F_v) 
    avg_M = np.mean(M_x) 
    avg_Mp = np.mean(M_p) 
    avg_Mv = np.mean(M_v)
    res = [avg_F, avg_Fp, avg_Fv, avg_M, avg_Mp, avg_Mv, ]
    return(res)

def yPlus(path):
    path = path
    yPlus = np.genfromtxt(path + '/postProcessing/yPlus/yPlus.dat',invalid_raise = False)    
    yPlus_1 = yPlus[::2]
    n_str = sum(1 for line in yPlus_1) -1
    y_min = yPlus_1[n_str,:][2]
    y_max = yPlus_1[n_str,:][3]
    y_avg = yPlus_1[n_str,:][4]
    y_data = [y_min, y_max, y_avg]
    return(y_data)

