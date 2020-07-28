# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 03:55:32 2015


@author: nasseh
"""

# This program is called data_set_conditioner. It gets rid noise and excessive data point.
# This version is base on ata_set_conditioner_V4 in the old version control system
import csv
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import polynomial as P
import pickle
import os
from os.path import exists, join

#Load setup file which contains information that master_fitter and fitter needs to work.
#each line in the set file is like this:
#filename,Aus_min_row,Aus_max_row,fer_min_row,fer_max_row,L0,cementite conditions,
#values must be separated by comma "," and all lines must end with a comma.
make_raw_data_plots = False
raw_data_dir = 'raw_dil_data'
file_list=[]
d=open('master_setup.txt','r')
for line in d:
    temp=line.split(',')
    temp[1:6]=map(float,temp[1:6]) # get rid of end character and turn numbers into float
    if temp[0][0]=='#' or temp[0][0]=='%': #this if statement ignors the lines in master_setup file that begin with '#' or '%'.
        pass
    else:
        file_list.append(temp)

current_dir=os.getcwd()
#A=os.listdir(current_dir)
if not exists('working_directory'):        
    os.makedirs('working_directory')

output = open(current_dir+'/working_directory/file_list_conditioned.pkl', 'wb')    
pickle.dump(file_list, output)      
output.close()  

time_master=[]
temp_master=[]
dil_master=[]
k=0

for var in file_list:
#    if sys.platform[0:5]=='linux':
#        try:
#            Data=open('/media/nasseh/FE0E012C0E00E00F/Science_Backup/Results/Python_Dil_analysis/Fitting_real_function/'+var[0],'rb')
#        except:
#            Data=open('/home/nasseh/Desktop/Python_Dil_analysis/Fitting_real_function/'+var[0],'rb')
#    else:
#        Data =  open('G:\\Science_Backup\\Results\\Python_Dil_analysis\\Fitting_real_function\\'+var[0], 'rb')
    
    Data=open(join(raw_data_dir,var[0]),'r')
    D= csv.reader(Data)
    t = []
    temp = []
    dil = []
    for row in D:
        t.append(row[0])
        temp.append(row[1])
        dil.append(row[2])
    del t[0],temp [0] ,dil [0]
    Data.close()

    t = np.array(t, dtype = 'float')
    temp = np.array(temp, dtype = 'float')
    dil =np.array(dil, dtype = 'float')*1E-6
    time_master.append(t)
    temp_master.append(temp)
    dil_master.append(dil)
    del t, temp, dil, row

reg_order=2 #The order of polynomial used in noise filtering. Higher sr better have higher reg-order
def regresser(bottom,center,top):
    time_reg=time_fit[center]
    coef1= P.polyfit(time_fit[bottom:top],temp_fit[bottom:top],reg_order)
    temp_reg=P.polyval(time_fit[center],coef1)
    coef1= P.polyfit(time_fit[bottom:top],dil_fit[bottom:top],reg_order)
    dil_reg=P.polyval(time_fit[center],coef1)
    regressed_point=[time_reg,temp_reg,dil_reg]
    return regressed_point

#def range_finder(T,dT)

start=0
time_cond=[] #conditioned time data
temp_cond=[]
dil_cond=[]

time_cond_master=[] #master file that has all of the conditioned time data
temp_cond_master=[]
dil_cond_master=[]

regressed_point_list=[]

for i in range(len(time_master)):
    print (i)
    time_cond.append(time_master[i][0])
    temp_cond.append(temp_master[i][0])
    dil_cond.append(dil_master[i][0])
    time_fit=time_master[i]
    dil_fit=dil_master[i]
    temp_fit=temp_master[i]
    points=[]
    for j in range(1,len(time_master[i])):
        try:        
            reged_point=regresser(j-20,j,j+20)
        except:
#            print"out"
            reged_point=[time_fit[j],temp_fit[j],dil_fit[j]]
        if abs(time_fit[j]-time_cond[-1])>=10:
            points.append(j)
            time_cond.append(time_fit[j])          
        elif abs(reged_point[1]-temp_cond[-1])>=3:
            points.append(j)
            temp_cond.append(temp_fit[j])
    regressed_point_list.append(points)   
time_cond=[]
temp_cond=[]
dil_cond=[]         
for i in range(len(time_master)):
    time_cond=[]
    temp_cond=[]
    dil_cond=[]    
    time_cond.append(time_master[i][0])
    temp_cond.append(temp_master[i][0])
    dil_cond.append(dil_master[i][0])
    time_fit=time_master[i]
    dil_fit=dil_master[i]
    temp_fit=temp_master[i]
    points=regressed_point_list[i]
    for j in range(1,len(points)-1):
        regressed=regresser(points[j-1],points[j],points[j+1])
        time_cond.append(regressed[0])
        temp_cond.append(regressed[1])
        dil_cond.append(regressed[2])
    time_cond_master.append(time_cond) 
    temp_cond_master.append(temp_cond)
    dil_cond_master.append(dil_cond)


output = open(current_dir+'/working_directory/time_conditioned.pkl', 'wb')    
pickle.dump(time_cond_master, output)      
output.close()        

output = open(current_dir+'/working_directory/temp_conditioned.pkl', 'wb')    
pickle.dump(temp_cond_master, output)      
output.close()

output = open(current_dir+'/working_directory/dil_conditioned.pkl', 'wb')    
pickle.dump(dil_cond_master, output)      
output.close()  
    
if make_raw_data_plots:    
    plt.figure(figsize=(8,6))
    for i in range(len(file_list)):
        plt.plot(time_master[i],1e6*dil_master[i], label="Raw "+file_list[i][0])
        
        
    for i in range(len(file_list)):
        plt.plot(time_cond_master[i],1e6*np.array(dil_cond_master[i]),".", label="Conditioned "+ file_list[i][0])
    plt.legend()
    plt.title("Original data vs conditioned")
    plt.xlabel("Time $(s)$")
    plt.ylabel("Dilation $(\mu m)$")
    
    plt.figure(figsize=(8,6))
    
    plt.title("Original data vs conditioned")
    plt.xlabel("Temperature $(\degree C)$")
    plt.ylabel("Dilation $(\mu m)$")
    
    for i in range(len(file_list)):
        plt.plot(temp_master[i],1e6*dil_master[i], label="Raw "+file_list[i][0])
        
    for i in range(len(file_list)):
        plt.plot(temp_cond_master[i],1e6*np.array(dil_cond_master[i]),".", label="Conditioned "+file_list[i][0])
    plt.legend()
    
    plt.figure(figsize=(8,6))
    plt.title("Original data vs conditioned")
    plt.xlabel("Time $(s)$")
    plt.ylabel("Temperature $(\degree C)$")
    
    for i in range(len(file_list)):
        plt.plot(time_master[i],temp_master[i],label=file_list[i][0])
        
    for i in range(len(file_list)):
        plt.plot(time_cond_master[i],temp_cond_master[i],"o", label=file_list[i][0])
    plt.legend()
    #plt.show()

