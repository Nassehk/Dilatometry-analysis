# -*- coding: utf-8 -*-
"""
Created on Thu May 14 08:52:35 2015

@author: nasseh
"""
import csv
#im1port sys
import numpy as np
#import glob
import matplotlib.pyplot as plt
import os
import pickle
import data_set_conditioner_V2

current_dir=os.getcwd()
file_list=pickle.load(open(current_dir+'/working_directory/file_list_conditioned.pkl','rb'))
#file_list=glob.glob('*.csv')
#d=open('master_setup.txt','r')
#for line in d:
#    temp=line.split(',')
#    temp[1:]=map(float,temp[1:])
#    file_list.append(temp)

#time_master=[]
#temp_master=[]
#dil_master=[]
#k=0

#for var in file_list:
#    
##    if sys.platform[0:5]=='linux':
##        try:
##            Data=open('/media/nasseh/FE0E012C0E00E00F/Science_Backup/Results/Python_Dil_analysis/Fitting_real_function/'+var[0],'rb')
##        except:
##            Data=open('/home/nasseh/Desktop/Python_Dil_analysis/Fitting_real_function/'+var[0],'rb')
##    else:
##        Data =  open('G:\\Science_Backup\\Results\\Python_Dil_analysis\\Fitting_real_function\\'+var[0], 'rb')
#
#    Data=open(var,'r')
#    D= csv.reader(Data)
#    t = []
#    temp = []
#    dil = []
#    for row in D:
#        t.append(row[0])
#        temp.append(row[1])
#        dil.append(row[2])
#    del t[0],temp [0] ,dil [0]
#    Data.close()
#    t = np.array(map(float,t))
#    temp = np.array(map(float,temp))
#    dil =np.array(map(float,dil))
#    time_master.append(t)
#    temp_master.append(temp)
#    dil_master.append(dil)
#    del t, temp, dil, row

current_dir=os.getcwd()
time_master=pickle.load(open(current_dir+'/working_directory/time_conditioned.pkl','rb'))
temp_master=pickle.load(open(current_dir+'/working_directory/temp_conditioned.pkl','rb'))    
dil_master=pickle.load(open(current_dir+'/working_directory/dil_conditioned.pkl','rb'))

for k in range(len(file_list)):
    plt.figure()
    plt.plot(temp_master[k],dil_master[k],'o')
    plt.title(file_list[k][0])
    plt.xlabel('Temperature (C)')
    plt.ylabel('Dilation (micron)')
    xy=zip(temp_master[k],dil_master[k])
    for i in xrange(200):    
        plt.annotate(i*len(xy)//200,(temp_master[k][i*len(xy)//200],dil_master[k][i*len(xy)//200]))
    