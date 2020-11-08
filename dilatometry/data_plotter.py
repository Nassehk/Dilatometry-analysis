# -*- coding: utf-8 -*-
"""
Created on Thu May 14 08:52:35 2015

@author: nasseh
"""
import numpy as np
#import glob
import matplotlib.pyplot as plt
import os
import pickle
import data_set_conditioner

current_dir=os.getcwd()
file_list=pickle.load(open(current_dir+'/working_directory/file_list_conditioned.pkl','rb'))

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
    xy=list(zip(temp_master[k],dil_master[k]))
    for i in range(200):    
        plt.annotate(i*len(xy)//200,(temp_master[k][i*len(xy)//200],dil_master[k][i*len(xy)//200]))
    
