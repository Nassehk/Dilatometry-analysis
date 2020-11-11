"""
Change log:

Master_fitter _v3.py
The code is ported to python 3.
TODO: update code with object orientation for phases.

Master_fitter_v2.52.py
Volume fraction calculations are done and results are stored in the excel file.
Martensite calculations is updated. The results with regard to martensite are more accurate.
results_xlsx_saver.py is update to allow easier manipulation of the module. The code has shrunk also.

Master_fitter_V2.51.py
Bug fix in plotting Figures 7 and 8. Improved Mf_dic generation for faster run.


Master_fitter_V2.50.py
Definition of martensite is handled better. Before this version code allowed for precipitation of cementite during martensitic transformation.
Multiple bug fixes.Plots are named properly with better axis titles.
User has the option to not make the plots.
User can tell the software to automatically open the results file externally.

To Do: Unify V10 and V10_plot, Volume of the phases at room temp.


Master_fitter_v2.48.py
Speeeeeeed! The program uses Memoizing and hash table to speed up the calculation.
About 20% increase in speed is achieved this way. Functions like Ms, Bs, L ...
which are shared between V10 and V10_plot are defined outside the fitter function which reduces the risk of 
V10 and V10_plot being different. Vcement is imported as a C library for speed.
Dataset conditioner is improved to V3 which is slightly more consistent in picking 
analysis points. 

Master_fitter_v2.44.py
Carbon Solubility in austenite can be introduced to master fitter using an excle file.

this version can calculate correction factor due to small differences in length of the samples. This can be done in two ways. 
The first method is done by adding a correction factor so that CTE of individual sample becomes equal to the average CTE.
The second method is done my adding a correction factor to length so that maximum strain during heating becomes equal to the average of maximum strain for all the samples.
Correction for L0 can be shut off through run settings. 

master_fitter_v2.3.py 
Speed Speed Speed
Mater fitter does parallel processing (you are welcome!) Best case, have one core per data file + 1 more to handle OS calls. 
Additional core will not improve speed. Two core functions a0_alpha and a0_bainite 
are written in C which provides orders of magnitude speed up. 
master_fitter_utils is the C shared library of these two functions which needs 
to be compiled for every computer system.

This version can optimize the weight factor of C on CTE of martensite and the weight 
factor of C on expansion of martensite along with some more which can be found in the code. 
Please refer to optimized_params for the full list.

Selective optimization of parameters is implemented.

Results are now stored in the working directory in the .XLSX format. 
A module named Results_xlsx_saver creates the results file

Program has gotten smarter. It can handle irregular data sampling rates during dilatometry. 
data_set_conditioner code takes care of this. 
This program runs every time before master_fitter runs to make sure that the changes in master_setup are taken into consideration.

A program called data_plotter is written to help with creating the master_setup_file. 
This part is still a little tricky. Follow these steps:
    1.Create master_setup.txt file by just writing the file names. Anything written after the file name is ignored by data_plotter. Suggestion edit the a working master_setup.txt that works
    2.Run the data_plotter.py program.
    3.In the plots created by data_plotter, Row numbers are annotated in the plots
    4.Use the plots and row numbers to identify the correct ranges for pure austenite and product.
    5.Fill in the blanks. Enter pure austenite start row then enter "," then enter pure austenite range end then enter ",". Do the same for product range.
    6.The last entry in each row is what you want the program to do about cementite. Read the code for more info on this.
    7.You have created your master_setup file!

convergence:
    The smaller the interval the more accurate the results are. However the accuracy comes at the price of speed. 
        Speed exponentially goes up as the interval goes down.
    The smallest interval is 1. Going from high to low interval has to be done gradually. 
        Sudden decrease in interval can make the program unstable because the 
        fit parameters may be too far from the answer.

Important notice:
    Before you run the code for the first time you need to compile the C shared library.
    These functions are all gathered in master_fitter_util.so which is created 
    using Cython tool from the source code master_fitter_utils.pyx. compile by typing 
    "python setup_master_fitter_utils.py build_ext --inplace" in a terminal that is 
    opened in the master_fitter_utils.py directory.
    

Sources of error:
    1. Master_setup.txt is not written according to the expected layout.
    2. When you get an error that referes to Vgama_param, most likely the initial guess of N_t is too big or small.
    3. "Cannot import Valpha, Vcement..." you may have to recompile the cython 
        tools using "python master_fitter_utils.py build_ext --inplace" in a terminal
        that is opened in the master_fitter_utils.py directory.
    

Send me your questions and inputs.
email the author at: khodaie@ualberta.ca or nassehk@gmail.com.
This code is distributed under MIT license. Update, upgrade, fix and share.
"""
__author__= "Nasseh Khodaie"
__version__=3.00
#import cProfile 
#def main():
import pandas as pd
import copy
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import polynomial as P
import scipy as scipy
from scipy import stats #for some reason it is necessory for optimize to have stats imported.
import os
import sys as sys
import time
import pickle
from numpy import sqrt
#from numba import jit, float64, vectorize 
#if sys.platform[0:3] != 'lin':
#    import winsound
from Results_xlsx_save import saver
#import data_set_conditioner_V2
import data_set_conditioner
from openpyxl import load_workbook as lw
from scipy.interpolate import interp1d
import multiprocessing as mp
from subprocess import call
from master_fitter_utils import Valpha,Vbainite,VCement
from functools import lru_cache
start=time.time()
clear = lambda: os.system('clear'if sys.platform[:5]==('linux' or 'darwi') else 'cls')
clear()


# global variable are solver settings that are shared between V10 and V10_plot functions.
global sr, reg_order, interval, Bs_model, end_fit_WF, overal_fit_WF, err_end_slop_WF
global err_maximum_transformation_WF, Bs_master_dic, Ms_master_dic, MF_dic

################################
#solver settings

optimize='yes'# It tells the program to optimize or just plot using the parameters provided.yes / no
show_plots='no'
open_excel_result='no'
interval=8

end_fit_WF=7
overal_fit_WF=7
err_end_slop_WF=6
err_maximum_transformation_WF=2.0
maximize_fraction_transformed='yes' #'yes' and 'y' will trigger maximizing transformation cost factor. Becareful with use of this. It may create incorrect lattice parameters.
#use maximization with caution. this is a tool to get close to the final solution. However if your final solution relies on maximization to be on, your solution might not be correct. 
sr=7
reg_order=3 


correct_L0='no'
L0_correction_method=2  # 1 for CTE matching which calculates correction to make all CTEs equal to average CTE.
                        # 2 for strain matching which matches the strain at austenitization for all of the samples to average starin. 
                        # 3 for normalizing all data to a strain at user defined temperature during cooling.                        
                        # Method 2 is a better method for solid samples as temperature gradient during cooling can effect CTE differently depending on the cooling rate 
normalizing_temp=895    # Only used for method 3 of L0_correction.
use_avr_CTE_product=0   # 1 for yes 0 for no. Setting to zero will ignore user provided values and will calculate  

Bs_model=1 # can be 1 for imperical and 2 for t zero method (model 2 only works for Nasseh's x80)
Ms_model=3 # can be 1 for Andrews1, 2 for Andrews 2, 3 for Capdevilla, 4 for interaction model Jiajun Wang
SSF_method="general" # can be user or genral. general assumes  linear increase of SSF from LBs to Ms, from SSF=0 to SSF=1, User defined is polynomial. 
use_XRD_for_a0_alpha='no' #yes overrides user defined value for a0_alpha and calculates a value based of literature.

#select which parameters you want to be optimized. 0 for no and 1 for yes.
optimization_method = 1                        
#                        1 for Nelder-Mead and 2 for differential evolution
#                        NM method searches for the best answer aroud an 
#                        starting point. It is faster but prone to pointing to a 
#                        local minimum not global. DE is more robust. 
#                        It searches for the best answer in a defined domain called bounds. It is however, slower than NM.

                       
optimized_param={'a0_gama':     0, 
                 'a0_alpha':    1,
                 'CTE_alpha_a': 0,
                 'CTE_alpha_b': 1,
                 'CTE_alpha_c': 1,
                 'c_wf_for_cte':1, 
                 'c_wf_for_a0': 1}

#bounds are used by differential evolution only
optimized_param_bounds={'a0_gama':     (3.63e-10*0.98,3.63e-10*1.02), 
                        'a0_alpha':    (2.86e-10*0.98,2.98e-10*1.02),
                        'CTE_alpha_a': (1.33e-5*0.9,1.33e-5*1.1),
                        'CTE_alpha_b': 0,
                        'CTE_alpha_c': 1,
                        'c_wf_for_cte':(1.96e-4,8.96e-3), 
                        'c_wf_for_a0': (4.53e-11,9.53e-12)}
#
                       
                        
#compile_C_library=1
#################################
#current_dir=os.getcwd()
#if compile_C_library==1:
#    if sys.platform[:3]=="win":
#
#        return_code = subprocess.call("python %s build_ext --inplace" %(current_dir+'\\Cython\\setup_master_fitter_utils.py'), shell=True) 
#    else:
#        return_code = subprocess.call("python current_dir=os.getcwd()+/Cython/setup_master_fitter_utils.py build_ext --inplace", shell=True) 


################################
#starting values of the fundamental parameters of the program.
#params are used by Nelde-Mead only 


#interval=2
a0_gama =     3.630354505175200875243544068702e-10
a0_alpha=     2.862496784948004707531459963016e-10
CTE_alpha_a=  0.000000000000000000000000000000e+00
CTE_alpha_b=  5.259052569948893497158337768937e-09
CTE_alpha_c=  1.303580639054562249575831833770e-05
c_wf_for_cte= 2.591159702736662369566833508117e-04
c_wf_for_a0=  2.996310242996993135602471934982e-12


#################################
result_master_time=[]
result_master_temp=[]
result_master_dil=[]
result_master_fraction=[]
result_master_ID=[]
result_master_FD_fraction=[]
result_master_C_in_FD_fraction=[] #fraction of ferrite derivative formed
result_master_C_in_FD_WP=[]
result_master_cem_fraction=[] #fraction of cementite formed
result_master_lattice_param_gamma_during_trans=[]
result_master_C_in_gama=[]
result_master_Vol_f_transformed=[]
result_master_FD_vol=[]#calculates volume of the phase at 20C
result_master_cem_vol=[]#calculates volume of the phase at 20C
result_master_aus_vol=[]#calculates volume of the phase at 20C
run_parameter={}

Bs_master_dic={}
MF_dic={}
Ms_master_dic={}
C_in_alpha_master_dic={}

Bs_temp_dic={}
MF_temp_dic={}
Ms_temp_dic={}
C_in_alpha_temp_dic={}

molar_weight={'fe':56,'mn':54.6,'cr':51.99, 'ni':78.7, 'mo':95.96, 'si':28.08, 'c':12, 'nb':92.9, 'ti':47.867, 'al':26.98, 'cu':63.546}
error_in_calcs=[]
elements=[]
amount0_wp=[]

if sr<interval: sr=interval
    
with open('chemistry.csv','r') as c:
    for line in c:
        temp=line.split(',')
        elements.append(str(temp[0].lower()))
        amount0_wp.append(float(temp[1]))
chemistry0=dict(zip(elements,amount0_wp))
chemistry0_temp=dict(zip(elements,amount0_wp))# this is used in some calculation and altered during the process. Should not be used for anything. use chemistry0 instead.
#def molar_fraction(elements,amounts_wp,req

for elm in molar_weight.keys():
    if elm not in chemistry0.keys():
        chemistry0.update({elm:0})


CSLA = pd.read_excel('C_solubility_limit_austenite.xlsx')
X_CSLA = CSLA.iloc[:,0].values
T_CSLA = CSLA.iloc[:,1].values
Solubility=interp1d(T_CSLA,X_CSLA,kind='linear')
CSLA_der=scipy.misc.derivative(Solubility,T_CSLA[2:-2],dx=0.1,n=2)


BB=np.unravel_index(CSLA_der.argmax(),CSLA_der.shape)[0] #BB is the index of maximum 2nd derivative of CSLA. this point is use for deviding the solubility data into two set to make the polynomial fit possible
BB=BB+(len(T_CSLA)-BB)//5

#def Solubility_aus_cement(T): # In mole frtaction, calculates the solubility limit of C in austenite with respect to cementite
CSLA_high=P.polyfit(T_CSLA[:BB],X_CSLA[:BB],4)
CSLA_low=P.polyfit(T_CSLA[BB:],X_CSLA[BB:],3)
#@jit(float64(float64))
def Solubility_aus_cement(T):#This function is not vectorizable because of the if statements.
    if T<T_CSLA[BB]:
        if T<T_CSLA.min():
            CSLA=0
        else:
            CSLA=P.polyval(T,CSLA_low)
    else:
        CSLA=P.polyval(T,CSLA_high)
    if CSLA<0:
        CSLA=0
    return CSLA
#plt.figure(87)
#plt.plot(T_CSLA,X_CSLA,'o')
CSLA=np.zeros(len(T_CSLA))
for i in range(len(T_CSLA)):
    CSLA[i]=Solubility_aus_cement(T_CSLA[i])

plt.figure(figsize=(8,6))
plt.title('Solubility of C in austenite')
plt.xlabel('Temperature $(\\degree C)$')
plt.ylabel('Solubility limit (mole fraction)')
plt.legend(loc='upper right')
#plt.plot((T_CSLA[2:-2],scipy.misc.derivative(Solubility,T_CSLA[2:-2],dx=0.1,n=2),'o',label='2nd derivative'))
plt.plot(T_CSLA,Solubility(T_CSLA),'o',label='Imported solubility data')
plt.plot(T_CSLA,CSLA,'--',label='Fitted function')
plt.legend(loc='upper left')


def a0_alpha_calculated(): #this calculates lattice parmater of ferrite at 20C using an equation found in letrature. 
    a0_alpha_wf={'mn':0.00006, 'si':-0.00003, 'ni':-0.00007, 'cr':0.00005, 
                 'p':-0.0001, 'ti':-0.0031} 
    a0_alpha_XRD=0.28664 # at room temp with no alloying element.
    for element in a0_alpha_wf.keys():
        if element in chemistry0.keys():
            a0_alpha_XRD=a0_alpha_XRD+a0_alpha_wf[element]*chemistry0[element]
    return a0_alpha_XRD


if use_XRD_for_a0_alpha.lower()=='yes':
    a0_alpha=a0_alpha_calculated()

        
def molar_fraction(elements,amounts_wp,required):
#    elements is a list of all alloying elements (excluding fe).
#    amounts is the coresponding weight pecent of the element in the alloy.
#    required is the name of the element that its molar fraction is required.   
    fe_wp=100-sum(amounts_wp)
    elements.append('fe')
    amounts_wp.append(fe_wp)
    element_mol=np.zeros(len(elements))
    for i in range(len(elements)):
        element_mol[i]=amounts_wp[i]/molar_weight[elements[i]]
    
    mole_fractions=element_mol/sum(element_mol)
    matched=dict(zip(elements,mole_fractions))
    return matched[required]
C0=molar_fraction(elements,amount0_wp,'c')

def MF_to_WP_generator(chem):
    local_chem = copy.deepcopy(chem)
    C_wp=np.linspace(chemistry0['c'],6.0,200)
    C_MF=[]
    for i in range(len(C_wp)):
        if 'fe' in local_chem.keys(): 
            local_chem.pop('fe')
        local_chem['c']=C_wp[i]
        elem,amo_wp=list(local_chem.keys()),list(local_chem.values())
        C_MF.append(molar_fraction(elem,amo_wp,'c'))
    MF_to_WP=P.polyfit(C_MF,C_wp,4)
    return MF_to_WP
MF_to_WP=MF_to_WP_generator(chemistry0)
            
def mf_to_wp(MF):
    try:#This will allow handelling numpy arrays
#        print "MF=", MF
        strMF=str(MF)[:(precision+2)] #precision of 0.0001. 
        try:#This will try finding strMF in the MF_dict.
            result=MF_dic[strMF] #this dictionary is created based on the results of all the other runs during optimization.
#            print "not bad"
        except:
            result=P.polyval(MF,MF_to_WP)
            MF_temp_dic[strMF]=result
#            print "updated Mf_temp_dic"                
        return result
#    else:
    except:
        return P.polyval(MF,MF_to_WP)

   
def CTE_alpha(T): # T in celsius
    return CTE_alpha_a*T*T+CTE_alpha_b*T+CTE_alpha_c


current_dir=os.getcwd()
file_list=pickle.load(open(current_dir+'/working_directory/file_list_conditioned.pkl','rb'))
time_master=pickle.load(open(current_dir+'/working_directory/time_conditioned.pkl','rb'))
temp_master=pickle.load(open(current_dir+'/working_directory/temp_conditioned.pkl','rb'))    
dil_master=pickle.load(open(current_dir+'/working_directory/dil_conditioned.pkl','rb'))
delta_L_tot=0
CTE_tot=0
strain_at_normalizing_temp_master=[]
aus_line_master=[]
CTE_product_master=[]
CTE_aus_master=[]

# This calculate correction factors for L0s
for i in range(len(file_list)):
    delta_L_tot=delta_L_tot+np.amax(dil_master[i])/file_list[i][5]
    row_min_aus=int(file_list[i][1])
    row_max_aus=int(file_list[i][2])
#    time_analysis_aus=time[row_min_aus:row_max_aus+1]
    dil_analysis_aus=dil_master[i][row_min_aus:row_max_aus+1]
    temp_analysis_aus=temp_master[i][row_min_aus:row_max_aus+1]
    aus_line=P.polyfit(temp_analysis_aus,dil_analysis_aus,1)
    CTE_aus_master.append(aus_line[1]/file_list[i][5])
    aus_line_master.append(aus_line)
    strain_at_normalizing_temp_master.append(P.polyval(normalizing_temp,aus_line)/file_list[i][5])
    CTE_gama=aus_line[1]/file_list[i][5]  
    CTE_tot=CTE_tot+CTE_gama

CTE_avr=CTE_tot/(i+1)
delta_L_avr=delta_L_tot/(i+1)
strain_at_normalizing_T_avr=np.average(strain_at_normalizing_temp_master)

CR_master=[]
for row in file_list:
    CR=row[0].split('.')[0][2:]
    
    if CR[0]=='0':
        CR_split=CR.split('0')
        i=0
        while CR_split[i].isdigit()==False:
            i+=1
#            at this point, i shows the number of zeros the are in front of the 
#            name of the file which shows the number of decimal point e.g. CC02.csv 
#            means cooling rate of 0.2 C/s
        CR=float(CR)*10**(-1*i)
    else:
        CR=float(CR)
    CR_master.append(CR)
        
plt.figure(figsize=(8,6))
plt.plot(CR_master,np.array(CTE_aus_master)*1e5,'-o')
plt.title('CTE of austenite vs cooling rate')
plt.xlabel('Cooling rate $(\degree C/s)$')
plt.ylabel('CTE $(\\times 10^{-5})$')    
   

def corr_fact_l0(x):
    
    if L0_correction_method==1:        
        err=(5+abs(current_dldT/(file_list[m][5]+x)-CTE_avr))**5-5**5
    if L0_correction_method==2:
        err=(5+abs(np.amax(dil_master[m])/(file_list[m][5]+x)-delta_L_avr))**5-5**5
    if L0_correction_method==3:
        err=(5+abs(strain_at_normalizing_temp_master[m]-strain_at_normalizing_T_avr+x/file_list[m][5]))**5-5**5
    else:
        print ("L0_correction_method not supported. Pick a supported number.")
    return err

L0_correction=[]
if correct_L0.lower()=='yes':
    for m in range(len(file_list)):
        row_min_aus=int(file_list[m][1])
        row_max_aus=int(file_list[m][2])
    #    time_analysis_aus=time[row_min_aus:row_max_aus+1]
        dil_analysis_aus=dil_master[m][row_min_aus:row_max_aus+1]
        temp_analysis_aus=temp_master[m][row_min_aus:row_max_aus+1]
        current_dldT=P.polyfit(temp_analysis_aus,dil_analysis_aus,1)[1] 
        
        x0=0#,a0_gama,CTE_0_gama,CTE_alpha_b,CTE_alpha_a]
        for i in range(5):
            res = scipy.optimize.minimize(corr_fact_l0,x0, method='Nelder-Mead')
            x0= res.x
            print (x0)
        L0_correction.append(x0[0])
    L0_correction=np.array(L0_correction)
    print (L0_correction)
else:
    L0_correction=np.zeros(len(file_list))
CTE_eq_order=1 #Highest power in the CTE polinomial

for i in range(len(file_list)):
    L0=file_list[i][5]+L0_correction[i]
    row_min_fer=int(file_list[i][3])
    row_max_fer=int(file_list[i][4])

    dil_analysis_fer=dil_master[i][row_min_fer:row_max_fer+1]
    temp_analysis_fer=temp_master[i][row_min_fer:row_max_fer+1]
    CTE_product_coef=P.polyder(P.polyfit(temp_analysis_fer,dil_analysis_fer,CTE_eq_order+1))/L0
    CTE_product_master.append(CTE_product_coef)

CTE_temp=[]
temperature_temp=[]

plt.figure(figsize=(8,6)) 
plt.title("CTE product") 
for i in range(len(file_list)):
    print (i , "   " , i)
    temp=np.arange(np.round(temp_master[i][int(file_list[i][4])]),np.round(temp_master[i][int(file_list[i][3])]))
    CTE=P.polyval(temp,CTE_product_master[i])     
    plt.plot(temp,1e5*CTE,'.',label=file_list[i][0])
    CTE_temp=CTE_temp+list(CTE)
    temperature_temp=temperature_temp+list(temp)
    
CTE_product_avr_coef=P.polyfit(temperature_temp,CTE_temp,CTE_eq_order)
plt.plot(range(100,500),1e5*P.polyval(list(range(100,500)),CTE_product_avr_coef),"--",label='Average CTE')
plt.legend()
plt.xlabel('Temperature $(\degree C)$')
plt.ylabel('CTE$ \\times 10^5$ $({\degree C}^{-1})$')
if use_avr_CTE_product==1:
    CTE_alpha_c=CTE_product_avr_coef[0]
    CTE_alpha_b=CTE_product_avr_coef[1]
k=0 #k is the number that shows location of the current data file in the file_list

Mn= chemistry0['mn']
Ni= chemistry0['ni']
Cr= chemistry0['cr']
Mo= chemistry0['mo']
Si= chemistry0['si']
Cu= chemistry0['cu']

precision=4 # This is a number of decimal points used in rounding carbon content for Ms and Bs calculation.

def Bs(C):# C is in mole fraction
#    print C
    strC=str(C)[:(precision+4)]    
    try:
        Bs=Bs_master_dic[strC] 
#        print "Bs from dic"
    except:
        C=mf_to_wp(C)
        if Bs_model==1: # Imperical equation from Yong cok Lee's paper
            Bs=745-110*C-59*Mn-39*Ni-68*Cr-106*Mo+17*Mn*Ni+6*Cr**2+29*Mo**2
        elif Bs_model==2:
            Bs=-23.3882*C**3 + 59.7411*C**2 - 359.987*C + 653.967 #works only for Nasseh's X80-N steel
    Bs_temp_dic[strC]=Bs
#        print "Ms_model must be a number between 1-4"
    return Bs

def Ms(C):# C is in mole fraction ref http://www.lucefin.com/en/siderurgia/area-tecnica/trattamenti-termici/
    strC=str(C)[:(precision+3)]    
    try:
        Ms=Ms_master_dic[strC]
    except:
        C=mf_to_wp(C)
        if Ms_model==1:
            Ms=539-423*C-30.4*Mn-17.7*Ni-12.1*Cr-7.5*Mo#-11.3*Cu-14.5*Si Ref Andrews equation       
        elif Ms_model==2:
            Ms=512-453*C-16.9*Ni-9.5*Mo+212*C**2-71.5*C*Mn+15*Cr+67.6*C*Cr
        elif Ms_model==3:
            Ms=764.2-302.6*C-30.6*Mn-16.6*Ni-8.9*Cr-2.4*Mo-11.3*Cu-14.5*Si-273.15
        elif Ms_model==4:
            Ms=540-584.9*C-23.1*Si-117.7*Mn-42.5*Cr+49.9*Mo-62.5*sqrt(C*Si)+178.3*sqrt(C*Mn)\
            -10*sqrt(C*Cr)+52.5*sqrt(C*Mo)+117.2*sqrt(Si*Mn)+50.9*sqrt(Si*Cr)-142.2*sqrt(Si*Mo)\
            -29.2*sqrt(Mn*Cr)-9.7*sqrt(Mn*Mo)+69.9*sqrt(Cr*Mo)
        else:
            pass
        Ms_temp_dic[strC]=Ms
    return Ms


def Bainite_ss_factor(T,C): #using relaxation observed in isothermal experiments, C is the amount of C in Austenite.
    if SSF_method.lower()=="general":
        SSF_coef=P.polyfit([Bs(C)-100,Ms(C)],[0,1],1) #this has to be fixed. lower Bs has to be used not just Bs which is for upper bainite. 100C is subtracted from Bs as an estimate for lower Bs start temp.
        BSSF=P.polyval(T,SSF_coef)
        if BSSF<0: BSSF=0 
        elif BSSF>1: BSSF=1
    elif SSF_method.lower()=="user":
        SSF_coef=np.array([ -4.82210399e+01,   2.19859596e-01,  -2.45490039e-04])#using relaxation observed in isothermal experiments, C is the amount of C in Austenite.       
        if T>512.5: BSSF=0
        elif T<450: BSSF=1
        else:
            BSSF=P.polyval(T,SSF_coef)
            if BSSF<0: BSSF=0
            if BSSF>1: BSSF=1
    else: 
        print ("SSF_method must be either general or custom")
    return BSSF

@lru_cache
def C_in_alpha(T):
    C= 1.4734491E-20*T**6 + 3.9638142E-17*T**5 - 1.1293268E-13*T**4 + 6.8406210E-11*T**3 - 9.3489472E-09*T**2 + 6.1810195E-07*T - 6.3920771E-06
    return C


def L(C): # L is an operator that changes C=N_C/(N_C+N_Fe) to N_C/N_Fe. 0<C<1
    return C/(1-C)

output = mp.Queue()

def Fitter(filename, output, a0_gama, CTE_alpha_a, CTE_alpha_b, CTE_alpha_c, a0_alpha,  c_wf_for_cte, Bs_master_dic, MF_dic, Ms_master_dic, C_in_alpha_master_dic):

    N_total=0
    k=0
    for var in file_list:
        if var[0]==filename:
            break
        else:
            k+=1
#    print "k=", k
    L0=file_list[k][5]+L0_correction[k]
#    time=np.array(time_master[k])
#    temperature=np.array(temp_master[k])
#    dilation=np.array(dil_master[k])

    time=time_master[k]
    temperature=temp_master[k]
    dilation=dil_master[k]

    # determine pure austenite range
    row_min_aus=int(file_list[k][1])
    row_max_aus=int(file_list[k][2])
#    time_analysis_aus=time[row_min_aus:row_max_aus+1]
    dil_analysis_aus=dilation[row_min_aus:row_max_aus+1]
    temp_analysis_aus=temperature[row_min_aus:row_max_aus+1]


#$$$$$$$$$$$$$$$$$$$$$$$$$    
    #lets calculate the CTE_0_gama using raw data
    CTE_actual_gama=P.polyfit(temp_analysis_aus,dil_analysis_aus,1)[1]/L0
    #CTE_0_gama is the CTE of austenite at zero mole fraction carbon
    CTE_0_gama=CTE_actual_gama+0.5E-6*C0*100  
#$$$$$$$$$$$$$$$$$$$$$$$$$ 
     
    row_min_fer=int(file_list[k][3])
    row_max_fer=int(file_list[k][4])

    #determin the range of interest for fraction transformed calculation.
    row_min= row_min_aus
    row_max=row_max_fer

    time_analysis=time[row_min:row_max+1]
    dil_analysis=dilation[row_min:row_max+1]
    temp_analysis=temperature[row_min:row_max+1]

    def Vgama (C,T):
        return ((a0_gama + 6.5e-4 * 1e-9 * C * 100) * (1 + (CTE_0_gama - 0.5E-6 * C * 100) * (T-726.85)))**3
        
#    def VCement(T): #T in centigrade. This function can be optimized for faster run
#        y=(1+(5.311E-6-1.942E-9*T+9.655E-12*T*T)*(T-20))
#        a=0.45234E-9*y
#        b=0.50883E-9*y
#        c=0.67426E-9*y
#        VCementite=a*b*c
#        return VCementite
        
#    @jit(nopython=True)   
    def C_in_gama(timestep,current_dn_alpha,dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn): #amount of C in gama on the timestep point
        C_total_in_ferrite=sum(dn_alpha[:timestep-1]*L(C_dn[:timestep-1]))+current_dn_alpha*L(C_dn[timestep-1]) #this calculates total number of C atoms in ferrite phase
        if ID_dn[timestep-1][3:6]=="Cem":
            current_dn_cem=((N_total-N_product[timestep-1])*(L(C_in_aus_molefraction[timestep-1])-L(Solubility_aus_cement(temp_fit[timestep])))-current_dn_alpha*(L(C_dn[timestep-1])-L(Solubility_aus_cement(temp_fit[timestep]))))/(L(0.25)-L(Solubility_aus_cement(temp_fit[timestep])))
        else:
            current_dn_cem=0
        C_total_in_cement=sum(dn_cement[:timestep-1]/3.0)+current_dn_cem/3.0
        N_gama=N_total-(N_product[timestep-1]+current_dn_alpha+current_dn_cem)
        C_in_gama_local=L(C0)*N_total-C_total_in_ferrite-C_total_in_cement #No. of C atoms in aus.
        return C_in_gama_local/(N_gama+C_in_gama_local)
    
    def fundamental_parameters(x):
        N_t=x
        modeled_l_gama=((N_t/4)*Vgama(C0,temp_analysis_aus))**(1/3.0)
        residual_gama=np.sum(np.abs(modeled_l_gama-(dil_analysis_aus+L0)))
        err=residual_gama
        return err

    x0=[3e22]#initial estimate of x0 (here it represents N_t) is very important. Too big of number will give error.
    for i in range(3):
        res = scipy.optimize.minimize(fundamental_parameters,x0, method='Nelder-Mead')
        x0= res.x
    N_total=x0[0]
    
    
    
    print ('N_total of %s=%s'%(filename,N_total))
#    print 'CTE_0_gama of %s=%s'%(filename,CTE_0_gama)
#    print
 
#    raw_input('Hit enter to begin fitting points')

    time_fit=np.zeros((row_max-row_min-2*sr)//interval)
    temp_fit=np.zeros((row_max-row_min-2*sr)//interval)
    dil_fit=np.zeros((row_max-row_min-2*sr)//interval)
#    dldt_fit=np.zeros((row_max-row_min-2*sr)//interval)
    num_of_analized_points=len(time_fit)

    n=sr
    
    # TODO all this should be done once out of the fitter. reduces amount of calculations.
    for i in range(num_of_analized_points):
        #print'i= ',i
        time_fit[i]=time_analysis[n]
        coef1= P.polyfit(time_analysis[n-sr:n+sr],temp_analysis[n-sr:n+sr],reg_order)
        temp_fit[i]=P.polyval(time_analysis[n],coef1)
        coef1= P.polyfit(time_analysis[n-sr:n+sr],dil_analysis[n-sr:n+sr],reg_order)
        dil_fit[i]=P.polyval(time_analysis[n],coef1)
        n=n+interval

        
    def expm_l_individual(x,i, return_type = 'modeled_L'):
        previ_alpha_v=0.0
        current_step_alpha=0.0
        previ_cem_v=0.0
        current_step_cem=0.0
        dn_cem=0.0
        for j in range(i-1): #This for loop calculates volume of all the previous steps. i-1 is correct
            if ID_dn[j][0:2]=="PF":
                previ_alpha_v=previ_alpha_v+dn_alpha[j]/2*Valpha(a0_alpha,C_dn[j],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)
            if ID_dn[j][0:2]==("BF"):
                previ_alpha_v=previ_alpha_v+dn_alpha[j]/2*Vbainite(a0_alpha,C_dn[j],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)
            
            if ID_dn[j][0:2]==("MS"):
                previ_alpha_v=previ_alpha_v+dn_alpha[j]/2*Vbainite(a0_alpha,C_dn[j],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)
      
            
            if ID_dn[j][3:6]=="Cem":
                previ_cem_v=previ_cem_v+dn_cement[j]/12*VCement(temp_fit[j])
        if ID_dn[i-1][0:2]=="PF":
            current_step_alpha=x/2*Valpha(a0_alpha,C_dn[i-1],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)
        if ID_dn[i-1][0:2]==("BF"):
            current_step_alpha=x/2*Vbainite(a0_alpha,C_dn[i-1],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)
        
        if ID_dn[i-1][0:2]==("MS"):
            current_step_alpha=x/2*Vbainite(a0_alpha,C_dn[i-1],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)
   
        
        if ID_dn[i-1][3:6]=="Cem":
            AA=L(Solubility_aus_cement(temp_fit[i]))
            dn_cem=((N_total-N_product[i-1])*(AA-L(C_in_aus_molefraction[i-1]))+x*(L(C_dn[i-1])-AA))/(-1/3.0+AA)            
            current_step_cem=dn_cem/12*VCement(temp_fit[i]) #12 iron atoms in one unit cell of cementite.
        current_gama_v=(N_total-sum(dn_alpha[:i-1])-x-sum(dn_cement[:i-1])-dn_cem)/4*Vgama(C_in_gama(i,x,dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn),temp_fit[i]) #Vgama should be fixed to get rid of x/N_total
        modeled_V=previ_alpha_v+previ_cem_v+current_step_alpha+current_step_cem+current_gama_v
        #print (raw_input('press any key to continue in individual loop'))
        if return_type == 'modeled_L':
            return (modeled_V)**(1/3.0)
        elif return_type == 'error':
            err=abs((dil_fit[i]+L0)**3-modeled_V)*10e19
            return err
        
    dn_alpha=np.ones(len(temp_fit))*10**(np.log10(N_total)-4) #this is the starting point of calculation which is close enough to zero to be considered zero fraction transformed, yet large enough to allow for optimization.
    dn_cement=np.zeros(len(temp_fit))
    #C_dn_alpha=abs((C_in_alpha((temp_fit[:-1]+temp_fistart=time.time()t[1:])/2))) #This is the C mole fraction associated with each dn
    C_dn=np.zeros(len(temp_fit))
    ID_dn=[]
    #print (raw_input('press any key to continue'))
    simulated_l=np.zeros(len(dn_alpha))
    C_in_aus=np.ones(len(dn_alpha))*C0
    C_in_aus_molefraction=np.ones(len(dn_alpha))*C0
#    CTE_aus_current=np.zeros(len(dn_alpha))
    N_product=np.zeros(len(temp_fit))
#    N_prod=np.zeros(len(temp_fit))

    for j in range(1):
        #print '******************'
#        i=1 # start filling dn_alpha from the third element since the first two should be zero. the first
        for i in range(1,num_of_analized_points):
            T=(temp_fit[i-1]+temp_fit[i])/2

#            print i, "out of", len(temp_fit)-1 ,"points completed"
            # at every i we calculated the dn[i-1]

            # First we shoud determine what is the type of the product forming between i-1 and i. 
            # We do that by assuming C_in_austenite during current step is equal to that of the last step.
            if i>1:
                if temp_fit[i-1]<Ms(C_in_aus[i-1]):
                    ID_dn.append("MS")                    
                elif temp_fit[i-1]<Bs(C_in_aus[i-1]):
#                    print "   Bainite forming according to empirical formula"
                    ID_dn.append("BF")
                else:
                    ID_dn.append("PF")

                if ID_dn[i-1][0:2]=="PF":
                    C_dn[i-1]=C_in_alpha(T)
                    if C_dn[i-1] > C_in_gama (i-1, dn_alpha[i-1], dn_alpha, C_dn, N_total, N_product, temp_fit, ID_dn):
                        C_dn[i-1]=C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)
                elif ID_dn[i-1][0:2]=="BF":
                    C_eq=C_in_alpha(T)
                    C_dn[i-1]=Bainite_ss_factor(T,C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn))*(C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)-C_eq)+C_eq
                    if C_dn[i-1] > C_in_gama(i-1, dn_alpha[i-1], dn_alpha, C_dn, N_total, N_product, temp_fit, ID_dn):
                        C_dn[i-1]=C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)
                elif ID_dn[i-1][0:2]=="MS":
                    C_dn[i-1]=C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)
                    
                if ('yes' or 'y') in file_list[k][6].lower() and ID_dn[i-1]!="MS" and C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)>Solubility_aus_cement(temp_fit[i-1]):
                    ID_dn[i-1]=ID_dn[i-1]+"+"+"Cem"
                elif (("user" in file_list[k][6].lower()) and ID_dn[i-1]!="MS" and (C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)>Solubility_aus_cement(temp_fit[i-1]))):
                    LCT=[float(s) for s in file_list[k][7].split() if s.isdigit()][0]#LCT stands for Lowest cementite tempereature
                    if T>LCT:
                        ID_dn[i-1]=ID_dn[i-1]+"+"+"Cem"

            elif i==1:
                ID_dn.append("PF")
                C_dn[i-1]=C_in_alpha((temp_fit[i-1]+temp_fit[i])/2)
            else:
                print ("There is a problem with ID_dn")
                quit()

            res=scipy.optimize.minimize(expm_l_individual, dn_alpha[i-1], args = (i, 'error'), method='nelder-mead', options={})
            K=res.x
            dn_alpha[i-1]=K[0]     
            
            #lets put a constraint on dn to control noise which result in incorrect negative dn            
            if ((dn_alpha[i-1])/N_total)<1e-3:
                dn_alpha[i-1]=0
                
            if ID_dn[i-1][3:6]=="Cem":
                AA=L(Solubility_aus_cement(temp_fit[i]))
                dn_cement[i-1]=((N_total-N_product[i-1])*(-L(C_in_aus_molefraction[i-1])+AA)+dn_alpha[i-1]*(L(C_dn[i-1])-AA))/(-1/3.0+AA)
                
#            if "MS" in ID_dn:
#                print dn_alpha[i-1]
#        
            N_product[i]=N_product[i-1]+dn_alpha[i-1]+dn_cement[i-1]
            #print (raw_input('press any key to go to next point in main loop'))
            simulated_l[i]= expm_l_individual(dn_alpha[i-1],i)

            C_in_aus[i]=C_in_gama(i,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)
            C_in_aus_molefraction[i]=C_in_gama(i,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)


    no_end_points=int((row_max_fer-row_min_fer)//interval)
    coef= P.polyfit(temp_fit[len(temp_fit)-no_end_points:],N_product[len(temp_fit)-no_end_points:],1)
    err_end_slope=(5+abs(coef[1]/N_total))**err_end_slop_WF-5**err_end_slop_WF
    # err_end_fit issues error when fit is not good at the ferrite end. sometimes ferrite end has fewer point that makes them less effective in overal error. 
    # separeting them allows for increasing there value in overal error 
    err_end_fit=(5+sum(abs((dil_fit[len(temp_fit)-no_end_points:]+L0)-simulated_l[len(temp_fit)-no_end_points:]))/L0)**end_fit_WF-5**end_fit_WF
    err_end_fraction=0# this issues error if end fraction goes above 1
    Numer_of_high_point=int(0)
    for i in range(len(temp_fit)//3,len(temp_fit)): 
#        print i
#        print "len(temp_fit)",len(temp_fit)
        if N_product[i]>N_total:
            Numer_of_high_point+=1
            err_end_fraction+=(5+((N_product[i]-N_total)/N_total))**20-5**20
    
    if (maximize_fraction_transformed.lower()==('yes'or 'y')):
        err_maximum_transformation=(5+sum(abs(N_product[:-3]/N_total-1)))**err_maximum_transformation_WF-5**err_maximum_transformation_WF
    else: err_maximum_transformation=0

    err_overal_fit=(5+sum(abs((dil_fit[1:-1]+L0)-simulated_l[1:-1]))/L0)**overal_fit_WF-5**overal_fit_WF
#    print 'err_pure_phase=' ,err_pure_phase
#    print "err_alpha_fit",err_alpha_fit
    total_cost= err_overal_fit+ err_end_slope + err_end_fraction + err_end_fit + err_maximum_transformation
#    err_overal_fit=err_alpha_fit


    output.put([total_cost,Bs_temp_dic,Ms_temp_dic,MF_temp_dic,C_in_alpha_temp_dic])#+err_pure_phase)#+err_gama_fit+err_alpha_fit)

def Fitter_plot(filename,output,a0_gama,CTE_alpha_a,CTE_alpha_b,CTE_alpha_c,a0_alpha,c_wf_for_cte, c_wf_for_a0):
    print ('----------V10_plot-------------')

    N_total=0
    k=0
    for var in file_list:
        if var[0]==filename:
            break
        else:
            k+=1
#    print "k=", k
    L0=file_list[k][5]+L0_correction[k]
#    time=np.array(time_master[k])
#    temperature=np.array(temp_master[k])
#    dilation=np.array(dil_master[k])

    time=time_master[k]
    temperature=temp_master[k]
    dilation=dil_master[k]
    # determine pure austenite range
    row_min_aus=int(file_list[k][1])
    row_max_aus=int(file_list[k][2])
#    time_analysis_aus=time[row_min_aus:row_max_aus+1]
    dil_analysis_aus=dilation[row_min_aus:row_max_aus+1]
    temp_analysis_aus=temperature[row_min_aus:row_max_aus+1]


#$$$$$$$$$$$$$$$$$$$$$$$$$    
    #lets calculate the CTE_0_gama using raw data
    CTE_gama=P.polyfit(temp_analysis_aus,dil_analysis_aus,1)[1]/L0
    #CTE_0_gama is the CTE of austenite at zero mole fraction carbon
    CTE_0_gama=CTE_gama+0.5e-6*C0*100   
#$$$$$$$$$$$$$$$$$$$$$$$$$ 


#    row_min_fer=file_list[k][3]
    row_min_fer=int(file_list[k][3])
    row_max_fer=int(file_list[k][4])

    dil_analysis_fer=dilation[row_min_fer:row_max_fer+1]
    temp_analysis_fer=temperature[row_min_fer:row_max_fer+1]
    CTE_product_coef=P.polyder(P.polyfit(temp_analysis_fer,dil_analysis_fer,3))/L0
    print (CTE_product_coef)
#    print "CTE_product for",file_list[k][0],'=', CTE_product  
    #determin the range of interest for fraction transformed calculation.
    row_min= row_min_aus
    row_max=row_max_fer

    # create experimentally measured dv/dt for every data point in analysis range

    time_analysis=time[row_min:row_max+1]
    dil_analysis=dilation[row_min:row_max+1]
    temp_analysis=temperature[row_min:row_max+1]

##
#    def C_in_alpha(T): #Mole fraction of C in alpha using thermocalc under para equilibrium
#        C= 1.4734491E-20*T**6 + 3.9638142E-17*T**5 - 1.1293268E-13*T**4 + 6.8406210E-11*T**3 - 9.3489472E-09*T**2 + 6.1810195E-07*T - 6.3920771E-06
#        return abs(C)

    def Vgama (C,T):
        return ((a0_gama + 6.5e-4 * 1e-9 * C * 100) * (1 + (CTE_0_gama - 0.5E-6 * C * 100) * (T-726.85)))**3
        
#    def Vgama_param(a0_gama,CTE_0_gama,T):
#        Vgama=((a0_gama+6.5e-4*1e-9*C0*100)*(1+(CTE_0_gama-0.5E-6*C0*100)*(T-726.85)))**3
#        return Vgama
#    
#    def Vgama_new (CTE_0_gama,C,T): # check if C is mole fraction
#        Vgama=((a0_gama+6.5e-4*1e-9*C*100)*(1+(CTE_0_gama-0.5E-6*C*100)*(T-726.85)))**3
#        return Vgama

    

#    
    def C_in_gama(timestep,current_dn_alpha,dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn): #amount of C in gama on the timestep point
        C_total_in_ferrite=sum(dn_alpha[:timestep-1]*L(C_dn[:timestep-1]))+current_dn_alpha*L(C_dn[timestep-1]) #this calculates total number of C atoms in ferrite phase
        if ID_dn[timestep-1][3:6]=="Cem":
            current_dn_cem=((N_total-N_product[timestep-1])*(L(C_in_aus_molefraction[timestep-1])-L(Solubility_aus_cement(temp_fit[timestep])))-current_dn_alpha*(L(C_dn[timestep-1])-L(Solubility_aus_cement(temp_fit[timestep]))))/(L(0.25)-L(Solubility_aus_cement(temp_fit[timestep])))
        else:
            current_dn_cem=0
        C_total_in_cement=sum(dn_cement[:timestep-1]/3.0)+current_dn_cem/3.0
        N_gama=N_total-(N_product[timestep-1]+current_dn_alpha+current_dn_cem)
        C_in_gama_local=L(C0)*N_total-C_total_in_ferrite-C_total_in_cement #No. of C atoms in aus.
        return C_in_gama_local/(N_gama+C_in_gama_local)
    
    def fundamental_parameters(x):
        N_t=x
        modeled_l_gama=((N_t/4)*Vgama(C0,temp_analysis_aus))**(1/3.0)
        residual_gama=np.sum(np.abs(modeled_l_gama-(dil_analysis_aus+L0)))
        err=residual_gama
        return err

    x0=[5e22]#,a0_gama,CTE_0_gama,CTE_alpha_b,CTE_alpha_a]
    for i in range(3):
        res = scipy.optimize.minimize(fundamental_parameters,x0, method='Nelder-Mead')
        x0= res.x
    N_total=x0[0]
    
    print ('N_total of %s=%s'%(filename,N_total))
    print ('CTE_0_gama of %s=%s'%(filename,CTE_0_gama))
#    print
#    input('Hit enter to begin fitting points')
############################################
# different from fitter


    simulated_aus=np.zeros(len(temp_analysis))
    for i in range (0, len(temp_analysis)):
        simulated_aus[i]= (N_total/4*Vgama(C0,temp_analysis[i]))**(1/3.0)

    plt.figure(7,figsize=(8,6))
    plt.plot(temp_analysis,1e3*(dil_analysis+L0),'.', label="experimental data -" + filename)
    plt.legend(loc=2)
    plt.title('Fit of austenite')
    plt.plot(temp_analysis,1e3*simulated_aus, '-',label="Fitted austenite - " + filename)
    plt.xlabel('Temperature$(\\degree C)$')
    plt.ylabel('Specimen length $(mm)$')
    plt.legend(loc=2)
    plt.grid()    
############################################
    
    
    time_fit=np.zeros((row_max-row_min-2*sr)//interval)
    temp_fit=np.zeros((row_max-row_min-2*sr)//interval)
    dil_fit=np.zeros((row_max-row_min-2*sr)//interval)
#    dldt_fit=np.zeros((row_max-row_min-2*sr)//interval)
    num_of_analized_points=len(time_fit)

    n=sr
    for i in range(num_of_analized_points):
        #print'i= ',i
        time_fit[i]=time_analysis[n]
        coef1= P.polyfit(time_analysis[n-sr:n+sr],temp_analysis[n-sr:n+sr],reg_order)
        temp_fit[i]=P.polyval(time_analysis[n],coef1)
        coef1= P.polyfit(time_analysis[n-sr:n+sr],dil_analysis[n-sr:n+sr],reg_order)
        dil_fit[i]=P.polyval(time_analysis[n],coef1)
        n=n+interval


    def expm_l_individual_err(x):
        previ_alpha_v=0.0
        current_step_alpha=0.0
        previ_cem_v=0.0
        current_step_cem=0.0
        dn_cem=0.0
        for j in range(i-1): #This for loop calculates volume of all the previous steps. i-1 is correct
            if ID_dn[j][0:2]=="PF":
                previ_alpha_v=previ_alpha_v+dn_alpha[j]/2*Valpha(a0_alpha,C_dn[j],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte,c_wf_for_a0, MF_to_WP)

            if ID_dn[j][0:2]==("BF"):
                previ_alpha_v=previ_alpha_v+dn_alpha[j]/2*Vbainite(a0_alpha,C_dn[j],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)

            if ID_dn[j][0:2]==("MS"):
                previ_alpha_v=previ_alpha_v+dn_alpha[j]/2*Vbainite(a0_alpha,C_dn[j],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)
  

            if ID_dn[j][3:6]=="Cem":
                previ_cem_v=previ_cem_v+dn_cement[j]/12*VCement(temp_fit[j])
        if ID_dn[i-1][0:2]=="PF":
            current_step_alpha=x/2*Valpha(a0_alpha,C_dn[i-1],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)

        if ID_dn[i-1][0:2]==("BF"):
            current_step_alpha=x/2*Vbainite(a0_alpha,C_dn[i-1],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)
       
        if ID_dn[i-1][0:2]==("MS"):
            current_step_alpha=x/2*Vbainite(a0_alpha,C_dn[i-1],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)

        if ID_dn[i-1][3:6]=="Cem":
            AA=L(Solubility_aus_cement(temp_fit[i]))
            dn_cem=((N_total-N_product[i-1])*(AA-L(C_in_aus_molefraction[i-1]))+x*(L(C_dn[i-1])-AA))/(-1/3.0+AA)
#            dn_cem=(1/(1/3.0-L(Solubility_aus_cement(temp_fit[i]))))*((dn_alpha[:i-1]*C_dn[:i-1])+(np.sum(dn_cem[:i-1])/3.0)+(N_total-N_product[i])*L(Solubility_aus_cement(temp_fit[i]))+x*(L(C_in_alpha((temp_fit[i]+temp_fit[i-1])/2.0)))-N_total*(L(C0)))
            current_step_cem=dn_cem/12*VCement(temp_fit[i]) #12 iron atoms in one unit cell of cementite.
        current_gama_v=(N_total-sum(dn_alpha[:i-1])-x-sum(dn_cement[:i-1])-dn_cem)/4 *Vgama(C_in_gama(i,x,dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn),temp_fit[i]) 
        modeled_V=previ_alpha_v+previ_cem_v+current_step_alpha+current_step_cem+current_gama_v
        err=abs((dil_fit[i]+L0)**3-modeled_V)*10e19
        return err
        
    def expm_l_individual(x,i):
        previ_alpha_v=0.0
        current_step_alpha=0.0
        previ_cem_v=0.0
        current_step_cem=0.0
        dn_cem=0.0
        for j in range(i-1): #This for loop calculates volume of all the previous steps. i-1 is correct
            if ID_dn[j][0:2]=="PF":
                previ_alpha_v=previ_alpha_v+dn_alpha[j]/2*Valpha(a0_alpha,C_dn[j],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)
            if ID_dn[j][0:2]==("BF"):
                previ_alpha_v=previ_alpha_v+dn_alpha[j]/2*Vbainite(a0_alpha,C_dn[j],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)
            
            if ID_dn[j][0:2]==("MS"):
                previ_alpha_v=previ_alpha_v+dn_alpha[j]/2*Vbainite(a0_alpha,C_dn[j],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)
      
            
            if ID_dn[j][3:6]=="Cem":
                previ_cem_v=previ_cem_v+dn_cement[j]/12*VCement(temp_fit[j])
        if ID_dn[i-1][0:2]=="PF":
            current_step_alpha=x/2*Valpha(a0_alpha,C_dn[i-1],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)
        if ID_dn[i-1][0:2]==("BF"):
            current_step_alpha=x/2*Vbainite(a0_alpha,C_dn[i-1],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)
        
        if ID_dn[i-1][0:2]==("MS"):
            current_step_alpha=x/2*Vbainite(a0_alpha,C_dn[i-1],temp_fit[i],CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)
   
        
        if ID_dn[i-1][3:6]=="Cem":
            AA=L(Solubility_aus_cement(temp_fit[i]))
            dn_cem=((N_total-N_product[i-1])*(AA-L(C_in_aus_molefraction[i-1]))+x*(L(C_dn[i-1])-AA))/(-1/3.0+AA)            
            current_step_cem=dn_cem/12*VCement(temp_fit[i]) #12 iron atoms in one unit cell of cementite.
        current_gama_v=(N_total-sum(dn_alpha[:i-1])-x-sum(dn_cement[:i-1])-dn_cem)/4*Vgama(C_in_gama(i,x,dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn),temp_fit[i]) #Vgama should be fixed to get rid of x/N_total
        modeled_V=previ_alpha_v+previ_cem_v+current_step_alpha+current_step_cem+current_gama_v
        #print (raw_input('press any key to continue in individual loop'))
        return (modeled_V)**(1/3.0)

      
    dn_alpha=np.ones(num_of_analized_points-1)*10**(np.log10(N_total)-4)
    dn_cement=np.zeros(num_of_analized_points-1)
    
    ###########################different from fitter#############
    lattice_param_gamma_during_trans=np.zeros(num_of_analized_points)
    #############################################################
    
    #C_dn_alpha=abs((C_in_alpha((temp_fit[:-1]+temp_fistart=time.time()t[1:])/2))) #This is the C mole fraction associated with each dn
    C_dn=np.zeros(num_of_analized_points-1)
    ID_dn=[]
    #print (raw_input('press any key to continue'))
    simulated_l=np.zeros(num_of_analized_points)
    C_in_aus=np.ones(num_of_analized_points)*C0
    C_in_aus_molefraction=np.ones(num_of_analized_points)*C0
#    CTE_aus_current=np.zeros(len(dn_alpha))
    N_product=np.zeros(num_of_analized_points)
#    N_prod=np.zeros(num_of_analized_points)


    ###########################different from fitter#############
    Vol_dn_alpha_at_20=np.zeros(num_of_analized_points-1)
    Vol_dn_cem_at_20=np.zeros(num_of_analized_points-1)
    Vol_transformed=np.zeros(num_of_analized_points-1)
    ###########################different from fitter#############
    
    for j in range(1):
        #print '******************'
#        i=1 # start filling dn_alpha from the third element since the first two should be zero. the first
        for i in range(1,num_of_analized_points):
            ii=i-1
            T=(temp_fit[ii]+temp_fit[i])/2
#            print i, "out of", num_of_analized_points-1 ,"points completed"
            # at every i we calculated the dn[i-1]
#            print (raw_input('press any key to continue in individual loop'))            # First we shoud determine what is the type of the product forming between i-1 and i. We do that by assuming C_in_austenite during current step is equal to that of the last step.
            if i>1:
                if temp_fit[i-1]<Ms(C_in_aus[i-1]):
                    ID_dn.append("MS")                    
                elif temp_fit[i-1]<Bs(C_in_aus[i-1]):
#                    print "   Bainite forming according to empirical formula"
                    ID_dn.append("BF")
                else:
                    ID_dn.append("PF")

                if ID_dn[i-1][0:2]=="PF":
                    C_dn[i-1]=C_in_alpha(T)
                    if C_dn[i-1]>C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn):
                        C_dn[i-1]=C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)
                elif ID_dn[i-1][0:2]=="BF":
                    C_eq=C_in_alpha(T)
                    C_dn[i-1]=Bainite_ss_factor(T,C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn))*(C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)-C_eq)+C_eq
                    if C_dn[i-1]>C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn):
                        C_dn[i-1]=C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)
                elif ID_dn[i-1][0:2]=="MS":
                    C_dn[i-1]=C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)
                    
                if ('yes' or 'y') in file_list[k][6].lower() and ID_dn[i-1]!="MS" and C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)>Solubility_aus_cement(temp_fit[i-1]):
                    ID_dn[i-1]=ID_dn[i-1]+"+"+"Cem"
                elif (("user" in file_list[k][6].lower()) and ID_dn[i-1]!="MS" and (C_in_gama(i-1,dn_alpha[i-1],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)>Solubility_aus_cement(temp_fit[i-1]))):
                    LCT=[float(s) for s in file_list[k][7].split() if s.isdigit()][0]#LCT stands for Lowest cementite tempereature
                    if T>LCT:
                        ID_dn[i-1]=ID_dn[i-1]+"+"+"Cem"

            elif i==1:
                ID_dn.append("PF")
                C_dn[i-1]=C_in_alpha((temp_fit[i-1]+temp_fit[i])/2)
            else:
                print ("There is a problem with ID_dn")
                quit()
            res=scipy.optimize.minimize(expm_l_individual_err, dn_alpha[ii], method='Nelder-Mead', options={})
            K=res.x
            dn_alpha[ii]=K[0]          
            #lets put a constraint on dn to control noise which result in incorrect negative dn            
            if ((dn_alpha[ii])/N_total)<1e-3:
                dn_alpha[ii]=0
                
            if ID_dn[ii][3:6]=="Cem":
                AA=L(Solubility_aus_cement(temp_fit[i]))
                dn_cement[ii]=((N_total-N_product[ii])*(-L(C_in_aus_molefraction[ii])+AA)+dn_alpha[ii]*(L(C_dn[ii])-AA))/(-1/3.0+AA)
                
            #Calculate volumes at 20C. These volumes will be used in volume fraction calculations.
            Vol_dn_alpha_at_20[ii]=dn_alpha[ii]/2*Valpha(a0_alpha,C_dn[ii],20,CTE_alpha_a,CTE_alpha_b,CTE_alpha_c, c_wf_for_cte, c_wf_for_a0, MF_to_WP)                
            
            Vol_dn_cem_at_20[ii]=dn_cement[ii]/12*VCement(20)

        
            N_product[i]=N_product[ii]+dn_alpha[ii]+dn_cement[ii]
            #print (raw_input('press any key to go to next point in main loop'))
            simulated_l[i]= expm_l_individual(dn_alpha[ii],i)
#            print i
            C_in_aus[i]=C_in_gama(i,dn_alpha[ii],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)
            lattice_param_gamma_during_trans[i]=(Vgama(C_in_aus[i],temp_fit[i]))**(1/3.0)
            C_in_aus_molefraction[i]=C_in_gama(i,dn_alpha[ii],dn_alpha,C_dn,N_total,N_product,temp_fit,ID_dn)
#            C_in_aus[i]= mf_to_wp(C)
            C_total_in_ferrite=sum(dn_alpha[:i]*L(C_dn[:i]))
            C_total_in_cementite=sum(dn_cement[:i]*(1/3.0))
            C_total_in_gama=(N_total-N_product[i])*L(C_in_aus_molefraction[i])
            Initial_C=N_total*L(C0)
            C_balance=(C_total_in_ferrite+C_total_in_cementite+C_total_in_gama)/Initial_C            
            print ("    T is" , temp_fit[i])
            print ("    Current product is ", ID_dn[ii])
            print ("    Fraction transformed is ", "%.3f" %(N_product[i]/N_total))
            print ("    C balance (should be 1)= ", C_balance)
#            print "    dil_fit = ", dil_fit[i]
            print ("----------------------------------")
    Vol_retained_aus_at_20=(N_total-N_product[-1])/4* Vgama(C_in_aus[i],20)
    Total_vol=np.sum(Vol_dn_alpha_at_20) + np.sum(Vol_dn_cem_at_20) + Vol_retained_aus_at_20
    for i in range(len(Vol_dn_alpha_at_20)):
        Vol_transformed[i]=np.sum(Vol_dn_alpha_at_20[0:i])+np.sum(Vol_dn_cem_at_20[0:i])
    Vol_f_transformed= Vol_transformed/Total_vol
    print ("volume fraction of retained austenite =", Vol_retained_aus_at_20/Total_vol)
    vol_f_retained_aus=Vol_retained_aus_at_20/Total_vol
    dataname=filename[:-4]
    f=np.zeros(len(time_fit))
    f=N_product/N_total
    
    
    
    
    plt.figure(8,figsize=(8,7))
    plt.subplot(211)
    plt.plot(temp_analysis,dil_analysis+L0,'.',label= 'Experimental '+ dataname)
    plt.xlabel('Temperature $(\\degree C)$')
    plt.ylabel('Sample length $(m)$')
    plt.grid()


    plt.subplot(211)
    plt.plot(temp_fit[1:-1],simulated_l[1:-1], 'o',label= 'Calculated ' + dataname)
    plt.legend(loc=2)
    plt.subplot(212)
    plt.plot(temp_fit,f,'o-', label= dataname)
    plt.xlabel('Temperature $(\\degree C)$')
    plt.ylabel('Ferrite mole fraction')
    plt.legend()
    plt.grid()

#    plt.figure(44,figsize=(8,7))
#    plt.plot(temp_analysis,dil_analysis+L0,'.',label= 'Experimental '+ dataname)
#    plt.xlabel('Temperature $(\\degree C)$')
#    plt.ylabel('Sample length $(m)$')
#    plt.grid()
##    C=np.polynomial.chebyshev.chebfit(temp_analysis[50:-50],dil_analysis[50:-50]+L0,5)
##    fitted=np.polynomial.chebyshev.chebval(temp_analysis[50:-50],C)
#    
##    C=np.polynomial.legendre.legfit(temp_analysis[10:-50],dil_analysis[10:-50]+L0,50)
##    fitted=np.polynomial.legendre.legval(temp_analysis[10:-50],C)
#    
##    C=np.polynomial.laguerre.lagfit(temp_analysis[10:-50],dil_analysis[10:-50]+L0,30)
##    fitted=np.polynomial.laguerre.lagval(temp_analysis[10:-50],C)
#    
#    C=np.polynomial.hermite.hermfit(temp_analysis[10:-50],dil_analysis[10:-50]+L0,30)
#    fitted=np.polynomial.hermite.hermval(temp_analysis[10:-50],C)
#    
#    plt.plot(temp_analysis[10:-50],fitted,'-',label= 'Chebyshev Approximation '+ dataname)
#    print ("end")
    
    result_master_time.append(time_fit)
    result_master_temp.append(temp_fit)
    result_master_dil.append(simulated_l)
    result_master_fraction.append(f)
    result_master_FD_fraction.append(dn_alpha/N_total)
    result_master_C_in_FD_fraction.append(C_dn)
    result_master_C_in_FD_WP.append(mf_to_wp(C_dn))
    result_master_C_in_gama.append(C_in_aus)
    result_master_cem_fraction.append(dn_cement/N_total)
    result_master_lattice_param_gamma_during_trans.append(lattice_param_gamma_during_trans)
    result_master_ID.append(ID_dn)
    result_master_Vol_f_transformed.append(Vol_f_transformed)
    result_master_FD_vol.append(Vol_dn_alpha_at_20/Total_vol)#calculates volume fraction of FD at 20C
    result_master_cem_vol.append(Vol_dn_cem_at_20/Total_vol)#calculates volume of the phase at 20C
    result_master_aus_vol.append(vol_f_retained_aus)#calculates volume of the phase at 20C    


def master_fiter(x):
    print ('-'*49)
    i=0
    for val in x: # I am a genius! 
        globals()[x0_param_names[i]]=val
        print (x0_param_names[i]+'=' ,"%.30e" %val)
        print()
        i=i+1


    processes = [mp.Process(target=Fitter, args=((file_list[k][0],output,a0_gama,CTE_alpha_a,CTE_alpha_b,CTE_alpha_c,a0_alpha,c_wf_for_cte,Bs_master_dic, MF_dic, Ms_master_dic,C_in_alpha_master_dic ))) for k in range(len(file_list))]

    for p in processes:
        p.start()
    for p in processes:
        p.join()
    results = [output.get() for p in processes]
#    print results
    err=0
    for result in results:
        err+=result[0]
        Bs_master_dic.update(result[1])
        Ms_master_dic.update(result[2])
        MF_dic.update(result[3])
        C_in_alpha_master_dic.update(result[4])
    print ('err=         ' ,"%.30e" %err)
    return err
x0_param_names=[]
x0_param_vals=[]
x0_param_bounds=[]
for key in optimized_param.keys():
    if optimized_param[key]==1:
        x0_param_names.append(key)
        x0_param_vals.append(eval(key))
        x0_param_bounds.append(optimized_param_bounds[key])       
x0_param_vals=np.array(x0_param_vals)

#x0=[a0_gama,CTE_alpha_c,CTE_alpha_b,CTE_alpha_a,a0_alpha,c_wf_for_cte]
##master_fiter(x0)
#res = scipy.optimize.minimize(master_fiter,x0, method='Nelder-Mead')
#B= res.x

if __name__=='__main__': 
#    global Bs_master_dic, MF_dic, Ms_master_dic
#    Bs_master_dic=mp.Manager().dict()
#    MF_dic=mp.Manager().dict()
#    Ms_master_dic=mp.Manager().dict()
    

    if optimize.lower()=='yes':
        if optimization_method==1:
            print ("Optimizing using the Nelder-Mead method")
            res = scipy.optimize.minimize(master_fiter,x0_param_vals, method='Nelder-Mead',options={'ftol': 0.0001})
        elif optimization_method==2: 
            print ("Optimizing using the differential evolution method")
            res = scipy.optimize.differential_evolution(master_fiter,x0_param_bounds, strategy='best1exp')
        else: 
            print ("Error. Choose 1 or 2 as optimization method!")
    
    #    print '&&&&&&&&&&&&&&& End of round one &&&&&&&&&&&&&&&'
        x0_param_vals=res.x
    
    i=0
    for val in x0_param_vals:
        locals()[x0_param_names[i]]=val
        i=i+1
    
    
    for k in range(len(file_list)):
        Fitter_plot(file_list[k][0],output,a0_gama,CTE_alpha_a,CTE_alpha_b,CTE_alpha_c,a0_alpha,c_wf_for_cte,c_wf_for_a0)
    print ("*-"*25)
    print ('a0_gama =    ' ,"%.30e" %a0_gama)
    print ('a0_alpha=    ' ,"%.30e" %a0_alpha)
    print ('CTE_alpha_a= ' ,"%.30e" %CTE_alpha_a)
    print ('CTE_alpha_b= ' ,"%.30e" %CTE_alpha_b)
    print ('CTE_alpha_c= ' ,"%.30e" %CTE_alpha_c)
    print ('c_wf_for_cte=' ,"%.30e" %c_wf_for_cte)
    print ('c_wf_for_a0= ' ,"%.30e" %c_wf_for_a0)
    plt.figure(figsize=(8,6))
    for i in range(len(result_master_time)):
        plt.plot(result_master_temp[i][10:-10],1e9*result_master_lattice_param_gamma_during_trans[i][10:-10], '-',label=file_list[i][0])
    plt.legend(loc=4)
    plt.xlabel('Temperature $(\\degree C)$')
    plt.ylabel('Austenite lattice parameter $(nm)$')
    plt.grid()      
    
    plt.figure(figsize=(8,6))
    for i in range(len(result_master_time)):
        plt.plot(result_master_C_in_gama[i][10:-10],1E9*result_master_lattice_param_gamma_during_trans[i][10:-10], '-',label=file_list[i][0])
    plt.legend(loc=4)
    plt.xlabel('C in austenite(mole fraction)')
    plt.ylabel('Austenite lattice parameter$(nm)$')
    plt.grid()     
    
    
    current_script_name= os.path.basename(__file__)
    L0_correction_report=""
    if correct_L0.lower()=='yes':
        if L0_correction_method==3:
            L0_correction_report="yes with method %d, normalizing temp=%d" %(L0_correction_method,normalizing_temp)
        else:
            L0_correction_report="yes with method %d" %L0_correction_method
    else:
        L0_correction_report=correct_L0
    if use_avr_CTE_product==0:
        initial_CTE_alpha_quess="no"
    else:
        initial_CTE_alpha_quess='yes'
    #this function saves the raw data and results as an excel file. 
    run_parameters={"sr":sr,"reg_order":reg_order,"interval":interval,"end_fit_WF":end_fit_WF,
    'Bs_model':Bs_model, "Ms_model":Ms_model, "overal_fit_WF":overal_fit_WF, "err_end_slope_WF":err_end_slop_WF,
    "err_maximum_transformation_WF":err_maximum_transformation_WF,"XRD for a0_alpha":use_XRD_for_a0_alpha,
    "a0_gama":a0_gama,"a0_alpha":a0_alpha, "CTE_alpha_a":CTE_alpha_a,
    "CTE_alpha_b":CTE_alpha_b, "CTE_alpha_c":CTE_alpha_c,"c_wf_for_cte":c_wf_for_cte,
    "c_wf_for_a0":c_wf_for_a0 , "script_name":current_script_name,'optimized':optimize,
    "L0_Correction":L0_correction_report, "use avr CTE product for initial quess of CTE_alpha":initial_CTE_alpha_quess}
    
    file_location=saver(file_list,time_master, temp_master, dil_master, result_master_time,result_master_temp, result_master_dil, result_master_fraction,result_master_ID,result_master_FD_fraction, result_master_C_in_FD_fraction, result_master_C_in_FD_WP,result_master_cem_fraction, result_master_Vol_f_transformed,result_master_FD_vol, result_master_cem_vol,result_master_aus_vol, run_parameters)
    end=time.time()
    if open_excel_result.lower()=='yes':
        call(['open', file_location])
    print ('run time=',end-start)
    if sys.platform[:3]==("lin" or "dar"):
        os.system("spd-say finished")
    #cProfile.run('main()')
    #main()
    #engine = pyttsx.init()
    #engine.say('program has finished running.')
    #engine.runAndWait()

    if show_plots.lower()=='yes':
        print ("plots are chosen to be shown")
        plt.show()

        
            

