#import numpy as np
cimport numpy as np
from numpy.polynomial import polynomial as P

#def L(np.ndarray[np.float_t, ndim=1] C):
#    return C/(1-C)

def VCement(float T): #T in centigrade.
#    cdef float y 
#    y=1+(5.311E-6-1.942E-9*T+9.655E-12*T*T)*(T-20)
    return 0.45234E-9*0.50883E-9*0.67426E-9*(1+(5.311E-6-1.942E-9*T+9.655E-12*T*T)*(T-20))**3

def Valpha(float a0_alpha, float C, float T, float CTE_alpha_a,float CTE_alpha_b, float CTE_alpha_c, float c_wf_for_cte, float c_wf_for_a0, np.ndarray[np.float_t, ndim=1] A): #C is in mole fraction, T is in celsius       
	#C= C_in_alpha(T) #This generates C in fole fraction. The next line will change it to weight percent     
    cdef float CTE_alpha_CC, y, A_alpha, C_alpha, Valpha 
    CTE_alpha_CC=CTE_alpha_c-C*c_wf_for_cte
#    C = 13.556037046*C*C*C + 16.897583563*C*C + 21.525922793*C - 0.0000000034672863191   
    C=P.polyval(C,A)
    y = CTE_alpha_a/3*(T*T*T-8000)+CTE_alpha_b/2*(T*T-400)+CTE_alpha_CC*(T-20) # this is the integration of the CTE equation. 
    if C<0.550:
        A_alpha=(a0_alpha+c_wf_for_a0*C)*(1+y)
        C_alpha=A_alpha
    else:
        A_alpha=(a0_alpha-0.013e-10*C)*(1+y)
        C_alpha=(a0_alpha+0.116e-10*C)*(1+y)
    return A_alpha*A_alpha*C_alpha

def Vbainite(float a0_alpha, float C, float T, float CTE_alpha_a,float CTE_alpha_b, float CTE_alpha_c, float c_wf_for_cte,float c_wf_for_a0, np.ndarray[np.float_t, ndim=1] A): #C is in mole fraction, T is in celsius       
	#C= C_in_alpha(T) #This generates C in fole fraction. The next line will change it to weight percent     
    cdef float CTE_alpha_CC, y, A_alpha, C_alpha, Valpha 
    CTE_alpha_CC=CTE_alpha_c-C*c_wf_for_cte
#    C = 13.556037046*C*C*C + 16.897583563*C*C + 21.525922793*C - 0.0000000034672863191   
    C=P.polyval(C,A)
    y = CTE_alpha_a/3*(T*T*T-8000)+CTE_alpha_b/2*(T*T-400)+CTE_alpha_CC*(T-20) # this is the integration of the CTE equation. 
    if C<0.55:
        A_alpha=(a0_alpha+c_wf_for_a0*C)*(1+y)
        C_alpha=A_alpha
    else:
        A_alpha=(a0_alpha-0.013e-10*C)*(1+y)
        C_alpha=(a0_alpha+0.116e-10*C)*(1+y)
    return A_alpha*A_alpha*C_alpha

