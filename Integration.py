import sys
sys.path.append('BackEnd/')
import PypyNumpy as np
import numpy as npp
import math
import Inputs as me
import Functions as funct 
NormPDF = funct.NormPDF
trapazoid = funct.trapazoid
Delta_Number = funct.Delta_Number
R0 = me.R0
        
Number = int(me.Number)                                                    # Number of neutrons to simulate
# 1 Produced radius shells
Neutron_Radius = npp.linspace(0,me.R0,me.NumberShells)
Neutron_Time = npp.linspace(-(me.TotalBurnSigma/2)*me.BurnSigma,me.TotalBurnSigma/2*me.BurnSigma,math.floor(me.BurnSigma*me.TotalBurnSigma/me.TimeStep))
   
def trapazoid_2_d(function,x_1,x_2,*other):
    I = 0
    del_x = (x_1[-1]-x_2[0])/(2*(len(x_1)-1))
    A = 2*npp.identity(len(x_1))
    A[0][0] , A[-1][-1] = 1 , 1    
    for t in x_2:
        y = function(np.vector(list(x_1)),t,*other)
        I += sum(A.dot(y))*del_x
    return I*((x_2[-1]-x_2[0])/len(x_2))
    
###################################################     
number_per_time ,number_per_r = [],[]
t_delta = Neutron_Time[1]-Neutron_Time[0]

''' Normalize the integral '''
Total = trapazoid_2_d(Delta_Number,Neutron_Radius,Neutron_Time,0)

def NORM(radius,time):
    return Delta_Number(radius,time,0)*Number/Total
'''Figure out how many are in each time '''

for time in Neutron_Time:
    them = trapazoid_2_d(NORM,Neutron_Radius,[time,time+t_delta])
    number_per_time.append(them)   
''' Now for that time, figure out how many are at each radius '''
Neutron_Radius = np.vector(list(Neutron_Radius))
Neutron_Time = np.vector(list(Neutron_Time))
index = 0

for number in number_per_time:
    Normeder = NormPDF(Neutron_Radius,Delta_Number,number,Neutron_Time[index],0)
    Number_In = trapazoid(Neutron_Radius,Normeder,'Discrete')
    number_per_r.append(Number_In)
    index+=1
f = open('BirthLocation.txt','w')
for i in number_per_r:
    for j in i:
        if j == i[-1]:          # If this is the last element in the row vector
            f.write('%f\n' % float(j)) # Make a new line after this element 
        else:
            f.write('%f ' % float(j))
f.close()
