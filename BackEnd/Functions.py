import math
from bisect import bisect_left
from random import random as ran
from Inputs import *
import PypyNumpy as np

'''-----------------------------------'''
''' Specific Functions for simulation '''
'''-----------------------------------'''

#################################################
############ Neutron Birth PFD's ################
#################################################
def ThetaPDF(x):
    ''' The neutron creation is assumed isotropic ''' 
    return np.vector(ones(N_theta))
def PhiPDF(x):
    ''' The neutron creation is assumed isotropic '''
    return np.vector(ones(N_theta))
skew , kurt = 0 , 0
def ballabio_shift(T_ion,reaction):
    ''' Inputs: T_ion [keV], Fusion Reaction [#]; Outputs: Mean energy shift [MeV] '''
    return alpha[reaction][0][0]*T_ion**(2/3)/(1+alpha[reaction][0][1]*T_ion**alpha[reaction][0][2])+alpha[reaction][0][3]*T_ion
def HatarikInfor(T_ion,mean,reaction):
    ''' Inputs: T_ion [keV], mean momentum [MeV], Fusion Reaction [#]; Outputs: width of neutron spectrum [MeV] '''
    return 2*Mn*(math.sqrt(mean[reaction]**2+Mn**2)-Mn)*(mean[reaction]**2+Mn**2)*(T_ion)/(mean[reaction]**2*(Mn+MHe[reaction]))    
def HatarikPDF(p,mue,sigma,skew,kurt):
    ''' Inputs: Momentum [MeV], Momentums of distribution; Outputs: Probability of neutron production at inputed momentum '''
    x = (p-mue)/sigma
    return ((2*math.pi*sigma)**(-1/2))*((-1*x**2/2).exp())*(Hermite(0,x)+skew/6*Hermite(3,x)+kurt/24*Hermite(4,x))
###################################
########### Reactivity ############
###################################
def Temp(r,time):
    ''' Inputs: Radius [m], Time [s]; Outputs: Temperature of the plasma [keV] '''
    ''' The .2 is the floor temperature that we allow the gas to be '''
    T = T_max*(1-(r/(1.0001*R0))**2)**(2/7)*math.exp(-.5*((time/BurnSigma))**2)
    if isinstance(r,type(np.vector(0))) == True:
        T = list(T)
        for i in range(len(T)):
            if T[i] < .2:            
                T[i] = .2 
        return np.vector(T)
    else:
        if T < .2:
            T = .2
        return T
def Delta_Number(radius,time,reaction):
    ''' Inputs: Radius [m], Time [s], Reaction [#]; Outputs: dN produced at inputted Radius and Temp '''
    constant = 1
    if reaction == 1 or reaction == 2: # if a DD or TT fusion reaction occurs the weighting goes down by a factor of 2 because of the 1/(1+delta)
        constant = .5
    return constant*(reactivity(Temp(radius,time),reaction)*N_i_frac[0][0]*N_i_frac[0][1]*4*math.pi*radius**2)
def Theta(T,reaction):
    return T/(1-((T*(c[reaction][1]+T*(c[reaction][3]+T*c[reaction][5])))/(1+T*(c[reaction][2]+T*(c[reaction][4]+T*c[reaction][6])))))
def Espi(T,reaction):
    return (B_g[reaction]**2/4/Theta(T,reaction))**(1/3)
def reactivity(T,reaction):
    return Theta(T,reaction)*c[reaction][0]*(Espi(T,reaction)/T**3/m_r[reaction])**(1/2)*(-3*Espi(T,reaction)).exp()
################################
### Birth Neutron Specs ######## 
################################  
def Emission():
    ''' Determines angle of emitted neutron '''
    def ThetaPick(): 
        return Choice(CDF_theta,x_theta)
    def PhiPick():
        return Choice(CDF_phi,x_phi)
    theta_momentum , phi_momentum = ThetaPick() ,PhiPick()  # The neutrons inital velocity is isotropically outwards
    return  theta_momentum,phi_momentum,theta_momentum,phi_momentum
################################
### Implision Velocity  ######## 
################################ 
def ImplosionVelocity(r,sigma):
    if isinstance(r,int)==True or isinstance(r,float)==True:
        return 2/math.sqrt(math.pi)*math.exp(-(r)**2/2/sigma**2)
    else:
        return (-1*(r**2)/2/sigma**2).exp()

'''--------------------------------------'''
''' General Functions used in simulation '''
'''--------------------------------------'''

####################################
##### Energy <-> Velocity ##########
#################################### 
def EnergyToVelocity(energy,m):
    ''' energy is energy of neutron in eV , mass in MeV'''
    return 2.999999E8*(math.sqrt((energy**2)+2*energy*m*1E6)/(energy+m*1E6)) # converts velocity to m/s
####################################
##### Momentum <-> Energy ##########
####################################    
def MomentumToEnergy(p,m):
    ''' momentum is in MeV, and m is in MeV '''
    return math.sqrt(m**2+p**2)-m    
def EnergyToMomentum(E,m):
    ''' momentum is in MeV, and m is in MeV '''
    return math.sqrt((E+m)**2-m**2)   
###################################
### General Purpose Functions #####
###################################
def linspace(start,stop,n):
    assert n > 1, 'N must be greater than one '''
    h = (stop-start)/(n-1)
    return [start+h*i for i in range(n)]
def ones(n):
    assert n > 0, 'N must be greater than zero '''
    return [1 for i in range(n)]     
def trapazoid(x,y,option):
    ''' Inputs: x data [type=list], y data [type = list], option [type=string] '''
    ''' Depending on what option is set to you will get different inputs '''
    ''' Option == 'Sum' returns to total integral sum
        Option == 'Cumulative' returns a list of length x containing the cumulative integral where the first element is zero and the last element is equal to the total integral 
        Option == 'Discrete' returns a list of length x contaning the integral for each porition of the dicrete intgral such that the sum of the list is the total integral '''
    if option == 'Sum':
        return (x[-1]-x[0])/(2*(len(x)-1))*(y[0]+2*sum(y[1:-1])+y[-1])
    index ,count,done= 0,0,[0]
    def area_trap(a,b,y_1,y_2):
        return (b-a)/2*(y_1+y_2)
    while index < len(x)-1:
        if option == 'Cumulative':
            count += area_trap(x[index],x[index+1],y[index],y[index+1])
        if option == 'Discrete':
            count = area_trap(x[index],x[index+1],y[index],y[index+1])
        done.append(count)
        index += 1
    return done
def interpol(x,x_data,y_data,option='quad'):
    ''' Inputs: x [type=float or int], x_data [type=list], y_data [type=list], option [type=string]; Output: Scalar interpolation approximation '''
    ''' option allows you to change the interpolation method but a quadratic polynomails is the default '''
    def fit(x,x_list,y_list):
        if option == 'quad':
            return y_list[0]*(x-x_list[1])*(x-x_list[2])/(x_list[0]-x_list[1])/(x_list[0]-x_list[2]) + y_list[1]*(x-x_list[0])*(x-x_list[2])/(x_list[1]-x_list[0])/(x_list[1]-x_list[2]) + y_list[2]*(x-x_list[0])*(x-x_list[1])/(x_list[2]-x_list[0])/(x_list[2]-x_list[1])
        elif option == 'linear':
            m = (y_list[1]-y_list[0])/(x_list[1]-x_list[0])
            return m*(x-x_list[0]) + y_list[0]
        elif option == 'log':
            y_list = [math.log(i) for i in y_list]
            m = (y_list[1]-y_list[0])/(x_list[1]-x_list[0])
            return m*(x-x_list[0]) + y_list[0]    
    pos = bisect_left(x_data,x)
    known = [pos-1,pos,pos+1]
    if pos == 0:
        known = [0,1,2]
    if pos == len(x_data)-1 or pos == len(x_data):
        known = [len(x_data)-3,len(x_data)-2,len(x_data)-1]
    if len(x_data) == 2:
        option = 'linear'
        return fit(x,[x_data[known[0]],x_data[known[1]]],[y_data[known[0]],y_data[known[1]]]) 
    else:
        return fit(x,[x_data[known[0]],x_data[known[1]],x_data[known[2]]],[y_data[known[0]],y_data[known[1]],y_data[known[2]]]) 
def interpolate(x_data,y_data,N,option='quad'):
    x_new = linspace(min(x_data),max(x_data),N)
    y_new = [interpol(i,x_data,y_data,option) for i in x_new]
    return x_new , y_new
###################################
##### Hermite Polynomials##########
###################################
''' First creates a table of each coefficant of the hermite polynomails using a recursion formula '''
order = 6         # What order polynomial do you need up to?
coef = [(order+1)*[0] for i in range(order)]
coef[1][0] , coef[2][1] = 1 , 1
for n in range(order-1):
    for counter in reversed(range(1,order-1)):
        coef[n+1][counter+1] += (coef[n][counter])
for n in range(order-1):
    for sort in reversed(range(1,order)):
        if coef[n][sort] != 0:
            break
    for counter in range(order):
        if counter < sort:
            coef[n+1][counter] += coef[n][counter-1]-(n-1)*coef[n-1][counter]
        else:
            break
def Hermite(n,x):
    H , n = 0 , n+1  
    for m in range(n):
        H += coef[n][m]*x**m
    return H   
###################################
######### Geometry ################
###################################
def Geometry(position):
    ''' Inputs: Position [type=Vector of length 3]; Output: Index of geometric region [type=int] '''    
    radius = math.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
    return bisect_left(Radius,radius)
#    def checker(x,check):
#        for i in check:
#            if i < x:
#                checker(x,check[:-check.index(i)])
#            else:
#                return check.index(i)
#        return len(check)-1
#    return checker(radius,Radius)
def GeometricDistance(position,velocity,region):
    if region == 0:
        return GeometricExpansionSphere(position,velocity,Radius[region])
    root1 = GeometricExpansionSphere(position,velocity,Radius[region])
    root2 = GeometricExpansionSphere(position,velocity,Radius[region-1])
    if type(root1) == str:
        if type(root2) == str:
            return 'Negative'
        else:
            return root2
    if type(root2) == str:
        if type(root1) == str:
            return 'Negative'
        else:
            return root1
    elif root1 < root2:
        return root1
    else: return root2   
def GeometricExpansionSphere(position,velocity,radius):
    ''' Inputs: position [type vector of length 3], velocity [type=vector of length 3], radius of spherical shell[type=float];
        Outputs: the minimum distance to the geometrical boudnary'''
    vel_unit = velocity/velocity.norm()
    coef_1 = vel_unit[0]**2 + vel_unit[1]**2 + vel_unit[2]**2
    coef_2 = 2 * (position[0]*vel_unit[0]+position[1]*vel_unit[1]+position[2]*vel_unit[2])
    coef_3 = position[0]**2 + position[1]**2 + position[2]**2- radius **2
    if coef_2**2 - 4*coef_1*coef_3 < 0:
        root1 , root2 = 'Imaginary' , 'Imaginary'
    else:
        root1  = (-coef_2 + math.sqrt(coef_2**2 - 4*coef_1*coef_3)) / 2*coef_1
        root2  = (-coef_2 - math.sqrt(coef_2**2 - 4*coef_1*coef_3)) / 2*coef_1   
    ''' Ensures that the vector is on the boundary, but round off errors make it negative, that we choose it anyways '''
    if root1 == 'Imaginary':
        return 'Imaginary'
    if abs(root1) < 1E-10:
        if abs(root2) < 1E-10:
            if abs(root1) < abs(root2):
                return root1
            else:
                return root2
        else:
            return abs(root1)
    else:       
        if root1 < 0:
            if root2 < 0:
                return 'Negative'
            else:
                return root2
        elif root2 < 0:
            if root1 < 0 :
                return 'Negative'
            else:
                return root1
        elif root1 < root2:
            return root1
        else: return root2  
###################################
### Creating CDF's ################
###################################
def Normalize(x,pdf,normalized_to,*other):
    ''' Inputs: x [type=list or vector], pdf [type = function or list or vector], normalized_to [type=float or int], *other are extra parameters that are needed for the pdf fuction used
        Output: a constant in which to multiply the pdf to in order to normalize it to the specificed normalization '''
    if callable(pdf) == True:
        y = pdf(x,*other)
    else:
        y = pdf
    summation = trapazoid(x,y,'Sum')
    if normalized_to == 0 or summation == 0:
        return 0
    else:
        return normalized_to/trapazoid(x,y,'Sum')   
def NormPDF(x,pdf,normalized_to,*other):
    ''' Inputs: x [type=list or vector], pdf [type = function or list or vector], normalized_to [type=float or int], *other are extra parameters that are needed for the pdf fuction used
        Output: a normalized pdf [type = function or vector] '''
    cont = Normalize(x,pdf,normalized_to,*other)
    if callable(pdf) == True:
        return pdf(x,*other)*cont  
    else:
        if isinstance(pdf,type(np.vector([0]))) != True:
            pdf = np.vector(pdf)
        return pdf*cont
def CreateCDF(pdf,a,b,N,normalized_to,*other):
    ''' Generalized CDF Creation from Normalized continious PDF'''
    ''' Inputs: pdf [type = function or vector], a = beginign of domain[type = float or int], b end of domain = [type = float or int], N = [type = int], normalized_to [type = float or int], 
        *other are extra parameters that are needed for the pdf fuction used. 
        Outputs: a list containing the cdf values '''
    if callable(pdf) == True :
        x = np.vector((linspace(a,b,N))) 
        return trapazoid(x,NormPDF(x,pdf,normalized_to,*other),'Cumulative'),x ##############################
    else:
        x = linspace(a,b,len(pdf)-1) 
        return trapazoid(x,NormPDF(x,pdf,normalized_to),'Cumulative'),x    
###################################
### Choosing from CDF #############
###################################
def Choice(cdf,x):
    ''' Inputs: cdf [type = list or vector], x values [type=list or vector]; Output: floating point value of choice '''
    ''' Returns the the chosen value (corresponding to the CDF) that are closest to the choice random value (ranging from 0 to 1)'''
    choice = ran()
    return interpol(choice,cdf,x,option='linear')

###################################
###### Choice Functions ###########
###################################
def DiscreteChoice(events):
    ''' Events is a list of the form [0,...information...,0]. So that information is the information about the weighting of each event you consider '''
    ep = ran()
    middle = ep*sum(events)
    index = 1
    for index in range(len(events)):    
        left  = sum(events[0:index])
        right = left + events[index] #sum(events[0:index+1])
        if (left < middle <= right) == True:
            return(index-1)
##################################
        
N_theta ,N_energy = 100 , 100          # Number of points to model PDF functions
theta_min , theta_max = 0 , 2*math.pi   # Birth theta limits
phi_min , phi_max = 0,math.pi           # Birth phi limits
CDF_theta, x_theta = CreateCDF(ThetaPDF,theta_min,theta_max,N_theta,1)
CDF_phi , x_phi = CreateCDF(PhiPDF,phi_min,phi_max,N_theta,1)
ImplosionVelocity,ImplosionRadius = CreateCDF(ImplosionVelocity,0,MaxRadius,100,PeakVelocity,Radius[2]/(2*math.sqrt(2*math.log(2))))
 
