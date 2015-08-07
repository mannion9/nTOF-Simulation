import sys
sys.path.append('BackEnd/')
import math,subprocess
import PypyNumpy as np
from time import time as tic
from math import cos as cos
from math import sin as sin
from random import random as ran
start = tic() 
'''This block of code imports functions and inputs from the child code'''
from Inputs import *
import Inputs as me
from Functions import *


###################################
##### Collisions ##################
###################################
def CrossSection(energy,region):
    ''' Energy is neutron energy in eV '''
    ''' Decides which atom to collide with, records it, calculates that cross section
       and returns it.'''
    cross_section = list(np.vector(Fits(energy,region))*Ratio[region])  # this line takes into account the % of substance in the material's effect on the choice of what atom is collided with
    cross_section.append(0),cross_section.insert(0,0)
    choice = DiscreteChoice(cross_section)
    xc = cross_section[choice+1]*BarnToCmSquare
    return(xc,choice) 
    
def ScatteringAngle(energy,choice,region):
    ''' Takes in the neutron energy and returns the scattering angle in the CM frame '''
    theta = DataBaseChoice(energy,choice,region) 
    phi = Choice(CDF_theta,x_theta)
    alpha = ((A[region][choice]-1)/(A[region][choice]+1))**2
    energy = .5*energy*((1+alpha)+(1-alpha)*theta)
    return theta,phi,energy
    
def DataBaseChoice(energy,choice,region):
    for i in range(len(Energy_dic_atoms)):
        if A[region][choice] == Energy_dic_atoms[i]:
            choice_index = i
    index_energy = min(range(len(Energy_dic[choice_index])),key=lambda i: abs(Energy_dic[choice_index][i]-energy))
    if Energy_dic[choice_index][index_energy] not in Energy_current[choice_index]:
        Energy_current[choice_index].append(Energy_dic[choice_index][index_energy])
        cosine_current = Cosines[choice_index][index_energy]
        probs_current = Probability[choice_index][index_energy]
        if A[region][choice] != 2 and A[region][choice] != 14:
            cosine_current.reverse()
            probs_current.reverse()
        cosine , probs= interpolate(cosine_current,probs_current,1000)         
        cdf , cosine = CreateCDF(probs,cosine[0],cosine[-1],1000,1)
        Angles[choice_index].append(cosine)
        Probs[choice_index].append(probs)
        CDF[choice_index].append(cdf)
    current_index = Energy_current[choice_index].index(Energy_dic[choice_index][index_energy])
    angle = Choice(CDF[choice_index][current_index],Angles[choice_index][current_index])
    return angle       
        

cross_sections_t = 0
###################################
##### Cross Section fits ##########
###################################
def Fits(energy,region):
    global cross_sections_t
    starting = tic()
    # Fits returing values
    cross_sections = []
    for i in A[region]:
        index = Energy_dic_atoms.index(i)
        cross_sections.append(interpol(energy,Cross_Energy[index],Cross[index],option='linear'))
    cross_sections_t += tic()-starting
    return cross_sections
   

Angles =[ [] for i in range(len(Energy_dic_atoms))]
Probs = [ [] for i in range(len(Energy_dic_atoms))]
CDF =[ [] for i in range(len(Energy_dic_atoms))]
Energy_current = [ [] for i in range(len(Energy_dic_atoms))]

    

COUNT = 0   
def Transport(position,energy,theta_momentum,phi_momentum,region,collision_counter):
    global COUNT
    steps = 0 
    velocity = EnergyToVelocity(energy,Mn)
    region = Geometry(position)
    velocity_direction = np.vector([cos(theta_momentum)*sin(phi_momentum),sin(theta_momentum)*sin(phi_momentum),cos(phi_momentum)])
    xc , atom = CrossSection(energy,region)
    mean_free_path = (1/NDensity[region]/xc)/100  # converts cross section into mean free path in meters
    geometric_distance = GeometricExpansionSphere(position,velocity_direction*velocity,Radius[region])
    scatter_atom = 0
    if geometric_distance <= 10E-9:
        ''' this is necessary because the GeometricExpansionSphere function find the roots of a polynomial, and 
            finding roots numerically is always difficult and leads to some amount of round off errors, this check
            is to ensure that if you are almost at a boundry there is a good reason to believe that if you are not
            excactly on the boundry, then you are supposed to be on the boundry but round off errors caused you to 
            stay just below the boundry, check accounts for this.'''
        region += 1
        geometric_distance = GeometricExpansionSphere(position,velocity_direction*velocity,Radius[region])       
        mean_free_path = (1/NDensity[region]/xc)/100
    if mean_free_path < geometric_distance:   
        COUNT += 1
        position = position+velocity_direction*mean_free_path
        tau = mean_free_path / velocity
        collision_counter+=1
        cos_theta_CMF ,phi, energy = ScatteringAngle(energy,atom,region)
        if abs(cos_theta_CMF) > 1:
            if cos_theta_CMF < 0:
                cos_theta_CMF = -1
            else:
                cos_theta_CMF = 1
        cos_theta_lab = (1 + A[region][atom]*cos_theta_CMF)/math.sqrt((A[region][atom]**2)+1+(2*A[region][atom]*cos_theta_CMF))
        theta_lab = math.acos(cos_theta_lab)
        theta_momentum += theta_lab
        scatter_atom = A[region][atom]
        steps += 1
    else:
        tau = 0
        while position.norm() < Radius[region]:
            steps += 1
            eps = ran()
            if region == 0:
                travel = eps*Radius[region]
            else:    
                travel = eps*(Radius[region]-Radius[region-1])  
            position += velocity_direction*travel
            if position.norm() > Radius[-2]:  
                travel -= position.norm()-Radius[-2]
                position = np.vector([Radius[-2],0,0]) # forces the neutron to the edge perfectly (I was finding many getting stuck just before the edge)
            tau += travel/EnergyToVelocity(energy,Mn)
            p_scatter = (Radius[region]/mean_free_path)*(1-position.norm()/mean_free_path)
            collision_prob = [0,p_scatter,1-p_scatter,0]
            choice = DiscreteChoice(collision_prob)
            if choice == 0:
                collision_counter += 1
                cos_theta_CMF ,phi, energy = ScatteringAngle(energy,atom,region)
                if abs(cos_theta_CMF) > 1:
                    if cos_theta_CMF < 0:
                        cos_theta_CMF = -1
                    else:
                        cos_theta_CMF = 1   
                cos_theta_lab = (1 + A[region][atom]*cos_theta_CMF)/math.sqrt((A[region][atom]**2)+1+(2*A[region][atom]*cos_theta_CMF))
                theta_lab = math.acos(cos_theta_lab)
                theta_momentum += theta_lab
                phi_momentum += phi
                scatter_atom = A[region][atom]
    return position , tau , region ,energy, collision_count,steps,collision_counter,theta_momentum,phi_momentum,velocity_direction,scatter_atom
    
    
    
####################################
#### Create neutron at radius ######
####################################
Number = int(Number)                                                    # Number of neutrons to simulate
# 1 Produced radius shells
Neutron_Radius = np.vector(linspace(0,me.R0,me.NumberShells))             # Radius to create neutrons at
Neutron_Time   = np.vector(linspace(-(TotalBurnSigma/2)*BurnSigma,TotalBurnSigma/2*BurnSigma,math.floor(BurnSigma*TotalBurnSigma/TimeStep)))   # Times to create neutrons at
N_tot = me.M_hs*me.N_A/(Ratio[0][0]*(me.MH2-me.MH3)+me.MH3)                              # Number Density of fuel


on = tic()
# Go To child process to integrate dN/dr/dt 
subprocess.call(['python3','Integration.py'])
# Open file the child process writes out for us
f = open('BirthLocation.txt','r')
f = f.read().splitlines()
number_per_r = []
TOTAL = 0
for line in f:
    row = []
    line = line.split(' ')      # Seperates each number
    for i in line:
        row.append(float(i))
    number_per_r.append(row)
    TOTAL += sum(row)
off = tic()
print('Time in Child_Integration:',tic()-on)
print('Estimated Number of neutrons:',TOTAL)


collisions = []
Timers = []
Stepers = 0
repeat=0
collision_count = 0
OUTPUT = []
Fusions = [0,0,0]  # Records the number of fusions per each reaction (indexing---> 0=DT,1=DD,2=TT)



trans = 0
birthing = 0
CDF_Creation = 0 

##################################
####### Transport neutrons #######
##################################
TRANSPORT = tic()
for time in number_per_r:   
    time_now = Neutron_Time[number_per_r.index(time)]
    for radius in time:
        radius_now = Neutron_Radius[time.index(radius)]
        neutron_count = 0
        T_ion = Temp(radius_now,time_now)
        
        CDF_t = tic()
        mean = [P_zero_temp[0],P_zero_temp[1]]

        ''' added '''
        mean = [mean[i] + EnergyToMomentum(ballabio_shift(T_ion,i),Mn)/1000 for i in range(len(mean))]
            
        sigma = [HatarikInfor(T_ion*1E-3,0),HatarikInfor(T_ion*1E-3,1)] # Tion from keV to MeV
        CDF_0 , x_0 = CreateCDF(HatarikPDF,mean[0]-5*sigma[0],mean[0]+5*sigma[0],N_energy,1,mean[0],sigma[0],skew,kurt) 
        CDF_1 , x_1 = CreateCDF(HatarikPDF,mean[1]-5*sigma[1],mean[1]+5*sigma[1],N_energy,1,mean[1],sigma[1],skew,kurt)
        CDF_2 , x_2 = CreateCDF(tt2n_spectrum[1],tt2n_spectrum[0][0],tt2n_spectrum[0][-1],100,1)
        CDF_Creation += tic()-CDF_t
        
       # TT2n
        index_temp = min(range(len(tt2n_react[0])),key=lambda i: abs(tt2n_react[0][i]-T_ion))
        tt2n_reactivity = tt2n_react[1][index_temp]
        
        reaction_list = [0,1,N_i_frac[0][0]*reactivity((np.vector([T_ion])),1)[0]/2/N_i_frac[0][1]/reactivity(np.vector([T_ion]),0)[0],N_i_frac[0][1]*tt2n_reactivity/2/N_i_frac[0][1]/reactivity(np.vector([T_ion]),0)[0],0]
        
        while neutron_count <= math.floor(radius):
            starterer = tic()
          
            neutron_count += 1
            if repeat == 1:
                Energy = Choice(CDF_2,x_2)*1E6  # eV
                repeat = 0 
            else:
                reaction = DiscreteChoice(reaction_list)  
                if reaction == 2:
                    CDF_energy, x_energy =  CDF_2 , x_2
                    Energy = Choice(CDF_energy,x_energy)*1E6  # eV
                    repeat = 1
                else:
                    if reaction == 0: 
                        CDF_momentum,x_momentum = CDF_0 , x_0 
                    if reaction == 1:
                        CDF_momentum,x_momentum = CDF_1 , x_1
                    Energy = MomentumToEnergy(Choice(CDF_momentum,x_momentum),Mn)*1E6 # eV
                    
            birthing += tic()-starterer
            
            Fusions[reaction] += 1
            BirthEnergy = Energy/1E6
            ThetaEmission , PhiEmission, ThetaMomentum, PhiMomentum = Emission()
            Theta_momentum = ThetaEmission
            Phi_momentum = PhiEmission
            Position = np.vector([cos(ThetaEmission)*sin(PhiEmission),sin(ThetaEmission)*sin(PhiEmission),cos(PhiEmission)])*radius_now
            Region = Geometry(Position)
            Timer = 0
            collisions_counter = 0
            scatters = []  
            
            tran_start = tic()
    
            while Position.norm() < Radius[-2]:
                collision_pre = collisions_counter
                Position , Tau , Region , Energy ,collision_count,steps,collisions_counter,Theta_momentum,Phi_momentum,Velocity_direction,scatter_atom = Transport(Position,Energy,Theta_momentum,Phi_momentum,Region,collisions_counter)                
                scatters.append(scatter_atom)                
                if collisions_counter > collision_pre:
                     collisions.append(1)
                Timer += Tau
                Stepers += steps 
                
            trans += tic()-tran_start    
            
            scatters.reverse()
            last_atom = 0
            for i in scatters:
                if i != 0:
                    last_atom = i

            Timers.append(Timer/1E-9)
            Velocity_direction = list(Velocity_direction)
            Output = [reaction,BirthEnergy,collisions_counter,int(last_atom),Energy/1E6,EnergyToVelocity(Energy,Mn),Timer/1E-9,Velocity_direction[0],Velocity_direction[1],Velocity_direction[2]]
            OUTPUT.append(Output)
print('In transport part of code:',tic()-TRANSPORT)

print('Pick Cross Section:',cross_sections_t)
print('Time Spent Choosing Fusion type:',birthing)
print('Time Spent Transporting:',trans)
print('Number of Mean Free Paths:',COUNT)

#################################
####### Write out results #######
#################################
f = open('AllData.txt','w')
# NumberCollisions Energy(MeV) Velocity(m/s), Time(ns), VelocityUnitVector
f.write('Reaction meaning, 0=DT,1=DD,2=TT\n')
f.write('Reaction, BirthEnergy (MeV), NumberOfCollisions, AtomicMassOfLastScatter, Energy(MeV), Velocity(m/s), Time(ns), X VelocityUnitVectorComponent, Y VelocityUnitVectorComponent, Z VelocityUnitVectorComponent\n')
for i in OUTPUT:
    for j in i:
        if j == i[0]:
            f.write('%d ' % j)
        elif j == i[-1]:          # If this is the last element in the row vector
            f.write('%f\n' % j) # Make a new line after this element 
        else:
            f.write('%f ' % j)
f.close()
stop = tic()

print('Time of simulation in seconds is',stop-start,'or',(stop-start)/60,'minutes with',len(Timers),'Neutrons')
print('The average time to hit detectors at 20 (m) is:',np.vector(Timers).average(),'nanoseconds')
print('The percent of neutron that under went a collisions was',sum(collisions)/len(Timers)*100)
print('DTn fusions:',Fusions[0]*100/sum(Fusions),'%   DDn fusions:',Fusions[1]*100/sum(Fusions),'% TT2n fusions:',Fusions[2]*100/sum(Fusions))
print(' ')