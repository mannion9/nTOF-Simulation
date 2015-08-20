import sys
sys.path.append('BackEnd/')   # Givvs Run.py access to files in BackEnd subdirectory
import math,subprocess
import PypyNumpy as np
from time import time as tic
from math import cos as cos
from math import sin as sin
from multiprocessing import Pool
from contextlib import closing

start = tic() 

from Inputs import *
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
    if abs(theta) > 1:
        if theta < 0:
            theta = -1
        else:
            theta = 1 
    phi = Choice(CDF_theta,x_theta)
    alpha = ((A[region][choice]-1)/(A[region][choice]+1))**2
    energy = .5*energy*((1+alpha)+(1-alpha)*theta)
    return theta,phi,energy
databasechoicetime = 0  
def DataBaseChoice(energy,choice,region):
    global databasechoicetime
    start = tic()
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
    databasechoicetime += tic()-start
    return angle       
        

###################################
##### Cross Section fits ##########
###################################
def Fits(energy,region):
    global cross_sections_t
    # Fits returing values
    cross_sections = []
    for i in A[region]:
        index = Energy_dic_atoms.index(i)
        cross_sections.append(interpol(energy,Cross_Energy[index],Cross[index],option='linear'))
    return cross_sections
   

Angles =[ [] for i in range(len(Energy_dic_atoms))]
Probs = [ [] for i in range(len(Energy_dic_atoms))]
CDF =[ [] for i in range(len(Energy_dic_atoms))]
Energy_current = [ [] for i in range(len(Energy_dic_atoms))]

def BeerLaw(x,mean_free_path):
    return   1-(-1*x/mean_free_path).exp()      

COUNT = 0 
count_list = []  

MaxImplosionRadius = max(ImplosionRadius)
    
    
####################################
#### Create neutron at radius ######
####################################
Number = int(Number)                                                    # Number of neutrons to simulate
# 1 Produced radius shells
Neutron_Radius = np.vector(linspace(0,R0,NumberShells))             # Radius to create neutrons at
Neutron_Time   = np.vector(linspace(-(TotalBurnSigma/2)*BurnSigma,TotalBurnSigma/2*BurnSigma,math.floor(BurnSigma*TotalBurnSigma/TimeStep)))   # Times to create neutrons at
N_tot = M_hs*N_A/(Ratio[0][0]*(MH2-MH3)+MH3)                              # Number Density of fuel


### Go To child process to integrate dN/dr/dt
subprocess.call(['python3','Integration.py'])
### Open file the child process writes out for us
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
print('Estimated Number of neutrons:',math.floor(TOTAL))

collisions = 0
Timers = []
Stepers = 0
repeat=0
collision_count = 0
OUTPUT = []
Fusions = [0,0,0]  # Records the number of fusions per each reaction (indexing---> 0=DT,1=DD,2=TT)
RUNNER = 0
OVERAL = 0
''' Neutron Attributes [Energy,Radius,Theta,Phi] '''
NeutronAttributes = [] 
##################################
####### Transport neutrons #######
##################################
for time in number_per_r: 
    time_now = Neutron_Time[number_per_r.index(time)]
    for radius in time:       
        
        radius_now = Neutron_Radius[time.index(radius)]
        neutron_count = 0
        T_ion = Temp(radius_now,time_now)
        mean = [P_zero_temp[0],P_zero_temp[1]]
        if Ballabio == 1:
            mean = [mean[i] + EnergyToMomentum(ballabio_shift(T_ion,i),Mn)/1000 for i in range(len(mean))]
            
        sigma = [HatarikInfor(T_ion*1E-3,mean,0),HatarikInfor(T_ion*1E-3,mean,1)] # Tion from keV to MeV
        CDF_0 , x_0 = CreateCDF(HatarikPDF,mean[0]-5*sigma[0],mean[0]+5*sigma[0],N_energy,1,mean[0],sigma[0],skew,kurt) 
        CDF_1 , x_1 = CreateCDF(HatarikPDF,mean[1]-5*sigma[1],mean[1]+5*sigma[1],N_energy,1,mean[1],sigma[1],skew,kurt)
        CDF_2 , x_2 = CreateCDF(tt2n_spectrum[1],tt2n_spectrum[0][0],tt2n_spectrum[0][-1],100,1)
       # TT2n
        index_temp = min(range(len(tt2n_react[0])),key=lambda i: abs(tt2n_react[0][i]-T_ion))
        tt2n_reactivity = tt2n_react[1][index_temp]
        
        reaction_list = [1,N_i_frac[0][0]*reactivity((np.vector([T_ion])),1)[0]/2/N_i_frac[0][1]/reactivity(np.vector([T_ion]),0)[0],N_i_frac[0][1]*tt2n_reactivity/2/N_i_frac[0][1]/reactivity(np.vector([T_ion]),0)[0]]
        reaction_list = [reaction_list[i]*fusion_list[i] for i in range(len(reaction_list))]
        reaction_list.append(0),reaction_list.insert(0,0)
       
        while neutron_count <= math.floor(radius):
            RUNNER += 1
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
            
            Fusions[reaction] += 1
            ThetaEmission , PhiEmission, Theta_momentum, Phi_momentum = Emission()

            NeutronAttributes.append((Energy,radius_now,ThetaEmission,PhiEmission,Theta_momentum,Phi_momentum,reaction))
            
def mini_data_func(position,velocity,energy,region,weighting):
    exponential = 0    
    ##vector_sum = position+detector_position
    velocity_vs_detector = velocity.dot(detector_position)/velocity.norm()/detector_position.norm()    
    ##velocity_difference = velocity.dot(vector_sum)/velocity.norm()/vector_sum.norm()
    #CM_needed = math.cos(velocity_difference)    
    CM_needed = velocity_vs_detector #math.cos(math.pi - velocity_vs_detector)    
#    print('Velocity_against detector:',math.acos(velocity_vs_detector)*180/math.pi)
#    print('Center of mass:',math.acos(CM_needed)*180/math.pi)
#    
    def find_prob(region):
        Probs = []
        for choice in range(len(A[region])):
            for i in range(len(Energy_dic_atoms)):
                if A[region][choice] == Energy_dic_atoms[i]:
                    #print('Lab scatter needed:',(1 + A[region][choice]*CM_needed)/math.sqrt((A[region][choice]**2)+1+(2*A[region][choice]*CM_needed)))                   
                    choice_index = i
                    index_energy = min(range(len(Energy_dic[choice_index])),key=lambda i: abs(Energy_dic[choice_index][i]-energy))   
                    scatter_index =  min(range(len(Probability[choice_index][index_energy])),key=lambda i: abs(Cosines[choice_index][index_energy][i]-CM_needed))
                    Probs.append(Probability[choice_index][index_energy][scatter_index])
        return Probs
    Probs = find_prob(region)
    
    while position.norm() < Radius[-2]:
        difference = detector_position-position
        region = Geometry(position)
        min_dis = GeometricDistance(position,difference,region)
        if min_dis <= 10E-9:
            ''' this is necessary because the GeometricExpansionSphere function find the roots of a polynomial, and 
                finding roots numerically is always difficult and leads to some amount of round off errors, this check
                is to ensure that if you are almost at a boundry there is a good reason to believe that if you are not
                excactly on the boundry, then you are supposed to be on the boundry but round off errors caused you to 
                stay just below the boundry, check accounts for this.'''
            position += 1E-9*difference/difference.norm()
            region = Geometry(position)
            difference = detector_position-position
            min_dis = abs(GeometricDistance(position,difference,region))            
        position += min_dis*(difference/difference.norm())
        xc = (np.vector(Fits(energy,region))*Ratio[region]).average()*BarnToCmSquare  # Average cross section in cm
        
        exponential += -min_dis*(xc*NDensity[region]/100)          
    weighting *= math.exp(-exponential)*sum(Probs)/len(Probs)
    return weighting

def Transport(position,energy,theta_momentum,phi_momentum,region,collision_counter,steps,reaction,BirthEnergy,Tau):
    def mini_transport(position,energy,theta_momentum,phi_momentum,region,collision_counter,steps,velocity_direction,tau,scatter_atom,weighting,reaction,Tau):
     
         # PUT MINI DATA FUNC HERE
        
        global COUNT  
        tag_transport = 0
        steps += 1
        
        if energy < 50E3 or steps > 50:
            region = 'Error'   
            position = np.vector([Radius[-2],0,0])
            return position , 0 , region ,energy,steps,collision_counter,theta_momentum,phi_momentum,velocity_direction,0,0,weighting,Tau
            
        velocity_direction = np.vector([cos(theta_momentum)*sin(phi_momentum),sin(theta_momentum)*sin(phi_momentum),cos(phi_momentum)])
        Velocity = EnergyToVelocity(energy,Mn)*velocity_direction # vector
        weighting = mini_data_func(position,Velocity,energy,region,weighting)
        xc = (np.vector(Fits(energy,region))*Ratio[region]).average()*BarnToCmSquare
        mean_free_path = (1/NDensity[region]/xc)/100  # converts cross section into mean free path in meters
        
        x = np.vector(linspace(0,5*mean_free_path,100))
        step = Choice(BeerLaw(x,mean_free_path),x)  

        geometric_distance = GeometricDistance(position,Velocity,region)        

        if geometric_distance <= 10E-9:
            ''' this is necessary because the GeometricExpansionSphere function find the roots of a polynomial, and 
                finding roots numerically is always difficult and leads to some amount of round off errors, this check
                is to ensure that if you are almost at a boundry there is a good reason to believe that if you are not
                excactly on the boundry, then you are supposed to be on the boundry but round off errors caused you to 
                stay just below the boundry, check accounts for this.'''
            position += 1E-9*Velocity/Velocity.norm()
            return position , 0 , region ,energy,steps,collision_counter,theta_momentum,phi_momentum,velocity_direction,0,0,weighting,Tau
        if step < geometric_distance:
            COUNT += 1
            collision_counter += 1
            position += step*Velocity/Velocity.norm()
            tau += step/Velocity.norm()
            xc , atom = CrossSection(energy,region)
            mean_free_path = (1/NDensity[region]/xc)/100  # converts cross section into mean free path in meters
            cos_theta_CMF ,phi, energy = ScatteringAngle(energy,atom,region)
            cos_theta_lab = (1 + A[region][atom]*cos_theta_CMF)/math.sqrt((A[region][atom]**2)+1+(2*A[region][atom]*cos_theta_CMF))
            theta_lab = math.acos(cos_theta_lab)
            theta_momentum += theta_lab
            phi_momentum += phi
            scatter_atom = A[region][atom]
        else:
            position += (geometric_distance)*Velocity/Velocity.norm()
            tau += geometric_distance/Velocity.norm()
            tag_transport = 1
        Tau += tau
        return position , tau , region ,energy,steps,collision_counter,theta_momentum,phi_momentum,velocity_direction,scatter_atom,tag_transport,weighting,Tau
    global COUNT 
    tau,weighting = 0,1/math.pi/4
    velocity_direction = np.vector([cos(theta_momentum)*sin(phi_momentum),sin(theta_momentum)*sin(phi_momentum),cos(phi_momentum)])
    Velocity = EnergyToVelocity(energy,Mn)*velocity_direction # vector
    if implosion == 1:
        implosion_velocity = interpol(position.norm(),ImplosionRadius,ImplosionVelocity)
        if steps == 0: # implosion velocity only changes the neutron at birth (aka before any transport so steps==0)
            if position.norm() == 0: # If its at the center
                implosion_direction = -1*velocity_direction
            else:
                implosion_direction = -1*position/radius
            Implosion_Velocity = implosion_velocity*implosion_direction # Vector
            Velocity = (Velocity+Implosion_Velocity)/(1+(Velocity.norm()*implosion_velocity)/3E16) # Add them relativistic
            #Velocity += Implosion_Velocity # Add them classically
            
            Momentum = Mn*(Velocity/3E8)/math.sqrt(1-(Velocity.dot(Velocity))/(3E8)**2)
            energy = MomentumToEnergy(Momentum.norm(),Mn)*1E6

    scatter_atom = 0
    if GeometricDistance(position,Velocity,region) <= 10E-9:
        ''' this is necessary because the GeometricExpansionSphere function find the roots of a polynomial, and 
            finding roots numerically is always difficult and leads to some amount of round off errors, this check
            is to ensure that if you are almost at a boundry there is a good reason to believe that if you are not
            excactly on the boundry, then you are supposed to be on the boundry but round off errors caused you to 
            stay just below the boundry, check accounts for this.'''
        #region += 1
        position += velocity_direction*1E-9
    Mini_Data = []   
    region = Geometry(position)
    if region == len(Radius)-1:
        return position ,0,     region ,energy,steps,collision_counter,theta_momentum,phi_momentum,velocity_direction,  0,  0, Mini_Data
    if region == 0:
       while position.norm() <= Radius[region]:
           position , tau , region ,energy,steps,collision_counter,theta_momentum,phi_momentum,velocity_direction,scatter_atom,tag_transport,weighting,Tau= mini_transport(position,energy,theta_momentum,phi_momentum,region,collision_counter,steps,velocity_direction,tau,scatter_atom,weighting,reaction,Tau)
           # YOUR DATA WRITTEN OUT
           Mini_Data.append([reaction,BirthEnergy/1E6,collision_counter,scatter_atom,energy/1E6,((detector_position.norm()-position.norm())/EnergyToVelocity(energy,Mn)+tau)/1E-9,weighting])           
           if tag_transport == 1 or region == 'Error':
               break
    else:
       while position.norm() <= Radius[region] and position.norm() >= Radius[region-1]:
           position , tau , region ,energy,steps,collision_counter,theta_momentum,phi_momentum,velocity_direction,scatter_atom,tag_transport,weighting,Tau = mini_transport(position,energy,theta_momentum,phi_momentum,region,collision_counter,steps,velocity_direction,tau,scatter_atom,weighting,reaction,Tau)
           # YOUR DATA WRITTEN OUT
           if region == len(Radius)-1:
               Mini_Data.append([reaction,BirthEnergy/1E6,collision_counter,scatter_atom,energy/1E6,Tau,weighting])
           else:
               Mini_Data.append([reaction,BirthEnergy/1E6,collision_counter,scatter_atom,energy/1E6,((detector_position.norm()-position.norm())/EnergyToVelocity(energy,Mn)+tau)/1E-9,weighting]) 
           if tag_transport == 1 or region == 'Error':
               break
    return position , tau , region ,energy,steps,collision_counter,theta_momentum,phi_momentum,velocity_direction,scatter_atom,weighting,Mini_Data
              
def Parallel(NeutronAttributes):
    global collisions,Stepers,Radius
    Energy,radius_now,ThetaEmission,PhiEmission,Theta_momentum,Phi_momentum,reaction = NeutronAttributes[0],NeutronAttributes[1],NeutronAttributes[2],NeutronAttributes[3],NeutronAttributes[4],NeutronAttributes[5],NeutronAttributes[6]
    Position = np.vector([cos(ThetaEmission)*sin(PhiEmission),sin(ThetaEmission)*sin(PhiEmission),cos(PhiEmission)])*radius_now
    BirthEnergy = Energy
    Region = 0 # For this simulation we assume all neutrons are born in the gas (region=0), the more general case would have be == Region=Geometry(Position)
    steps,collisions_counter = 0,0
    scatters = []  
    Weights = []
    MINI_DATA = []
    Tau_counter = 0
    while Region < len(Radius)-1:
        collision_pre = collisions_counter
        Position , Tau , Region , Energy ,steps,collisions_counter,Theta_momentum,Phi_momentum,Velocity_direction,scatter_atom,Weight,Mini_Data = Transport(Position,Energy,Theta_momentum,Phi_momentum,Region,collisions_counter,steps,reaction,BirthEnergy,Tau_counter)                
        MINI_DATA.append(Mini_Data)        
        if type(Region) == int:                
            Region = Geometry(Position)   
            Weights.append(Weight)
        else:
            break
        scatters.append(scatter_atom)                
        if collisions_counter > collision_pre:
             collisions += 1
        Tau_counter += Tau
        Stepers += steps 
    if Region != 'Error':
        return MINI_DATA        
    else:
        return ['Error','Error','Error','Error','Error','Error','Error','Error','Error','Error','Error','Error']
Parallel_start = tic()
if __name__ == '__main__':
    with closing(Pool()) as p:
        OUTPUT = p.map(Parallel,NeutronAttributes)
print('Time to spent in parallel:',tic()-Parallel_start)
stop = tic()


Timer = []
Final_Tally = []
NumberScatters = []
LastScatterAtom = []
for neutrons in OUTPUT:
    for i in reversed(range(len(neutrons))):
        if len(neutrons[i]) != 0:
            if type(neutrons[i][0]) != str:
                Final_Tally.append(neutrons[i])
                if neutrons[i][0][3] == 0 and neutrons[i][0][2] != 0:
                    for j in reversed(range(len(neutrons)-2)):
                        if neutrons[j][0][3] != 0:
                            LastScatterAtom.append(neutrons[j][0][3])
                            NumberScatters.append(neutrons[j][0][2])
                            break
                        else:
                            LastScatterAtom.append(neutrons[j][0][3])
                            NumberScatters.append(neutrons[j][0][2])
                            break
                else:
                    LastScatterAtom.append(neutrons[i][0][3])
                    NumberScatters.append(neutrons[i][0][2])
            else:
                break
            
# Put all the scatters in the data file of the neutron that is transported the whole way (AllData.txt)
for index in range(len(Final_Tally)):
    quant = Final_Tally[index][0]
    if quant[2] == 0:
        Final_Tally[index][0][2] = NumberScatters[index]
        Final_Tally[index][0][3] = LastScatterAtom[index]
        
ditched = 0 
#################################
####### Write out results #######
#################################
f = open('FullHistory.txt','w')
# NumberCollisions Energy(MeV) Velocity(m/s), Time(ns), VelocityUnitVector
f.write('Reaction meaning, 0=DT,1=DD,2=TT\n')
f.write('Radius,time of flight,')
f.write('Reaction, BirthEnergy (MeV), NumberOfCollisions, AtomicMassOfLastScatter, Energy(MeV), Time(ns), Weight\n')
for neutrons in OUTPUT:    
    for history_transports in neutrons:
        for history_collisions in history_transports:
            for j in range(len(history_collisions)):
                if type(history_collisions[j]) == str:
                    ditched += 1
                    break
                if j == 0 or j == 2 or j ==3:
                    f.write('%d ' % history_collisions[j])
                elif j == len(history_collisions)-1:          # If this is the last element in the row vector
                    f.write('%f\n' % history_collisions[j]) # Make a new line after this element 
                else:
                    f.write('%f ' % history_collisions[j])
f.close()
Timer = [col[0][5] for col in Final_Tally]
f = open('AllData.txt','w')
f.write('Reaction meaning, 0=DT,1=DD,2=TT\n')
f.write('Reaction, BirthEnergy (MeV), NumberOfCollisions, AtomicMassOfLastScatter, Energy(MeV), Time(ns), Weight\n')
for row in Final_Tally:
    if type(row[0][0])  == str:
        continue
    for j in range(len(row[0])):
        if j == len(row[0])-1:
            if type(row[0][j]) == float:
                f.write('%f\n' % row[0][j]) # Make a new line after this element 
            else:
                f.write('%d\n' * row[0][j])
        else:
            if type(row[0][j]) == float:
                f.write('%f ' % row[0][j])
            elif type(row[0][j]) == int:
                    f.write('%d ' % row[0][j])
f.close()
print('Time of simulation in seconds is',stop-start,'or',(stop-start)/60,'minutes with',len(OUTPUT),'Neutrons')
print('The number of neutrons that were below 50 KeV and were not recorded was:',ditched)
print('The average time to hit detectors at 20 (m) is:',sum(Timer)/len(Timer),'nanoseconds')
print('DTn fusions:',Fusions[0]*100/sum(Fusions),'%   DDn fusions:',Fusions[1]*100/sum(Fusions),'% TT2n fusions:',Fusions[2]*100/sum(Fusions))
print(' ')
