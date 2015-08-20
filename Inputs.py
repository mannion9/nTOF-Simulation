import sys
sys.path.append('BackEnd/')
import PypyNumpy as np
import math
from DataBaseReadInFunctions import elastic_read_in,scatter_read_in

###################################
############ Inputs ###############
###################################
Number = 1E6                                        # Number of neutrons to simulate
T_max = 5                                          # Maximum temperature of fuel in keV
BurnWidth = 150E-12                                 # FWHM of burn in seconds
TotalBurnSigma = 10                                 # Number of sigma to create nuetrons at (centered at time=0)
TimeStep = 1E-12                                    # Time step for creating neutrons
fusion_list = [1,1,1]                               # Set to zero to shut off reaction [DT,DD,TT]

Ballabio = 1                                        # Set to 0 for brysk and 1 for ballabio
implosion = 1                                       # Set to 1 for implosion or 0 for static
PeakVelocity = 150E3 
M_hs = 15E-6                                        # Total mass of hot spot in grams

M_ice = 171E-6-M_hs                                 # Mass of ice in grams
CH_rho_r = 200E-3                                   # Rho*R for ablator in g/cm^2

###################################
######## Aparatus Features ########
###################################
''' DT/DT/CH/VACCUM'''
Radius = [25E-6,50E-6,50E-6+70E-6,20,30]                                               # Radius of each region in meters
Density = np.vector([M_hs/(4*math.pi*(Radius[0]*1E2)**3/3),M_ice/(4*math.pi*((Radius[1]*1E2)**3-(Radius[0]*1E2)**3)/3),CH_rho_r/(Radius[2]*1E2-Radius[1]*1E2),1E-10,1E-10]) #g/cm^3 
A = [[2,3],[2,3],[12,1],[1],[1]]                                               # Atomic mass of each substance with that region
Ratio = np.vector([[.5,.5],[.5,.5],[(1/(1+1.3)),(1.3/(1+1.3))],[1],[1]])       # Ratio of each substance that make up material     
''' DT/DT/CH/VACCUM/AL/VACCUM/AIR/W '''
#Radius = [33E-6,53E-6,53E-6+70E-6,.5E-2,.5E-2+.5E-3,11,19.98,20,30]                                               # Radius of each region in meters
#Density = np.vector([M_hs/(4*math.pi*(Radius[0]*1E2)**3/3),M_ice/(4*math.pi*((Radius[1]*1E2)**3-(Radius[0]*1E2)**3)/3),CH_rho_r/(Radius[2]*1E2-Radius[1]*1E2),1E-10,2.7,1E-10,.0012,18.3,1E-10]) #g/cm^3      
#A = [[2,3],[2,3],[12,1],[1],[27],[1],[16,14,40],[183],[1]]                                               # Atomic mass of each substance with that region
#Ratio = np.vector([[.5,.5],[.5,.5],[(1/(1+1.3)),(1.3/(1+1.3))],[1],[1],[1],[.19,.8,.01],[1],[1]])       # Ratio of each substance that make up material
''' DT/DT/CH/VACCUM/AL/VACCUM/AIR '''
#Radius = [33E-6,53E-6,53E-6+70E-6,.5E-2,.5E-2+.5E-3,11,20,30]                                               # Radius of each region in meters
#Density = np.vector([M_hs/(4*math.pi*(Radius[0]*1E2)**3/3),M_ice/(4*math.pi*((Radius[1]*1E2)**3-(Radius[0]*1E2)**3)/3),CH_rho_r/(Radius[2]*1E2-Radius[1]*1E2),1E-10,2.7,1E-10,.0012,1E-10]) #g/cm^3      
#A = [[2,3],[2,3],[12,1],[1],[27],[1],[16,14,40],[1]]                                               # Atomic mass of each substance with that region
#Ratio = np.vector([[.5,.5],[.5,.5],[(1/(1+1.3)),(1.3/(1+1.3))],[1],[1],[1],[.19,.8,.01],[1]]) 

R0 = Radius[0]                                                                                  # Edge radius of neutron production (largest radius at which a neuutron may be born)
MaxRadius = Radius[2] # Where does this implosion stop?
NumberShells = 50                                                                               # number of shells to break capsul into for temperature gradient
detector_position = np.vector([20,0,0])
###################################
######## Database Inputs  #########
###################################
Energy_dic_atoms = [2,3,12,1,16,183,14,40,27]  # This is the databases index for each atom. When adding new elements ensure that the index matches that of this list (rember start at zero)
Energy_dic = [ [] for i in range(len(Energy_dic_atoms))] # Holds all the energy of data that we have for differnetial cross section 
Cosines = [ [] for i in range(len(Energy_dic_atoms))]
Probability = [ [] for i in range(len(Energy_dic_atoms))]
''' Elastic Scattering Cross Section '''
Cross_Energy,Cross = [ [] for i in range(len(Energy_dic_atoms))],[ [] for i in range(len(Energy_dic_atoms))]      
Cross_Energy[0] , Cross[0] = elastic_read_in('Database/Elastic/Hydrogen-2-ElasticXC.txt')
Cross_Energy[1] , Cross[1] = elastic_read_in('Database/Elastic/Hydrogen-3-ElasticXC.txt')
Cross_Energy[2] , Cross[2] = elastic_read_in('DataBase/Elastic/Carbon-12-ElasticXC.txt')
Cross_Energy[3] , Cross[3] = elastic_read_in('DataBase/Elastic/Hydrogen-1-ElasticXC.txt')
Cross_Energy[4] , Cross[4] = elastic_read_in('DataBase/Elastic/Oxygen-16-ElasticXC.txt')
Cross_Energy[5] , Cross[5] = elastic_read_in('Database/Elastic/Tun-183-ElasticXC.txt')
Cross_Energy[6] , Cross[6] = elastic_read_in('Database/Elastic/Nitrogen-14-ElasticXC.txt')
Cross_Energy[7] , Cross[7] = elastic_read_in('Database/Elastic/Argon-40-ElasticXC.txt')
Cross_Energy[8] , Cross[8] = elastic_read_in('Database/Elastic/Aluminium-27-ElasticXC.txt')
''' Differnetial Scattering Data in CM frame'''
Energy_dic[0] , Cosines[0] , Probability[0] = scatter_read_in('DataBase/Scattering/Hydrogen-2-ElasticScatterCrossSection.txt',1)
Energy_dic[1] , Cosines[1] , Probability[1] = scatter_read_in('DataBase/Scattering/Hydrogen-3-ElasticScatterCrossSection.txt',1)
Energy_dic[2] , Cosines[2] , Probability[2] = scatter_read_in('DataBase/Scattering/Carbon-12-ElasticScatterCrossSection.txt',1)
Energy_dic[3] , Cosines[3] , Probability[3] = scatter_read_in('DataBase/Scattering/Hydrogen-1-ElasticScatterCrossSection.txt',1)
Energy_dic[4] , Cosines[4] , Probability[4] = scatter_read_in('DataBase/Scattering/Oxygen-16-ElasticScatterCrossSection.txt',1)
Energy_dic[5] , Cosines[5] , Probability[5] = scatter_read_in('DataBase/Scattering/Tungsten-183-ElasticScatterCrossSection.txt',1)
Energy_dic[6] , Cosines[6] , Probability[6] = scatter_read_in('DataBase/Scattering/Nitrogen-14-ElasticScatterCrossSection.txt',1)
Energy_dic[7] , Cosines[7] , Probability[7] = scatter_read_in('DataBase/Scattering/Argon-40-ElasticScatterCrossSection.txt',1)
Energy_dic[8] , Cosines[8] , Probability[8] = scatter_read_in('DataBase/Scattering/Aluminium-27-ElasticScatterCrossSection.txt',1)
''' Birth Spectrum of T+T->2n+alpha reaction '''
f = open('Database/tt2n-spec.dat',"r")
f = f.read().splitlines()
data = [line.replace('\t',' ') for line in f]
data = [line.split('  ') for line in data]
data = [ [float(i) for i in data[j]] for j in range(len(data)) ]
tt2n_spectrum = [[col[0] for col in data],[col[1] for col in data]]
''' Reactivity of T+T->2n+alpha reaction '''
f = open('Database/tt2n-reac.txt','r')
f = f.read().splitlines()
tt2n_react = [row.split(',') for row in f]
tt2n_react = [[float(col[0]) for col in tt2n_react],[float(col[1]) for col in tt2n_react]]

###################################
#########  Constants ##############
###################################
BurnSigma = BurnWidth/(2*math.sqrt(2*math.log(2)))
N_A = 6.022E23                                              # Avagadros number
BarnToCmSquare = 1E-24                                      # Conversion factor from barns to m^2 
Amu = 931.494028                                            # 1 Atomic mass unit in MeV
Mn = 1.00866491597*Amu                                      # Neutron mass
MH2 = 2.013553212724*Amu                                    # Deuterium mass
MH3 = 3.0155007134*Amu                                      # Tritium mass
MHe3 = 3.0149322473*Amu                                     # Helium 3 mass
MHe4 = 4.001506179127*Amu                                   # Helium 4 mass

MassFraction = [[mass*fraction for mass,fraction in zip(A[i],Ratio[i])] for i in range(len(A))]
NDensity=[Density[i]*N_A/sum(MassFraction[i]) for i in range(len(Radius))]
N_i_frac = [[NDensity[i]*Ratio[i][j] for j in range(len(Ratio[i]))] for i in range(len(Ratio))]
MHe = [MHe4,MHe3]
Q_values = [(MH2+MH3)-(MHe4+Mn),((MH2+MH2)-(Mn+MHe3))]                                 # Q values
T_n_zero_temp = [(MHe4+Q_values[0]/2)*Q_values[0]/(MH2+MH3),(MH3+Q_values[1]/2)*Q_values[1]/(MH2+MH2)]
P_zero_temp = [math.sqrt(T_n_zero_temp[0]**2+2*Mn*T_n_zero_temp[0]),math.sqrt(T_n_zero_temp[1]**2+2*Mn*T_n_zero_temp[1])]
alpha = [[[5.30509, 2.4736e-3, 1.84, 1.3818],[5.1068e-4, 7.6223e-3, 1.78, 8.7691e-5]],[[4.69515,-0.040729,0.47,0.81844],[1.7013E-3,0.16888,0.49,7.9460E-4]]]

# parameters from bosh and hale paper
# [DTn,DDn]
m_r = [1124656,937814]                                               # reduced mass of reactants in keV
B_g = [34.3827,31.3970]
c = [[1.17302E-9 , 1.51361E-2 , 7.51886E-2 , 4.60643E-3 , 1.35E-2 , -1.0675E-4,1.366E-5] ,[5.65718E-12,3.41267E-3,1.99167E-3,0,1.05060E-5,0,0]]

      

