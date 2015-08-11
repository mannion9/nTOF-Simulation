'''Reaction meaning, 0=DT,1=DD'''
'''Reaction, BirthEnergy (MeV), NumberOfCollisions, AtomicMassOfLastScatter, Energy(MeV), Velocity(m/s), Time(ns), X VelocityUnitVectorComponent, Y VelocityUnitVectorComponent, Z VelocityUnitVectorComponent'''
Events = []
f = open('AllData.txt','r')
f = f.read().splitlines()
f = f[2:]  # Gets rid of first line of text that explains columns 
for line in f:
    row = []
    line = line.split(' ')      # Seperates each number
    for i in line:
        row.append(float(i))
    Events.append(row)
#Events_1 = []
#f = open('/Users/mannion2/Desktop/Codes/Development/8-7/AllData.txt','r')
#f = f.read().splitlines()
#f = f[2:]  # Gets rid of first line of text that explains columns 
#for line in f:
#    row = []
#    line = line.split(' ')      # Seperates each number
#    for i in line:
#        row.append(float(i))
#    Events_1.append(row)
# Time
T_max = 1200
BirthTime = [col[6] for col in Events if col[6] < T_max]
Time      = [col[6] for col in Events if col[6] < T_max]
CarbonT   = [col[6] for col in Events if col[6] < T_max and col[3] == 12]
HydroT    = [col[6] for col in Events if col[6] < T_max and col[3] == 1 ]
DT        = [col[6] for col in Events if col[6] < T_max and col[3] == 2 ]
TT        = [col[6] for col in Events if col[6] < T_max and col[3] == 3 ]
OT        = [col[6] for col in Events if col[6] < T_max and col[3] == 16]
WT        = [col[6] for col in Events if col[6] < T_max and col[3] == 183]
AT        = [col[6] for col in Events if col[6] < T_max and col[4] == 40]
NT        = [col[6] for col in Events if col[6] < T_max and col[3] == 14]
ALT       = [col[6] for col in Events if col[6] < T_max and col[3] == 27]
# Energy
BirthEnergy = [col[1] for col in Events]
Energy      = [col[4] for col in Events]
#####Energy_1    = [col[4] for col in Events_1]
CarbonE     = [col[4] for col in Events if col[3] == 12]
HydroE      = [col[4] for col in Events if col[3] == 1]
DE          = [col[4] for col in Events if col[3] == 2]
TE          = [col[4] for col in Events if col[3] == 3]
OE          = [col[4] for col in Events if col[3] == 16]
WE          = [col[4] for col in Events if col[3] == 183]
AE          = [col[4] for col in Events if col[3] == 40]
NE          = [col[4] for col in Events if col[3] == 14]
ALE         = [col[4] for col in Events if col[3] == 27]


import matplotlib.pyplot as plt
binner=500
# Energy
plt.figure(10)
BirthEnergy = plt.hist(BirthEnergy,bins=binner)
Energy = plt.hist(Energy,bins=binner)
#####Energy_1 = plt.hist(Energy_1,bins=binner)
CarbonE = plt.hist(CarbonE,bins=binner)
HydroE = plt.hist(HydroE,bins=binner)
DE = plt.hist(DE,bins=binner)
TE = plt.hist(TE,bins=binner)
OE = plt.hist(OE,bins=binner)
NE = plt.hist(NE,bins=binner)
if len(AE) <= 1:
    AE = [[0],[0,0]] 
else:
    AE = plt.hist(AE,bins=binner)
if len(WE) <= 1:
    WE = [[0],[0,0]] 
else:
    WE = plt.hist(WE,bins=binner)
if len(ALE) <= 1:
    ALE = [[0],[0,0]]
else:
    ALE = plt.hist(ALE,bins=binner)

# Time 
Time = plt.hist(Time,bins=binner)
CarbonT = plt.hist(CarbonT,bins=binner)
HydroT = plt.hist(HydroT,bins=binner)
DT = plt.hist(DT,bins=binner)
TT = plt.hist(TT,bins=binner)
OT = plt.hist(OT,bins=binner)
NT = plt.hist(NT,bins=binner)
if len(AT) <= 1:
    AT = [[0],[0,0]] 
else:
    AT = plt.hist(AT,bins=binner)
if len(WT) <= 1:
    WT = [[0],[0,0]] 
else:
    WT = plt.hist(WT,bins=binner)
if len(ALT) <= 1:
    ALT = [[0],[0,0]]
else:
    ALT = plt.hist(ALT,bins=binner)

plt.close(10)

plt.figure(1)
Hydro_Spectrum, = plt.plot(HydroE[1][:-1],HydroE[0],color='black',label='Hydrogen Scatters')
DE_Spectrum, = plt.plot(DE[1][:-1],DE[0],color='red',label='Deutirium Scatters',linestyle='--')
TE_Spectrum, = plt.plot(TE[1][:-1],TE[0],color='black',label='Tritium Scatters',linestyle='-.')
OE_Spectrum, = plt.plot(OE[1][:-1],OE[0],color='yellow',label='Oxygen Scatters')
NE_Spectrum, = plt.plot(NE[1][:-1],NE[0],color='green',label='Nitrogen Scatters',linestyle='-.',linewidth=3)
AE_Spectrum, = plt.plot(AE[1][:-1],AE[0],color='red',label='Argon Scatters',linestyle='--',linewidth=5)
WE_Spectrum, = plt.plot(WE[1][:-1],WE[0],color='blue',label='Tungsten Scatters',linestyle='--')
ALE_Spectrum, = plt.plot(ALE[1][:-1],ALE[0],color='red',label='Aluminium Scatters',linewidth=1)
Carbon_Spectrum, = plt.plot(CarbonE[1][:-1],CarbonE[0],color='blue',label='Carbon Scatters',linestyle='-.',linewidth=3)
Full_Spectrum, = plt.plot(Energy[1][:-1],Energy[0],color='black',label='Full Spectrum',linewidth=1)
Birth_Spectrum, = plt.plot(BirthEnergy[1][:-1],BirthEnergy[0],color='blue',label='Birth Spectrum',linewidth=4)
######dude, = plt.plot(Energy_1[1][:-1],Energy_1[0],linewidth=1,color='red')
plt.xlabel('Energy (MeV)'),plt.ylabel('Count'),plt.title('Neutron Energy Spectrum' ),plt.yscale('log')
plt.legend(handles=[Birth_Spectrum,Full_Spectrum,DE_Spectrum,TE_Spectrum,OE_Spectrum,Hydro_Spectrum,Carbon_Spectrum,NE_Spectrum,AE_Spectrum,WE_Spectrum,ALE_Spectrum],loc='upper left')
plt.tight_layout()
plt.xlim([min(Energy[1][:-1]),max(Energy[1][:-1])])
# Time
plt.figure(2)
HydroT_Spectrum, = plt.plot(HydroT[1][:-1],HydroT[0],color='black',label='Hydrogen Scatters')
DT_Spectrum, = plt.plot(DT[1][:-1],DT[0],color='red',label='Deutirium Scatters',linestyle='--')
TT_Spectrum, = plt.plot(TT[1][:-1],TT[0],color='black',label='Tritium Scatters',linestyle='-.')
OT_Spectrum, = plt.plot(OT[1][:-1],OT[0],color='yellow',label='Oxygen Scatters')
NT_Spectrum, = plt.plot(NT[1][:-1],NT[0],color='green',label='Nitrogen Scatters',linestyle='-.',linewidth=3)
AT_Spectrum, = plt.plot(AT[1][:-1],AT[0],color='red',label='Argon Scatters',linestyle='--',linewidth=5)
WT_Spectrum, = plt.plot(WT[1][:-1],WT[0],color='blue',label='Tungsten Scatters',linestyle='--')
ALT_Spectrum, = plt.plot(ALT[1][:-1],ALT[0],color='red',label='Aluminium Scatters',linewidth=1)
CarbonT_Spectrum, = plt.plot(CarbonT[1][:-1],CarbonT[0],color='green',label='Nitrogen Scatters',linestyle='-.',linewidth=3)
FullTime_Spectrum, = plt.plot(Time[1][:-1],Time[0],color='black',label='Full Spectrum',linewidth=1)
plt.xlabel('Time (ns)'),plt.ylabel('Count'),plt.title('Neutron Time Of Flight' ),plt.yscale('log')
plt.legend(handles=[Full_Spectrum,CarbonT_Spectrum,HydroT_Spectrum,DT_Spectrum,TT_Spectrum,OT_Spectrum,NT_Spectrum,AT_Spectrum,ALT_Spectrum,WT_Spectrum],loc='upper right')
plt.tight_layout()
plt.xlim([min(Time[1][:-1]),max(Time[1][:-1])])
plt.show(1)
plt.show(2)

