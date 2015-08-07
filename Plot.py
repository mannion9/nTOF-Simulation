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
# Energy
BirthEnergy = [col[1] for col in Events]
Energy      = [col[4] for col in Events]
CarbonE     = [col[4] for col in Events if col[3] == 12]
HydroE      = [col[4] for col in Events if col[3] == 1]
DE          = [col[4] for col in Events if col[3] == 2]
TE          = [col[4] for col in Events if col[3] == 3]
OE          = [col[4] for col in Events if col[3] == 16]
WE          = [col[4] for col in Events if col[3] == 183]
AE          = [col[4] for col in Events if col[3] == 40]
NE          = [col[4] for col in Events if col[3] == 14]

import matplotlib.pyplot as plt
binner=500
# Energy
plt.figure(10)
BirthEnergy = plt.hist(BirthEnergy,bins=binner)
Energy = plt.hist(Energy,bins=binner)
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

# Time 
Time = plt.hist(Time,bins=binner)
CarbonT = plt.hist(CarbonT,bins=binner)
HydroT = plt.hist(HydroT,bins=binner)
DT = plt.hist(DT,bins=binner)
TT = plt.hist(TT,bins=binner)
OT = plt.hist(OT,bins=binner)
NT = plt.hist(NT,bins=binner)
if len(AT) == 0:
    AT = [[0],[0,0]] 
else:
    AT = plt.hist(AE,bins=binner)
if len(WT) == 0:
    WT = [[0],[0,0]] 
else:
    WE = plt.hist(WE,bins=binner)

plt.close(10)

plt.figure(1)
Hydro_Spectrum, = plt.plot(HydroE[1][:-1],HydroE[0],color='black',label='Hydrogen Scatters')
OE_Spectrum, = plt.plot(OE[1][:-1],OE[0],color='yellow',label='Oxygen Scatters')
NE_Spectrum, = plt.plot(NE[1][:-1],NE[0],color='green',label='Nitrogen Scatters',linestyle='-.',linewidth=3)
WE_Spectrum, = plt.plot(WE[1][:-1],WE[0],color='blue',label='Tungsten Scatters',linestyle='--')
DE_Spectrum, = plt.plot(DE[1][:-1],DE[0],color='red',label='Deutirium Scatters',linestyle='--')#,marker='^')
TE_Spectrum, = plt.plot(TE[1][:-1],TE[0],color='black',label='Tritium Scatters',linestyle='-.')#,marker='^')
AE_Spectrum, = plt.plot(AE[1][:-1],AE[0],color='red',label='Argon Scatters',linestyle='--',linewidth=5)#,marker='^')
Birth_Spectrum, = plt.plot(BirthEnergy[1][:-1],BirthEnergy[0],color='blue',label='Birth Spectrum',linewidth=4)
Full_Spectrum, = plt.plot(Energy[1][:-1],Energy[0],color='black',label='Full Spectrum',linewidth=4)
Carbon_Spectrum, = plt.plot(CarbonE[1][:-1],CarbonE[0],color='green',label='Carbon Scatters')
plt.xlabel('Energy (MeV)'),plt.ylabel('Count'),plt.title('Neutron Energy Spectrum' ),plt.yscale('log')
plt.legend(handles=[Birth_Spectrum,Full_Spectrum,DE_Spectrum,TE_Spectrum,OE_Spectrum,Hydro_Spectrum,Carbon_Spectrum,NE_Spectrum,AE_Spectrum,WE_Spectrum],loc='upper left')
#plt.legend(handles=[Birth_Spectrum,Full_Spectrum,DE_Spectrum,TE_Spectrum,OE_Spectrum,Hydro_Spectrum,Carbon_Spectrum],loc='upper left')
#plt.legend(handles=[Birth_Spectrum,Full_Spectrum,DE_Spectrum,TE_Spectrum,Carbon_Spectrum],loc='upper left')
plt.tight_layout()
plt.xlim([min(Energy[1][:-1]),max(Energy[1][:-1])])
#plt.legend(handles=[Birth_Spectrum,Full_Spectrum,Carbon_Spectrum,Hydro_Spectrum,DE_Spectrum,TE_Spectrum,WE_Spectrum,OE_Spectrum],loc='upper left')

# Time
plt.figure(2)
HydroT_Spectrum, = plt.plot(HydroT[1][:-1],HydroT[0],color='black',label='Hydrogen Scatters')
DT_Spectrum, = plt.plot(DT[1][:-1],DT[0],color='red',label='Deutirium Scatters',linestyle='-.')#,marker='^')
TT_Spectrum, = plt.plot(TT[1][:-1],TT[0],color='blue',label='Tritium Scatters',linestyle='-.')#,marker='^')
OT_Spectrum, = plt.plot(OT[1][:-1],OT[0],color='yellow',label='Oxygen Scatters')#,marker='^')
NT_Spectrum, = plt.plot(NT[1][:-1],NT[0],color='green',label='Nitrogen Scatters')#,marker='^')
AT_Spectrum, = plt.plot(AT[1][:-1],AT[0],color='red',label='Argon Scatters',linewidth=4)#,marker='^')
FullTime_Spectrum, = plt.plot(Time[1][:-1],Time[0],color='black',label='Full Spectrum',linewidth=4)
CarbonT_Spectrum, = plt.plot(CarbonT[1][:-1],CarbonT[0],color='green',label='Nitrogen Scatters',linestyle='-.',linewidth=3)
plt.xlabel('Time (ns)'),plt.ylabel('Count'),plt.title('Neutron Time Of Flight' ),plt.yscale('log')
#plt.legend(handles=[Full_Spectrum,DT_Spectrum,TT_Spectrum],loc='upper right')
#plt.legend(handles=[Birth_Spectrum,Full_Spectrum,Carbon_Spectrum,Hydro_Spectrum,DE_Spectrum,TE_Spectrum,OT_Spectrum],loc='upper right')
plt.legend(handles=[Birth_Spectrum,Full_Spectrum,Carbon_Spectrum,Hydro_Spectrum,DE_Spectrum,TE_Spectrum,OT_Spectrum,NT_Spectrum,AT_Spectrum],loc='upper right')
plt.tight_layout()
plt.xlim([min(Time[1][:-1]),max(Time[1][:-1])])
plt.show(1)
plt.show(2)

#plt.figure(1)
#scattering = plt.hist(Times,bins=1000)
#plt.xlabel('Time (ns)')
#plt.ylabel('Count')
#plt.yscale('log')
#plt.title('nTOF to reach 20m for %d neutrons, T_max = %d keV with Hatarik' %(len(Times),10))
#plt.show()
#
#plt.figure(2)
#energies = plt.hist(Energy,bins=10000)
#plt.xlabel('Energy (MeV)')
#plt.ylabel('Count')
#plt.yscale('log')
#plt.title('Energies when passing 20m for %d neutrons, T_max = %d keV with Hatarik' %(len(Times),10))
#plt.show()
#
#plt.figure(4)
#last_scatter = plt.hist(BirthEnergy,bins=10000)
#plt.yscale('log')
#plt.xlabel('Energy (MeV)'),plt.ylabel('Count')
#plt.title('Birth Energy Spectrum with DT/DD/TT reactions')
#plt.show()
