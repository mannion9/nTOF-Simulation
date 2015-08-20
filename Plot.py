'''Reaction meaning, 0=DT,1=DD'''
'''Reaction, BirthEnergy (MeV), NumberOfCollisions, AtomicMassOfLastScatter, Energy(MeV), Time(ns), Weight'''
Events = []
f = open('FullHistory.txt','r')
#f = open('AllData.txt','r')
f = f.read().splitlines()
f = f[2:]  # Gets rid of first line of text that explains columns 
for line in f:
    row = []
    line = line.split(' ')      # Seperates each number
    for i in line:
        row.append(float(i))
    Events.append(row)
    

# read in indexs 0
ai = 3 #scatter atom index
ti = 5 #time index
bi = 1 #birth energy index
ei = 4 #energy index
wi = 6 #weighting index
# Time
T_max = 1200
Time      = [col[ti] for col in Events if col[ti] < T_max]
CarbonT   = [col[ti] for col in Events if col[ti] < T_max and col[ai] == 12]
HydroT    = [col[ti] for col in Events if col[ti] < T_max and col[ai] == 1 ]
DT        = [col[ti] for col in Events if col[ti] < T_max and col[ai] == 2 ]
TT        = [col[ti] for col in Events if col[ti] < T_max and col[ai] == 3 ]
OT        = [col[ti] for col in Events if col[ti] < T_max and col[ai] == 16]
WT        = [col[ti] for col in Events if col[ti] < T_max and col[ai] == 183]
AT        = [col[ti] for col in Events if col[ti] < T_max and col[ai] == 40]
NT        = [col[ti] for col in Events if col[ti] < T_max and col[ai] == 14]
ALT       = [col[ti] for col in Events if col[ti] < T_max and col[ai] == 27]
# Energy
BirthEnergy = [col[bi] for col in Events]
Energy      = [col[ei] for col in Events]
CarbonE     = [col[ei] for col in Events if col[ai] == 12]
HydroE      = [col[ei] for col in Events if col[ai] == 1]
DE          = [col[ei] for col in Events if col[ai] == 2]
TE          = [col[ei] for col in Events if col[ai] == 3]
OE          = [col[ei] for col in Events if col[ai] == 16]
WE          = [col[ei] for col in Events if col[ai] == 183]
AE          = [col[ei] for col in Events if col[ai] == 40]
NE          = [col[ei] for col in Events if col[ai] == 14]
ALE         = [col[ei] for col in Events if col[ai] == 27]
# Weighting
FullWE   = [col[wi] for col in Events]
CarbonWE = [col[wi] for col in Events if col[ai] == 12]
HydroWE  = [col[wi] for col in Events if col[ai] == 1]
DWE      = [col[wi] for col in Events if col[ai] == 2]
TWE      = [col[wi] for col in Events if col[ai] == 3]
OWE      = [col[wi] for col in Events if col[ai] == 16]
WWE      = [col[wi] for col in Events if col[ai] == 183]
AWE      = [col[wi] for col in Events if col[ai] == 40]
NWE      = [col[wi] for col in Events if col[ai] == 14]
ALWE     = [col[wi] for col in Events if col[ai] == 27] 
 
FullWT   = [col[wi] for col in Events if col[ti] < T_max]
CarbonWT = [col[wi] for col in Events if col[ti] < T_max and col[ai] == 12]
HydroWT  = [col[wi] for col in Events if col[ti] < T_max and col[ai] == 1]
DWT      = [col[wi] for col in Events if col[ti] < T_max and col[ai] == 2]
TWT      = [col[wi] for col in Events if col[ti] < T_max and col[ai] == 3]
OWT      = [col[wi] for col in Events if col[ti] < T_max and col[ai] == 16]
WWT      = [col[wi] for col in Events if col[ti] < T_max and col[ai] == 183]
AWT      = [col[wi] for col in Events if col[ti] < T_max and col[ai] == 40]
NWT      = [col[wi] for col in Events if col[ti] < T_max and col[ai] == 14]
ALWT     = [col[wi] for col in Events if col[ti] < T_max and col[ai] == 27]    

def emptyCheck(x,weighting):
    if len(x) <= 1:
        return [[0],[0,0]]
    else:
        return plt.hist(x,bins=binner,weights=weighting)
        
import matplotlib.pyplot as plt
binner=500
# Energy
plt.figure(10)
BirthEnergy = plt.hist(BirthEnergy,bins=binner,weights=FullWE)
Energy = plt.hist(Energy,bins=binner,weights=FullWE)
CarbonE = emptyCheck(CarbonE,CarbonWE)
HydroE = emptyCheck(HydroE,HydroWE)
DE = emptyCheck(DE,DWE)
TE = emptyCheck(TE,TWE)
OE = emptyCheck(OE,OWE)
NE = emptyCheck(NE,NWE)
AE = emptyCheck(AE,AWE)
WE = emptyCheck(WE,WWE)
ALE = emptyCheck(ALE,ALWE)

# Time 
Time = plt.hist(Time,bins=binner,weights=FullWT)
CarbonT = emptyCheck(CarbonT,CarbonWT)
HydroT = emptyCheck(HydroT,HydroWT)
DT = emptyCheck(DT,DWT)
TT = emptyCheck(TT,TWT)
OT = emptyCheck(OT,OWT)
NT = emptyCheck(NT,NWT)
AT = emptyCheck(AT,AWT)
WT = emptyCheck(WT,WWT)
ALT = emptyCheck(ALT,ALWT)
plt.close(10)

# Energy
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
Birth_Spectrum, = plt.plot(BirthEnergy[1][:-1],BirthEnergy[0],color='blue',label='Birth Spectrum',linewidth=3)
Full_Spectrum, = plt.plot(Energy[1][:-1],Energy[0],color='black',label='Full Spectrum',linewidth=4)

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
CarbonT_Spectrum, = plt.plot(CarbonT[1][:-1],CarbonT[0],color='green',label='Carbon Scatters',linestyle='-.',linewidth=3)
FullTime_Spectrum, = plt.plot(Time[1][:-1],Time[0],color='black',label='Full Spectrum',linewidth=1)
plt.xlabel('Time (ns)'),plt.ylabel('Count'),plt.title('Neutron Time Of Flight' ),plt.yscale('log')
plt.legend(handles=[Full_Spectrum,CarbonT_Spectrum,HydroT_Spectrum,DT_Spectrum,TT_Spectrum,OT_Spectrum,NT_Spectrum,AT_Spectrum,ALT_Spectrum,WT_Spectrum],loc='upper right')
plt.tight_layout()
plt.xlim([min(Time[1][:-1]),max(Time[1][:-1])])
plt.show(1)
plt.show(2)

