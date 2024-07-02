from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

allsignal = []
bc = []
bc2 = []
names = []
start = []
start1 = []
start2 = []
m_start = []
stop = []
m_sum = []
k = 0

boron = Path("./hists/all_boron/boron")
background = Path("./hists/all_boron/background")
background2 = Path("./hists/all_boron/background2")
muon = Path("./hists/all_muons/result")

for file in boron.rglob("*.txt"):
    k+=1
    data = np.loadtxt(file)
    sum = np.sum(data.flatten())
    allsignal.append(sum)
    name = str(str((str(file).split('_', 2))[2]).split('.')[0])
    names.append(name.split('_')[0]+'\n-\n'+name.split('_')[1])
    start.append(name.split('_')[0])
    stop.append(name.split('_')[1])


for file in background.rglob("*.txt"):
    k+=1
    data = np.loadtxt(file)
    name = str(str((str(file).split('_', 2))[2]).split('.')[0])
    start1.append(name.split('_')[0])
    sum = np.sum(data.flatten())
    bc.append(sum)  
    

for file in background2.rglob("*.txt"):
    k+=1
    data = np.loadtxt(file)
    name = str(str((str(file).split('_', 2))[2]).split('.')[0])
    start2.append(name.split('_')[0])
    sum = np.sum(data.flatten())
    bc2.append(sum)

for file in muon.rglob("*.txt"):
    k+=1
    data = np.loadtxt(file)
    sum = np.sum(data.flatten())
    m_sum.append(sum)
    name = str(str((str(file).split('/'))[-1]).split('.')[0])
    m_start.append(name.split('_')[1])

#m_start = np.array(m_start).astype('int')

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

INTENSITY = {}
INTENSITY["Boron"] = allsignal
INTENSITY["Time"] = names
INTENSITY["Start"] = np.array(start).astype('int')
INTENSITY["Stop"] = np.array(stop).astype('int')

INTENSITY = pd.DataFrame(INTENSITY)
INTENSITY = INTENSITY.sort_values(by=['Start'])


I1 = {}
I1["Background"] = bc
I1["Start"] = np.array(start1).astype('int')
I1 = pd.DataFrame(I1)

INTENSITY = INTENSITY.merge(I1, left_on='Start', right_on='Start', how='inner', copy=False)


I2 = {}
I2["Background2"] = bc2
I2["Start"] = np.array(start2).astype('int')
I2 = pd.DataFrame(I2) 

INTENSITY = INTENSITY.merge(I2, left_on='Start', right_on='Start', how='inner', copy=False)


I3 = {}
I3["Muons"] = m_sum
I3["Start"] = np.array(m_start).astype('int')
I3 = pd.DataFrame(I3) 

INTENSITY = INTENSITY.merge(I3, left_on='Start', right_on='Start', how='inner', copy=False)

#print(I1, I2, I3)
INTENSITY["SumBackground"] = (INTENSITY["Background"]+INTENSITY["Background2"])/2
INTENSITY["Signal"] = INTENSITY["Boron"]-INTENSITY["SumBackground"]
INTENSITY["RelativeSelection"] = INTENSITY["Signal"]/INTENSITY["Muons"]



#print(INTENSITY[INTENSITY["RelativeSelection"] > 0.05] 
INTENSITY = INTENSITY[INTENSITY["RelativeSelection"] < 0.05]                  

INTENSITY["ErrSignal"] = np.sqrt(INTENSITY["Signal"]*(1-INTENSITY["RelativeSelection"]))    #poka ne uchla pogreshnost' fona
INTENSITY["ErrMuons"] = np.sqrt(INTENSITY["Muons"])
INTENSITY["ErrRelativeSelection"] = INTENSITY["RelativeSelection"] * np.sqrt((INTENSITY["ErrSignal"]/INTENSITY["Signal"])**2 + (INTENSITY["ErrMuons"]/INTENSITY["Muons"])**2)

INTENSITY["Time"] = (INTENSITY["Start"] + INTENSITY["Stop"])/2

data = np.loadtxt('./stat_pos_2023.txt', dtype='U')
#print(data)

FILE = {}
FILE["Start"] = data[1:, 1].astype('int')
FILE["Stop"] = data[1:, 2].astype('int')
FILE["Eff"] = data[1:, -1].astype('float')
FILE = pd.DataFrame(FILE)
FILE = FILE[FILE["Start"] > 19850]
FILE = FILE[FILE["Stop"] < 110200]

#perenormirovka po otrezkam dliny tau
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

tau = 2000

NORMINTENSITY = {}
NORMINTENSITY["Time"] = (np.array(INTENSITY["Stop"]) + np.array(INTENSITY["Start"]))/2
NORMINTENSITY["Stop"] = INTENSITY["Stop"]
NORMINTENSITY["Start"] = INTENSITY["Start"]
NORMINTENSITY["Signal"] = np.array(INTENSITY["Signal"])/(np.array(INTENSITY["Stop"]) - np.array(INTENSITY["Start"]))
NORMINTENSITY["Muons"] = np.array(INTENSITY["Muons"])/(np.array(INTENSITY["Stop"]) - np.array(INTENSITY["Start"]))
NORMINTENSITY["RelationSelection"] = np.array(INTENSITY["RelativeSelection"])

NORMINTENSITY["ErrMuons"] = NORMINTENSITY["Muons"]*INTENSITY["ErrMuons"]/INTENSITY["Muons"]
NORMINTENSITY["ErrSignal"] = NORMINTENSITY["Signal"]*INTENSITY["ErrSignal"]/INTENSITY["Signal"]
NORMINTENSITY["ErrRelativeSelection"] = np.array(INTENSITY["ErrRelativeSelection"])
NORMINTENSITY = pd.DataFrame(NORMINTENSITY)

INTENSITY1 ={}
x=[]
y=[]
z=[]
erry = []
errz = []

for t1 in range(20000, 110000, tau):
    t2 = t1+tau
    SET1 = NORMINTENSITY[NORMINTENSITY["Start"] < t2]
    SET = SET1[SET1["Stop"]> t1]
    SET[SET["Start"] < t1]["Start"] = t1
    SET[SET["Stop"] > t2]["Stop"] = t2
    l = np.sum(np.array(SET["Stop"]) - np.array(SET["Start"]))
    rho = np.sum((np.array(SET["Stop"]) - np.array(SET["Start"]))/l * np.array(SET["Signal"]))
    rho1 = np.sum((np.array(SET["Stop"]) - np.array(SET["Start"]))/l * np.array(SET["Muons"]))
    err = np.sqrt(np.sum((np.array(SET["Stop"]) - np.array(SET["Start"]))/l * (np.array(SET["ErrSignal"])/np.array(SET["Signal"]))**2)) * rho
    err1 = np.sqrt(np.sum((np.array(SET["Stop"]) - np.array(SET["Start"]))/l * (np.array(SET["ErrMuons"])/np.array(SET["Muons"]))**2)) * rho1

    y.append(rho)
    z.append(rho1)
    x.append(t1+tau/2)
    erry.append(err)
    errz.append(err1)

'''
plt.scatter(x, y)
plt.show()
'''

INTENSITY1 = {}

INTENSITY1["Signal"] = np.array(y)*tau
INTENSITY1["Muons"] = np.array(z)*tau
INTENSITY1["Time"] = x
INTENSITY1["ErrSignal"] = np.array(erry)*tau
INTENSITY1["ErrMuons"] = np.array(errz)*tau
INTENSITY1 = pd.DataFrame(INTENSITY1)
INTENSITY1["RelativeSelection"] = INTENSITY1["Signal"]/INTENSITY1["Muons"]
INTENSITY1["ErrRelativeSelection"] = np.sqrt((INTENSITY1["ErrSignal"]/INTENSITY1["Signal"])**2 + (INTENSITY1["ErrMuons"]/INTENSITY1["Muons"])**2)*INTENSITY1["RelativeSelection"]

average_e = INTENSITY1["RelativeSelection"].mean()                               # average efficiency value from experiment

#peremormirovka fila
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

FILE["Time"] = (FILE["Start"] + FILE["Stop"])/2

x = []
y = []

for t1 in range(20000, 110000, tau):
    t2 = t1+tau
    SETT1 = FILE[FILE["Start"] < t2]
    SETT = SETT1[SETT1["Stop"]> t1]
    SETT[SETT["Start"] < t1]["Start"] = t1
    SETT[SETT["Stop"] > t2]["Stop"] = t2
    l = np.sum(np.array(SETT["Stop"]) - np.array(SETT["Start"]))
    eff = np.sum((np.array(SETT["Stop"]) - np.array(SETT["Start"]))/l * np.array(SETT["Eff"]))
    
    y.append(eff)
    x.append(t1+tau/2)

FILE1 = {}
FILE1["Time"] = x
FILE1["Eff"] = y
FILE1 = pd.DataFrame(FILE1)

# plt.scatter(x, y)
# plt.scatter(FILE["Time"], FILE["Eff"])
# plt.show()

average_f = FILE["Eff"].mean()                                                  # average efficiency value from file

# plots:
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''
plt.plot(np.array(INTENSITY[INTENSITY["Start"]%1000 == 0]["Start"]), np.ones(len(INTENSITY[INTENSITY["Start"]%1000 == 0])), '-o')
plt.show()
'''
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(20,12))
# graph of the number of selected boron versus conventional time
ax.scatter(INTENSITY1["Time"], INTENSITY1["Signal"])
ax.errorbar(INTENSITY1["Time"], INTENSITY1["Signal"], yerr=INTENSITY1["ErrSignal"],fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
plt.xlabel('time', fontsize=20)
plt.ylabel('number of ivents', fontsize=20)
plt.grid()
plt.title("Borum decays", fontsize=20)
#plt.legend(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
#plt.show()

fig.savefig('re_select_boron_norm.png')

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(20,12))
# graph of the number of stopped muons versus conventional time
ax.scatter(INTENSITY1["Time"], INTENSITY1["Muons"], alpha=0.75, s=10)
ax.errorbar(INTENSITY1["Time"], INTENSITY1["Muons"], yerr=INTENSITY1["ErrMuons"],fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
plt.xlabel('time', fontsize=20)
plt.ylabel('number of ivents', fontsize=20)
plt.grid()
plt.title("Stopped muons", fontsize=20)
#plt.legend(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
#plt.show() 

fig.savefig('re_stopped_norm.png')

del plt
import matplotlib.pyplot as plt

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''
fig, ax = plt.subplots(figsize=(20,12))
# graph of the ratio of selected boron versus conventional time (efficiency)
ax.scatter(INTENSITY1["Time"], INTENSITY1["RelativeSelection"], alpha=0.75, s=10)
ax.errorbar(INTENSITY1["Time"], INTENSITY1["RelativeSelection"], yerr=INTENSITY1["ErrRelativeSelection"],fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
plt.xlabel('time', fontsize=20)
plt.ylabel('effectiveness', fontsize=20)
plt.grid()
plt.title("Effectiveness", fontsize=20)
#plt.legend(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
#plt.show()

fig.savefig('re_relative_hist_boron_exp.png')

del plt
import matplotlib.pyplot as plt
'''
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(20,12))
# graph of the ratio of selected boron to stopped muons and efficiency (from file) versus conventional time
ax.scatter(INTENSITY1["Time"], INTENSITY1["RelativeSelection"]/average_e, alpha=0.75, s=10)
ax.errorbar(INTENSITY1["Time"], INTENSITY1["RelativeSelection"]/average_e, yerr=INTENSITY1["ErrRelativeSelection"]/average_e,fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
plt.scatter(FILE1["Time"], FILE1["Eff"]/average_f, alpha=0.75, s=10)
plt.xlabel('time', fontsize=20)
plt.ylabel('effectiveness', fontsize=20)
plt.grid()
plt.title("Effectiveness", fontsize=20)
#plt.legend(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
#plt.show()

fig.savefig('re_relative_hist_norm_boron.png')

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(20,12))
# graph of the ratio of selected boron to stopped muons and efficiency (from file) versus conventional time
ax.scatter(INTENSITY1["Time"], INTENSITY1["RelativeSelection"]/average_e, alpha=0.75, s=10)
ax.errorbar(INTENSITY1["Time"], INTENSITY1["RelativeSelection"]/average_e, yerr=INTENSITY1["ErrRelativeSelection"]/average_e,fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
plt.scatter(FILE1["Time"], FILE1["Eff"]/average_f, alpha=0.75, s=10)
plt.xlabel('time', fontsize=20)
plt.ylabel('effectiveness', fontsize=20)
plt.xlim([19000, 110000])
plt.ylim([0.5, 2.00])
plt.grid()
plt.title("Effectiveness", fontsize=20)
#plt.legend(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
#plt.show()

fig.savefig('re_relative_hist_norm_boron_cut.png')

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(20,12))
ax.scatter(INTENSITY1["Time"], INTENSITY1["RelativeSelection"]/average_e - FILE1["Eff"]/average_f, alpha=0.75, s=10)
ax.errorbar(INTENSITY1["Time"], INTENSITY1["RelativeSelection"]/average_e - FILE1["Eff"]/average_f, yerr=INTENSITY1["ErrRelativeSelection"]/average_e,fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
plt.xlabel('time', fontsize=20)
plt.ylabel('effectiveness', fontsize=20)
plt.grid()
plt.title("Effectiveness", fontsize=20)
#plt.legend(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
#plt.show()

fig.savefig('comparison_difference.png')

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(20,12))
ax.scatter(INTENSITY1["Time"], (INTENSITY1["RelativeSelection"]/average_e) / (FILE1["Eff"]/average_f), alpha=0.75, s=10)
ax.errorbar(INTENSITY1["Time"], (INTENSITY1["RelativeSelection"]/average_e) / (FILE1["Eff"]/average_f), yerr=INTENSITY1["ErrRelativeSelection"]/average_e/(FILE1["Eff"]/average_f),fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
plt.xlabel('time', fontsize=20)
plt.ylabel('effectiveness', fontsize=20)
plt.grid()
plt.title("Effectiveness", fontsize=20)
#plt.legend(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
#plt.show()

fig.savefig('comparison_relation.png')

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

'''
fig, ax = plt.subplots(figsize=(20,12))
# graph of the ratio of selected boron to stopped muons and efficiency (from file) versus conventional time (normalized)
ax.scatter(NEWINTENSITY["Start"], NEWINTENSITY["RelativeSelection"]/average_e, alpha=0.75, s=10)
ax.errorbar(NEWINTENSITY["Start"], NEWINTENSITY["RelativeSelection"]/average_e, yerr=NEWINTENSITY["ErrRelativeSelection"],fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
plt.scatter(FILE["Start"], FILE["Eff"]/average_f, alpha=0.75, s=10)
plt.xlabel('time', fontsize=20)
plt.ylabel('effectiveness', fontsize=20)
plt.grid()
plt.title("Effectiveness", fontsize=20)
#plt.legend(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
#plt.show()

fig.savefig('relative_hist_norm_boron.png')

del plt
import matplotlib.pyplot as plt
'''
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
fig, ax = plt.subplots(figsize=(20,12))
# graph of the ratio of selected boron to stopped muons and efficiency (from file) versus conventional time
ax.scatter(INTENSITY["Start"], INTENSITY["RelativeSelection"], alpha=0.75, s=10)
ax.errorbar(INTENSITY["Start"], INTENSITY["RelativeSelection"], yerr=INTENSITY["ErrRel"],fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
plt.scatter(FILE["Start"], FILE["Eff"]/average_f, alpha=0.75, s=10)
plt.xlabel('time', fontsize=20)
plt.ylabel('effectiveness', fontsize=20)
plt.xlim([19000, 110000])
#plt.ylim([0.9, 1.1])
plt.grid()
plt.title("Effectiveness", fontsize=20)
#plt.legend(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
#plt.show()

fig.savefig('relative_hist_boron.png')
"""
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''
colors = plt.cm.viridis(INTENSITY["Signal"]/np.max(INTENSITY["Signal"]))
colors = sns.color_palette(colors, as_cmap=True)

fig, axs = plt.subplots(figsize=(len(names), 40), constrained_layout=True)
axs = sns.barplot(x="Time", y="Signal", data=INTENSITY, palette=colors, dodge=False)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid(True)
fig.savefig('boron_signal.png')

colors = plt.cm.viridis(INTENSITY["Relative selection"]/np.max(INTENSITY["Relative selection"]))
colors = sns.color_palette(colors, as_cmap=True)

fig, axs = plt.subplots(figsize=(len(names), 40), constrained_layout=True)
axs = sns.barplot(x="Time", y="Relative selection", data=INTENSITY, palette=colors, dodge=False)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid(True)
fig.savefig('relation_boron_signal.png')
'''
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

'''
fig, axs = plt.subplots(sharey=True, figsize=(40, 20))
sns.boxplot(x=INTENSITY["Boron"], flierprops={"marker": "x"}, boxprops={"facecolor": (.4, .6, .8, .5)}, medianprops={"color": "coral"}, linewidth=5)
axs.grid(True)
plt.tick_params(axis='both', which='major', labelsize=20)
fig.savefig('distribution_boron.png')

fig, axs = plt.subplots(sharey=True, figsize=(40, 20))
sns.boxplot(x=INTENSITY["SumBackground"], flierprops={"marker": "x"}, boxprops={"facecolor": (.4, .6, .8, .5)}, medianprops={"color": "coral"}, linewidth=5)
axs.grid(True)
plt.tick_params(axis='both', which='major', labelsize=20)
fig.savefig('distribution_background.png')
'''