from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from statistics import mean

def convert_date(Ns, data_dates):
    ans = []
    for N in Ns:
        for idx, item in enumerate(data_dates['f_stop']):
            if item > N:
                break
        ans.append(data_dates['date'][idx])
    return np.array(ans, dtype='datetime64')

sums = []
names = []
start = []
stop = []

p = Path("./hists/all_muons/result")

for file in p.rglob("*.txt"):
    data = np.loadtxt(file)
    sum = np.sum(data.flatten())
    sums.append(sum)
    name = str(str((str(file).split('/'))[-1]).split('.')[0])
    names.append(name.split('_')[0]+'\n-\n'+name.split('_')[1])
    start.append(name.split('_')[0])
    stop.append(name.split('_')[1])

start = np.array(start).astype('int')
stop = np.array(stop).astype('int')


re_sums = []
re_start = []
k = 0

r = Path("./hists/all_muons/reader_result")

for file in r.rglob("*.txt"):
    k+=1
    data = np.loadtxt(file)
    sum = np.sum(data.flatten())
    re_sums.append(sum)
    name = str(str((str(file).split('/'))[-1]).split('.')[0])
    re_start.append(name.split('_')[1])

re_start = np.array(re_start).astype('int')


INTENSITY = {}
INTENSITY["Sum"] = sums
# INTENSITY["Time"] = names
INTENSITY["Start"] = start
INTENSITY['Stop'] = stop
INTENSITY["Time"] = 1/2*(stop+start)


INTENSITY = pd.DataFrame(INTENSITY)
INTENSITY = INTENSITY.sort_values(by=['Start'])

I = {}
I["Select"] = re_sums
I["Start"] = re_start

I = pd.DataFrame(I)

INTENSITY = INTENSITY.merge(I, left_on='Start', right_on='Start', how='inner', copy=False)

INTENSITY = INTENSITY[INTENSITY['Sum'] != 0]
INTENSITY["Ratio"] = INTENSITY["Select"]/INTENSITY["Sum"]

# INTENSITY["Err_sel"] = np.sqrt(INTENSITY["Select"])
INTENSITY["Err_sum"] = np.sqrt(INTENSITY["Sum"])
# INTENSITY["Rel_err"] = np.sqrt((INTENSITY["Err_sel"]/INTENSITY["Sum"])**2 + (INTENSITY["Err_sum"]*INTENSITY["Select"]/(INTENSITY["Sum"]**2))**2)

INTENSITY["Err_select"] = np.sqrt(INTENSITY["Select"]*(1-(INTENSITY["Select"]/INTENSITY["Sum"])))
INTENSITY["Err_ratio"] = INTENSITY["Err_select"]/INTENSITY["Sum"]
	
data = np.loadtxt('./stat_pos_2023.txt', dtype='U',skiprows=1)
#print(data)

FILE = {}
FILE["Start"] = data[1:, 1].astype('int')
FILE["Stop"] = data[1:, 2].astype('int')
FILE["Time"] = 1/2*(data[1:, 1].astype('int') + data[1:, 2].astype('int'))
FILE["Eff"] = data[1:, -1].astype('float')
#print(INTENSITY["Sum"])
#print(INTENSITY["Err_sel"].mean())
FILE = pd.DataFrame(FILE)
average_f = mean(FILE["Eff"])                                       # average efficiency value from file
average_e = mean(INTENSITY["Ratio"])                   # average efficiency value from experiment

dates = []
f_stops = []
f_starts = []
for i in range(len(data)):
    line = data[i]
    f_starts.append(float(line[1]))
    f_stops.append(float(line[2]))
    x = str(line[0]).split("_", 1)[1]
    year = ('20'+ x.split(".")[2])
    month = x.split(".")[1]
    day = x.split(".")[0]
    dates.append(np.datetime64(year + '-' + month + '-' + day))
    
    
data_dates = pd.DataFrame({'f_start':f_starts, 'f_stop': f_stops, 'date': dates})

# print(len(convert_date(np.array(INTENSITY["Start"]), data_dates)), len(INTENSITY["Start"]))


# plots:
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# fig, ax = plt.subplots(figsize=(20,12))
# # graph of the number of selected decayed muons versus conventional time
# plt.scatter(convert_date(np.array(INTENSITY["Start"]), data_dates), INTENSITY["Select"])
# plt.errorbar(convert_date(np.array(INTENSITY["Start"]), data_dates), INTENSITY["Select"], yerr=INTENSITY["Err_select"],fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
# plt.xlabel('time', fontsize=20)
# plt.ylabel('number of events', fontsize=20)
# plt.grid()
# plt.title("Decayed muons", fontsize=30)
# #plt.legend(fontsize=20)
# plt.tick_params(axis='both', which='major', labelsize=20)
# #plt.show()

# fig.savefig('select_hist_data.png')

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# fig, ax = plt.subplots(figsize=(20,12))
# # graph of the number of stopped muons versus conventional time
# plt.scatter(convert_date(np.array(INTENSITY["Start"]), data_dates), INTENSITY["Sum"], alpha=0.75, s=10)
# plt.errorbar(convert_date(np.array(INTENSITY["Start"]), data_dates), INTENSITY["Sum"], yerr=INTENSITY["Err_sum"],fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
# plt.xlabel('time', fontsize=20)
# plt.ylabel('number of events', fontsize=20)
# plt.grid()
# plt.title("Stopped muons", fontsize=30)
# #plt.legend(fontsize=20)
# plt.tick_params(axis='both', which='major', labelsize=20)
# #plt.show()

# fig.savefig('stopped_hist_data.png')

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# fig, ax = plt.subplots(figsize=(20,12))
# # graph of the ratio of selected decayed muons versus conventional time (efficiency)
# plt.scatter(convert_date(np.array(INTENSITY["Start"]), data_dates), INTENSITY["Ratio"], alpha=0.75, s=10)
# plt.errorbar(convert_date(np.array(INTENSITY["Start"]), data_dates), INTENSITY["Ratio"], yerr=INTENSITY["Err_ratio"],fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
# plt.xlabel('time', fontsize=20)
# plt.ylabel('efficiency', fontsize=20)
# plt.grid()
# plt.title("Efficiency", fontsize=30)
# plt.legend(fontsize=20)
# plt.tick_params(axis='both', which='major', labelsize=20)
# #plt.show()

# fig.savefig('relative_hist_exp_data.png')

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# fig, ax = plt.subplots(figsize=(20,12))
# # graph of the ratio of detached decayed muons to stopped muons and efficiency (from file) versus conventional time (normalized)
# plt.scatter(convert_date(np.array(INTENSITY["Start"]), data_dates), INTENSITY["Ratio"]/average_e, alpha=0.75, s=10)
# plt.errorbar(convert_date(np.array(INTENSITY["Start"]), data_dates), INTENSITY["Ratio"]/average_e, yerr=INTENSITY["Err_ratio"],fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
# plt.scatter(FILE["Start"], FILE["Eff"]/average_f, alpha=0.75, s=10)
# plt.xlabel('time', fontsize=20)
# plt.ylabel('efficiency', fontsize=20)
# plt.grid()
# plt.title("Efficiency", fontsize=30)
# plt.legend(fontsize=20)
# plt.tick_params(axis='both', which='major', labelsize=20)
# #plt.show()

# fig.savefig('relative_hist_data.png')

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

FILE = FILE[FILE['Start']>20000]
FILE = FILE[FILE['Start']<110000]


INTENSITY = INTENSITY[INTENSITY["Ratio"] > 0.05]  
average_f = mean(FILE["Eff"])                                       
average_e = mean(INTENSITY["Ratio"])    

tims1 = np.array(INTENSITY['Time'])
tims2 = np.array(FILE['Time'])
indices1 = []
indices2 = []
for j, t in enumerate(tims1):
    find = False
    k = 0
    for i in range(len(tims2)):
        if (np.array(FILE['Stop'])[i] > t and np.array(FILE['Start'])[i] <= t) or (np.array(FILE['Stop'])[i] >= t and np.array(FILE['Start'])[i] < t):
            find = True
            k = i
    
    if not find:
        print(t, 'POLUNDRAAA')
    else:
        indices1.append(j)
        indices2.append(k)

indices1 = np.array(indices1)
indices2 = np.array(indices2)

fig, ax = plt.subplots(figsize=(20,13))
# graph of the ratio of detached decayed muons to stopped muons and efficiency (from file) versus conventional time (normalized)
plt.scatter(convert_date(np.array(INTENSITY["Time"])[indices1], data_dates), np.array(INTENSITY["Ratio"])[indices1]/average_e, alpha=0.75, s=10)
plt.errorbar(convert_date(np.array(INTENSITY["Time"])[indices1], data_dates), np.array(INTENSITY["Ratio"])[indices1]/average_e, yerr=np.array(INTENSITY["Err_ratio"])[indices1]/average_e,fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
plt.scatter(convert_date(np.array(FILE["Time"])[indices2], data_dates), np.array(FILE["Eff"])[indices2]/average_f, alpha=0.75, s=10)
plt.xlabel('time', fontsize=40)
plt.ylabel('efficiency', fontsize=40)
# plt.xlim([19000, 110000])
# plt.ylim([0.95, 1.08])
# plt.legend(fontsize=10)
plt.grid()
plt.title("Efficiency", fontsize=50)
plt.legend(['experimental data (muons)', 'Monte-Carlo'], fontsize=30, loc=1)
plt.tick_params(axis='both', which='major', labelsize=30)

plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
         rotation_mode="anchor")

fig.savefig('relative_hist_cut_data.png')

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# fig, ax = plt.subplots(figsize=(20,12))
# # graph of the ratio of detached decayed muons to stopped muons and efficiency (from file) versus conventional time (normalized)
# plt.scatter(convert_date(np.array(FILE["Start"]), data_dates), FILE["Eff"], alpha=0.75, s=10)
# plt.xlabel('time', fontsize=40)
# plt.ylabel('efficiency', fontsize=40)
# plt.grid()
# plt.title("Efficiency(time)", fontsize=50)
# plt.legend(fontsize=20)
# plt.tick_params(axis='both', which='major', labelsize=30)
# #plt.show()

# fig.savefig('Efficiency_data.png')

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# date1 = []
# dates = data[1:, 0].astype('str')
# #print(len(dates))

# for i in range(len(dates)):
#     x = dates[i].split("_", 1)[1]
#     year = ('20'+ x.split(".")[2])
#     month = x.split(".")[1]
#     day = x.split(".")[0]
#     date1.append(np.datetime64(year + '-' + month + '-' + day))
# #print(date1)

# FILE['Date'] = date1 
'''
fig, ax = plt.subplots(figsize=(20,12))
plt.scatter(FILE["Date"], FILE["Start"], alpha=0.75, s=10)
plt.xlabel('time', fontsize=40)
plt.ylabel('file number', fontsize=40)
plt.grid()
#plt.title("Effectiveness(time)", fontsize=50)
plt.legend(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=30)
plt.show()

fig.savefig('Effectiveness.png')
'''
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(20,13))
# graph of the ratio of detached decayed muons to stopped muons and efficiency (from file) versus conventional time (normalized)
plt.scatter(convert_date(np.array(INTENSITY["Time"])[indices1], data_dates), np.array(INTENSITY["Ratio"])[indices1]/average_e/np.array(FILE["Eff"])[indices2]*average_f, alpha=0.75, s=10)
plt.errorbar(convert_date(np.array(INTENSITY["Time"])[indices1], data_dates), np.array(INTENSITY["Ratio"])[indices1]/average_e/np.array(FILE["Eff"])[indices2]*average_f, yerr=np.array(INTENSITY["Err_ratio"])[indices1]/average_e/np.array(FILE["Eff"])[indices2]*average_f,fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
plt.xlabel('time', fontsize=40)
plt.ylabel('Relative efficiency', fontsize=40)
# plt.xlim([19000, 110000])
# plt.ylim([0.95, 1.08])
# plt.legend(fontsize=10)
plt.grid()
plt.title("Efficiency", fontsize=50)
#plt.legend(['experimental data (muons)', 'Monte-Carlo'], fontsize=30, loc=1)
plt.tick_params(axis='both', which='major', labelsize=30)

plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
         rotation_mode="anchor")

fig.savefig('relative_efficiency_norm.png')

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

INTENSITY = INTENSITY.iloc[indices1.tolist()]
FILE = FILE.iloc[indices2.tolist()]
INTENSITY['Eff_f'] = np.array(FILE["Eff"])

INTENSITY1 = INTENSITY[INTENSITY['Time']<54312]
INTENSITY1 = INTENSITY1[INTENSITY1['Time']>24776]
FILE1 = FILE[FILE['Time']>=24776]
FILE1 = FILE1[FILE1['Time']<=54312]
average_f1 = mean(FILE1["Eff"])                                       
average_e1 = mean(INTENSITY1["Ratio"])  
INTENSITY1.to_csv('intensity1.csv')

INTENSITY2 = INTENSITY[INTENSITY['Time']<87233]
INTENSITY2 = INTENSITY2[INTENSITY2['Time']>56627]
FILE2 = FILE[FILE['Time']>=56627]
FILE2 = FILE2[FILE2['Time']<=87233]
average_f2 = mean(FILE2["Eff"])                                       
average_e2 = mean(INTENSITY2["Ratio"]) 
INTENSITY2.to_csv('intensity2.csv')


INTENSITY3 = INTENSITY[INTENSITY['Time']<117657]
INTENSITY3 = INTENSITY3[INTENSITY3['Time']>89360]
FILE3 = FILE[FILE['Time']>=89360]
FILE3 = FILE3[FILE3['Time']<=117657]
average_f3 = mean(FILE3["Eff"])                                       
average_e3 = mean(INTENSITY3["Ratio"]) 
INTENSITY3.to_csv('intensity3.csv')


tims1 = np.array(INTENSITY3['Time'])
tims2 = np.array(FILE3['Time'])
indices1 = []
indices2 = []
for j, t in enumerate(tims1):
    find = False
    k = 0
    for i in range(len(tims2)):
        if (np.array(FILE3['Stop'])[i] > t and np.array(FILE3['Start'])[i] <= t) or (np.array(FILE3['Stop'])[i] >= t and np.array(FILE3['Start'])[i] < t):
            find = True
            k = i
    
    if not find:
        print(t, 'POLUNDRAAA')
    else:
        indices1.append(j)
        indices2.append(k)

indices1 = np.array(indices1)
indices2 = np.array(indices2)

# fig, ax = plt.subplots(figsize=(20,13))
# # graph of the ratio of detached decayed muons to stopped muons and efficiency (from file) versus conventional time (normalized)
# plt.scatter(convert_date(np.array(INTENSITY3["Time"])[indices1], data_dates), np.array(INTENSITY3["Ratio"])[indices1]/average_e3/np.array(FILE3["Eff"])[indices2]*average_f3, alpha=0.75, s=10)
# plt.errorbar(convert_date(np.array(INTENSITY3["Time"])[indices1], data_dates), np.array(INTENSITY3["Ratio"])[indices1]/average_e3/np.array(FILE3["Eff"])[indices2]*average_f3, yerr=np.array(INTENSITY3["Err_ratio"])[indices1]/average_e3/np.array(FILE3["Eff"])[indices2]*average_f3,fmt='none', marker='o', markersize=4, linestyle='none', ecolor='k', elinewidth=1, capsize=5, capthick=1)
# plt.xlabel('time', fontsize=40)
# plt.ylabel('Relative efficiency', fontsize=40)
# # plt.xlim([19000, 110000])
# # plt.ylim([0.95, 1.08])
# # plt.legend(fontsize=10)
# plt.grid()
# plt.title("Efficiency", fontsize=50)
# #plt.legend(['experimental data (muons)', 'Monte-Carlo'], fontsize=30, loc=1)
# plt.tick_params(axis='both', which='major', labelsize=30)

# plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
#          rotation_mode="anchor")

# fig.savefig('relative_efficiency3.png')

























#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

'''
colors = plt.cm.viridis(INTENSITY["Ratio"]/np.max(INTENSITY["Ratio"]))
colors = sns.color_palette(colors, as_cmap=True)

fig, axs = plt.subplots(figsize=(len(sums), 20), constrained_layout=True)
axs = sns.barplot(x="Time", y="Ratio", data=INTENSITY, palette=colors, dodge=False)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.ylim([0.145,0.165])
plt.grid(True)
fig.savefig('final_hist.pdf')
'''

'''
colors = plt.cm.viridis(INTENSITY["Select"]/np.max(INTENSITY["Select"]))
colors = sns.color_palette(colors, as_cmap=True)

fig, axs = plt.subplots(figsize=(len(sums), 20), constrained_layout=True)
axs = sns.barplot(x="Time", y="Select", data=INTENSITY, palette=colors, dodge=False)
plt.tick_params(axis='both', which='major', labelsize=20)
fig.savefig('select_hist.pdf')
'''

'''
fig, axs = plt.subplots(sharey=True, figsize=(40, 20))
INTENSITY = INTENSITY[INTENSITY["Ratio"] > 0.145]
#sns.boxplot(x=INTENSITY["Ratio"], flierprops={"marker": "x"}, boxprops={"facecolor": (.4, .6, .8, .5)}, medianprops={"color": "coral"}, linewidth=5)
sns.distplot(INTENSITY["Ratio"], hist=True, kde=False, bins=int(50), color = 'blue',hist_kws={'edgecolor':'black'})
plt.tick_params(axis='both', which='major', labelsize=30)
axs.grid(True)
fig.savefig('distribution_relation.pdf')
'''

# raf_file = "/home/clusters/02/n_skrobova/pavlyuk/reader_sumevents.txt"
# data = np.loadtxt(raf_file, dtype='U')
# raf_sums = data[::2].astype('int')
# raf_names = data[1::2]
# print(raf_sums)
# print(raf_names)