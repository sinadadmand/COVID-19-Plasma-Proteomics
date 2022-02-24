import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats

path_healthy="/home/ali/Desktop/combined_data/combined_filtered_healthy.xlsx"
path_severe_infectious="/home/ali/Desktop/combined_data/combined_filtered_severe.xlsx"
path_moderate_infectious="/home/ali/Desktop/combined_data/combined_filtered_moderate.xlsx"
path_moderate_severe_recovery="/home/ali/Documents/projects/COVID_projects/Combined_data/combined_filtered_moderate_severe_recovery.xlsx"
path_critical_infectious="/home/ali/Desktop/combined_data/combined_filtered_critical.xlsx"
path_critical_recovery="/home/ali/Desktop/combined_data/combined_filtered_critical_recovery.xlsx"


healthy=pd.read_excel(path_healthy)
severe=pd.read_excel(path_severe_infectious)
moderate=pd.read_excel(path_moderate_infectious)
critical=pd.read_excel(path_critical_infectious)
severe_moderate_recovery=pd.read_excel(path_moderate_severe_recovery)
critical_recovery=pd.read_excel(path_critical_recovery)

merged=pd.merge(healthy,severe,how="outer",on="Accession")
merged=pd.merge(merged,moderate,how="outer",on="Accession")
merged=pd.merge(merged,critical,how="outer",on="Accession")
merged=pd.merge(merged,severe_moderate_recovery,how="outer",on="Accession")
merged=pd.merge(merged,critical_recovery,how="outer",on="Accession")
print(severe_moderate_recovery)


z=merged[["Accession","genes_1","genes_3","genes_2","genes_4","genes_1R","genes1_RR","H1","H2","H3","H4","H5",
                        "M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2",
                        "S1_1","S1_2","S1_3","S2_1","S2_2","S2_3","S3_1","S3_2",
                          "C1R","C2R","C3R","C4R","C5R","C1_1","C1_2",
                            "C2_1","C2_2","C2_3","C3_1","C3_2"
                            ,"C4_1","C4_2","C5_1","C5_2","C6_1", "C6_2",
                            "C3RR","C4RR","C5RR","C1RR","C2RR"]]


z["genes_1"].fillna(z["genes_2"],inplace=True)
z["genes_1"].fillna(z["genes_3"],inplace=True)
z["genes_1"].fillna(z["genes_4"],inplace=True)
z["genes_1"].fillna(z["genes_1R"],inplace=True)
z["genes_1"].fillna(z["genes1_RR"],inplace=True)


print(z)


z=z[["Accession","genes_1","H1","H2","H3","H4","H5",
                        "M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2",
                        "S1_1","S1_2","S1_3","S2_1","S2_2","S2_3","S3_1","S3_2",
                          "C1R","C2R","C3R","C4R","C5R","C1_1","C1_2",
                            "C2_1","C2_2","C2_3","C3_1","C3_2"
                            ,"C4_1","C4_2","C5_1","C5_2","C6_1", "C6_2",
                            "C3RR","C4RR","C5RR","C1RR","C2RR"]]


# Moderate-Severe Early --> M1_1, M4_1, S1_1, S2_1

# Critical Early --> C3_1, C4_1, C5_1, C6_1


# earlyms-> 6, 9, 11, 14

# latems ->  7, 8, 10, 12, 13, 15, 16



# print(z)

aa= pd.DataFrame({"a":["CLEC3B","MB","S100A9","ITIH2","MST1"]})

myPSM_filter_df=z[z["genes_1"].isin(aa["a"])]
#
myPSM_filter_df.fillna(0.0001,inplace=True)

print(myPSM_filter_df)


myPSM_filter_df.to_excel("/home/ali/Desktop/all_genes.xlsx")


which=4


# myPSM_filter_df.to_excel("/home/ali/Desktop/overall.xlsx")

average=myPSM_filter_df[["H1","H2","H3","H4","H5","M1_1","M4_1","S1_1","S2_1","C3_1","C4_1","C5_1","C6_1"]].mean(axis=1)
std=myPSM_filter_df[["H1","H2","H3","H4","H5","M1_1","M4_1","S1_1","S2_1","C3_1","C4_1","C5_1","C6_1"]].std(axis=1)
print(average)
print(std)


#### POST HOC TEST

healthy=myPSM_filter_df[["H1","H2","H3","H4","H5"]].iloc[0].to_list()


# max_num_h=max(healthy)
# healthy.remove(max_num_h)

ms=myPSM_filter_df[["M1_1","M4_1","S1_1","S2_1"]].iloc[0].to_list()
# max_num=max(ms)
# ms.remove(max_num)
# print(ms)

cr=myPSM_filter_df[["C3_1","C4_1","C5_1"]].iloc[0].to_list()
# min_num_cr=min(cr)
# cr.remove(min_num_cr)


total=healthy+ms+cr
samples=["Healthy"]*5+["Moderate-Severe"]*4+["Critical"]*3


print(len(healthy))
print(len(ms))
print(len(cr))

print(len(total))

print(len(samples))

my_Data=pd.DataFrame({"samples":samples,"total":total})
import scikit_posthocs as sp

# a= sp.posthoc(my_Data, val_col='total',group_col="samples")
# print(a)


healthy=myPSM_filter_df[["H1","H2","H3","H4","H5"]].iloc[4].to_list()
ms=myPSM_filter_df[["M1_1","M4_1","S1_1","S1_2"]].iloc[4].to_list()
cr=myPSM_filter_df[["C3_1","C4_1","C5_1","C6_1"]].iloc[4].to_list()


total=healthy+ms+cr

samples=["Healthy"]*5+["Moderate-Severe"]*4+["Critical"]*4


my_Data=pd.DataFrame({"samples":samples,"total":total})






import scikit_posthocs as sp


a= sp.posthoc_ttest(my_Data, val_col='total',group_col="samples")
print(a)

 

 # print(my_Data)

myPSM_filter_df["H1"]=(myPSM_filter_df["H1"]-average)/std
myPSM_filter_df["H2"]=(myPSM_filter_df["H2"]-average)/std
myPSM_filter_df["H3"]=(myPSM_filter_df["H3"]-average)/std
myPSM_filter_df["H4"]=(myPSM_filter_df["H4"]-average)/std
myPSM_filter_df["H5"]=(myPSM_filter_df["H5"]-average)/std
# myPSM_filter_df["H6"]=(myPSM_filter_df["H6"]-average)/std

myPSM_filter_df["M1_1"]=(myPSM_filter_df["M1_1"]-average)/std
myPSM_filter_df["M4_1"]=(myPSM_filter_df["M4_1"]-average)/std
myPSM_filter_df["S1_1"]=(myPSM_filter_df["S1_1"]-average)/std
myPSM_filter_df["S2_1"]=(myPSM_filter_df["S2_1"]-average)/std

myPSM_filter_df["C3_1"]=(myPSM_filter_df["C3_1"]-average)/std
myPSM_filter_df["C4_1"]=(myPSM_filter_df["C4_1"]-average)/std
myPSM_filter_df["C5_1"]=(myPSM_filter_df["C5_1"]-average)/std
myPSM_filter_df["C6_1"]=(myPSM_filter_df["C6_1"]-average)/std


# H5

# "C6_1"

healthy=myPSM_filter_df[["H1","H2","H3","H4","H5"]]

severe_moderate_early=myPSM_filter_df[["M1_1","M4_1","S1_1","S1_2"]]
critical_early=myPSM_filter_df[["C3_1","C4_1","C5_1"]]




healthy=healthy.iloc[which].to_list()
severe_moderate_early=severe_moderate_early.iloc[which].to_list()
critical_early=critical_early.iloc[which].to_list()

list_healthy=["Healthy"]*5
list_msearly=["Moderate-Severe Early"]*4
list_crearly=["Critical Early"]*3


list_of_x=list_healthy+list_msearly+list_crearly
numbers=healthy+severe_moderate_early+critical_early

# # import random
# # my_l=[]
# # for each in numbers:
# #     if each==0:
# #         num= random.uniform(0.0,0.5)
# #         my_l.append(num)
# #     else:
# #         my_l.append(each)
# # print(my_l)
# #
#
data_frame=pd.DataFrame({"log2(FC)":numbers,"Stages":list_of_x})
data_frame.to_excel("/home/ali/Desktop/only_gene.xlsx")

a,b=plt.subplots()
import matplotlib

font = {'family' : 'normal',
        'size'   : 20}

matplotlib.rc('font', **font)

#,palette={"Critical":"brown","Moderate-Severe":"blue"}

print(data_frame)
ax = sns.stripplot(x="Stages", y="log2(FC)",
                 data=data_frame,jitter=True, edgecolor="gray",palette={"Healthy":"gray","Critical Early":"brown","Moderate-Severe Early":"blue"})
# ax = sns.boxplot(x="Stages", y="log2(FC)",
#                  data=data_frame,palette={"Healthy":"gray","Critical Early":"brown","Moderate-Severe Early":"blue"})

plt.legend(bbox_to_anchor=(0., 0.48, 1.6, 0.3),fontsize=20)

from statannot import add_stat_annotation

a=["Early Infection","Late Infection","Post Infection"]
# add_stat_annotation(ax,data=data_frame,x="Stages", y="log2(FC)",hue="Patients",
#                     box_pairs=[(("Early","Critical"),("Early","Moderate-Severe")),
#                                (("Late","Critical"),("Late","Moderate-Severe")),
#                                (("Recovery","Critical"),("Recovery","Moderate-Severe"))],
#                     test='Mann-Whitney', text_format='star', loc='inside', verbose=2)


# import scipy
#
# first_analysis= scipy.stats.mannwhitneyu(A1,A3)
# print(first_analysis)
#
# for i in range(0,3):
#     f_ = data_frame[data_frame["Stages"] == a[i]]
#     first_early = f_[f_["Patients"] == "Moderate-Severe"]["log2(FC)"].to_list()
#     print(first_early)
#
#     s_ = data_frame[data_frame["Stages"] == a[i]]
#     second_early = s_[s_["Patients"] == "Critical"]["log2(FC)"].to_list()
#     print(second_early)
#     k = scipy.stats.kruskal(first_early, second_early)[1]
#     print(k)



# add_stat_annotation(ax,data=data_frame,x="Stages", y="log2(FC)",
#                     box_pairs=[("Healthy","Critical Early"),("Healthy","Moderate-Severe Early"),("Moderate-Severe Early","Critical Early")],
#                     test='Kruskal', text_format='star', loc='inside', verbose=2)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel("z-scores",fontsize=20)
ax.set(xlabel=None)


plt.text(0.9,4.5,"MB",fontsize=20)


# plt.text(0.5,1.8,"ns",fontsize=15)
# plt.text(1.48,2.1,"**",fontsize=15)
# plt.text(0.92,2.43,"**",fontsize=15)


# plt.axvline(x=0,ymin=0.76,ymax=0.79,color="black")
# plt.axvline(x=1,ymin=0.76,ymax=0.79,color="black")
# plt.axhline(y=1.71,xmin=0.167,xmax=0.5,color="black")
#
# plt.axvline(x=+1,ymin=0.83,ymax=0.86,color="black")
# plt.axvline(x=2,ymin=0.83,ymax=0.86,color="black")
# plt.axhline(y=2.05,xmin=0.5,xmax=0.833,color="black")
#
# plt.axvline(x=0.01,ymin=0.9,ymax=0.93,color="black")
# plt.axvline(x=+2,ymin=0.9,ymax=0.93,color="black")
# plt.axhline(y=2.38,xmin=0.17,xmax=0.833,color="black")


plt.legend(loc="upper center", bbox_to_anchor=(0.47,1.13), ncol= 2,frameon=False,fontsize=18)
plt.subplots_adjust(left=0.107,right=0.729,top=0.857)
# plt.text(0.95,3.9,"MB",fontsize=20)

plt.ylim(-2,4.3)

plt.show()
