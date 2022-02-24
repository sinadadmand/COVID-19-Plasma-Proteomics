import pandas as pd
import numpy as np
from functools import reduce
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
my_pathlist=["/home/ali/Desktop/combined_data/combined_filtered_healthy.xlsx",
             "/home/ali/Documents/projects/COVID_projects/Combined_data/combined_filtered_moderate.xlsx",
             "/home/ali/Documents/projects/COVID_projects/Combined_data/combined_filtered_severe.xlsx",
             "/home/ali/Desktop/combined_data/combined_filtered_critical.xlsx",
             "/home/ali/Documents/projects/COVID_projects/Combined_data/combined_filtered_moderate_severe_recovery.xlsx",
             "/home/ali/Desktop/combined_data/combined_filtered_critical_recovery.xlsx"
             ]

# healty=pd.read_excel(my_pathlist[0])
# critical=pd.read_excel(my_pathlist[1])
# recovery=pd.read_excel(my_pathlist[2])

healty=pd.read_excel(my_pathlist[0])
moderate=pd.read_excel(my_pathlist[1])
severe=pd.read_excel(my_pathlist[2])
recovery=pd.read_excel(my_pathlist[4])
recovery_c=pd.read_excel(my_pathlist[5])
critical=pd.read_excel(my_pathlist[3])



merged=pd.merge(healty,critical,how="outer",on="Accession")
merged=pd.merge(merged,recovery_c,how="outer",on="Accession")


z=merged[["Accession","genes_1R","genes_1","genes_2","H1","H2","H3","H4","H5",
                        "C1_1","C1_2","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2",
                          "C5_1","C5_2","C6_1","C6_2",
                         "C1R","C2R",
                         "C3R","C4R","C5R"]]

z["genes_1"].fillna(z["genes_2"],inplace=True)
# z["genes_1"].fillna(z["genes1R"],inplace=True)
z["genes_1"].fillna(z["genes_1R"],inplace=True)
# z["genes_1"].fillna(z["genes_3"],inplace=True)
# z["genes_1"].fillna(z["genes_4"],inplace=True)

z=z[["Accession","genes_1","H1","H2","H3","H4","H5",
                        "C1_1","C1_2","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2",
                        "C5_1","C5_2","C6_1","C6_2",
                         "C1R","C2R",
                         "C3R","C4R","C5R"]]



z.to_excel("/home/ali/Desktop/convert_it.xlsx")
myPSM_filter_df= z
myPSM_filter_df.fillna(0.0001,inplace=True)
myPSM_filter_df.reset_index(inplace=True)
z=myPSM_filter_df[myPSM_filter_df.columns[1:]]

list_of_proteins=["CLEC3B","MST1","MB","S100A9","ITIH2"]



print(len(list_of_proteins))

z=z[z["genes_1"].isin(list_of_proteins)]

print(z)

print(len(z))

z["Healthy"]=z[["H1","H2","H3","H4","H5"]].mean(axis=1)
z["Early Critical"]=z[["C3_1","C4_1",
                        "C5_1","C6_1"]].mean(axis=1)

z["Late Critical"]=z[["C1_1","C1_2","C2_1","C2_2","C2_3",
                        "C3_2","C4_2",
                        "C5_2","C6_2"]].mean(axis=1)

z["Recovery"]=z[["C1R","C2R","C3R","C4R","C5R"]].mean(axis=1)

z["Early"]=np.log2(z["Early Critical"]/z["Healthy"])
z["Late"]=np.log2(z["Late Critical"]/z["Healthy"])
z["Recover"]=np.log2(z["Recovery"]/z["Healthy"])


names= z["genes_1"].to_list()

z.to_excel("/home/ali/Desktop/birbakacanmi.xlsx")


earl= z["Early"].to_list()
lat=z["Late"].to_list()
rec=z["Recover"].to_list()

fig,ax=plt.subplots(figsize=(10,6))

# print(names)

print(z)

#
#
print(len(earl))
ax.plot(["Uninfected", "Early", "Inflammatory", "Recovery"], [0,earl[0], lat[0], rec[0]], color="darkgreen", alpha=1)
ax.plot(["Uninfected", "Early", "Inflammatory", "Recovery"], [0,earl[1], lat[1], rec[1]], color="blue", alpha=1)
ax.plot(["Uninfected", "Early", "Inflammatory", "Recovery"], [0,earl[2], lat[2], rec[2]], color="brown", alpha=1)
ax.plot(["Uninfected", "Early", "Inflammatory", "Recovery"], [0,earl[3], lat[3], rec[3]], color="orange", alpha=1)
ax.plot(["Uninfected", "Early", "Inflammatory", "Recovery"], [0,earl[4], lat[4], rec[4]], color="m", alpha=1)




plt.ylabel("log2 Relative Fold Change",fontsize=20)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)


ax.annotate(names[0],(3.1,rec[0]),fontsize=14,color="darkgreen")
ax.annotate(names[1],(3.1,rec[1]),fontsize=14,color="blue")
ax.annotate(names[2],(3.1,rec[2]),fontsize=14,color="brown")
ax.annotate(names[3],(3.1,rec[3]-0.2),fontsize=14,color="orange")
ax.annotate(names[4],(3.1,rec[4]+0.2),fontsize=14,color="m")
# ax.annotate(names[5],(3.1,rec[5]),fontsize=14,color="gray")
# ax.annotate(names[6],(3.1,rec[6]+0.1),fontsize=14,color="purple")
# ax.annotate(names[7],(3.1,rec[7]-0.15),fontsize=14,color="black")
# ax.annotate(names[8],(3.1,rec[8]+0.1),fontsize=14,color="black")
# ax.annotate(names[9],(3.1,rec[9]),fontsize=14,color="darkgreen")
# ax.annotate(names[10],(3.1,rec[10]+0.15),fontsize=14,color="brown")
# ax.annotate(names[11],(3.1,rec[11]),fontsize=14,color="indigo")
# ax.annotate(names[12],(3.1,rec[12]+0.2),fontsize=14,color="darkblue")
# ax.annotate(names[13],(3.1,rec[13]),fontsize=14,color="steelblue")

plt.xlim(0,3.5)
plt.ylim(-3,7)
# plt.legend(bbox_to_anchor=(0., 0.47, 1.3, 0.3),fontsize=20)

plt.show()
