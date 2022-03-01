import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
from functools import reduce
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statannot import add_stat_annotation
import matplotlib

# Let me filter the combined data of severe, healthy and critical, then merge the data!

path_healthy="/home/ali/Desktop/combined_data/combined_filtered_healthy.xlsx"
path_severe_i="/home/ali/Desktop/combined_data/combined_filtered_severe.xlsx"
path_mild_i="/home/ali/Desktop/combined_data/combined_filtered_moderate.xlsx"
path_critical_i="/home/ali/Desktop/combined_data/combined_filtered_critical.xlsx"
path_ms_recovery="/home/ali/Documents/projects/COVID_projects/Combined_data/combined_filtered_moderate_severe_recovery.xlsx"
path_c_recovery="/home/ali/Documents/projects/COVID_projects/Combined_data/combined_filtered_critical_recovery.xlsx"


healthy=pd.read_excel(path_healthy)
severe= pd.read_excel(path_severe_i)
moderate= pd.read_excel(path_mild_i)
critical= pd.read_excel(path_critical_i)
mode_severe=pd.read_excel(path_ms_recovery)
critical_rec=pd.read_excel(path_c_recovery)




merged=pd.merge(healthy,moderate, on="Accession", how="outer")
merged=pd.merge(merged,severe, on="Accession", how="outer")
merged=pd.merge(merged,critical, on="Accession", how="outer")
merged=pd.merge(merged,mode_severe, on="Accession", how="outer")
merged=pd.merge(merged,critical_rec, on="Accession", how="outer")



merged=merged[["Accession","genes_1","genes_2","genes_3","genes_4","genes_1R","genes1_RR","H1","H2","H3","H4","H5","M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2","S1_1","S1_2","S1_3","S2_1"
     ,"S2_2","S2_3","S3_1","S3_2","C1_1","C1_2","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1","C5_2",
     "C6_1","C6_2","C1RR","C2RR","C3RR","C4RR","C5RR","C1R","C2R","C3R","C4R","C5R"]]





z=merged
z["genes_1"].fillna(z["genes_2"],inplace=True)
z["genes_1"].fillna(z["genes_3"],inplace=True)
print(z["genes_4"])

z["genes_1"].fillna(z["genes_4"],inplace=True)
z["genes_1"].fillna(z["genes_1R"],inplace=True)
z["genes_1"].fillna(z["genes1_RR"],inplace=True)


myPSM_filter_df=z
myPSM_filter_df=myPSM_filter_df[["Accession","genes_1","H1","H2","H3","H4","H5","M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2","S1_1","S1_2","S1_3","S2_1"
     ,"S2_2","S2_3","S3_1","S3_2","C1_1","C1_2","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1","C5_2",
     "C6_1","C6_2","C1RR","C2RR","C3RR","C4RR","C5RR","C1R","C2R","C3R","C4R","C5R"]]



zo= myPSM_filter_df[["M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2","S1_1","S1_2","S1_3","S2_1"
     ,"S2_2","S2_3","S3_1","S3_2","C1_1","C1_2","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1","C5_2",
     "C6_1","C6_2","C1RR","C2RR","C3RR","C4RR","C5RR","C1R","C2R","C3R","C4R","C5R"]]

path_moderate="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2A/moderate_vs_healthy/moderate_healthy.xlsx"
path_severe="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2B/severe_vs_healthy/severe_healthy.xlsx"
path_critical="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2C/critical_vs_healthy/critical_healthy.xlsx"


moderate=pd.read_excel(path_moderate)
severe=pd.read_excel(path_severe)
critical=pd.read_excel(path_critical)

concated=pd.concat([moderate,severe,critical])

concated.reset_index(inplace=True)
concated=concated[concated.columns[2:]]


concated=list(set(concated["Accession"].to_list()))

print(len(concated))

myPSM_filter_df=myPSM_filter_df[myPSM_filter_df["Accession"].isin(concated)]


print(len(myPSM_filter_df))
myPSM_filter_df.to_excel("/home/ali/Desktop/result_bakis.xlsx")



myPSM_filter_df=myPSM_filter_df[myPSM_filter_df.columns[0:]]

print(myPSM_filter_df)



myPSM_filter_df["E1"]=myPSM_filter_df["M1_1"]
myPSM_filter_df["E2"]=myPSM_filter_df["M4_1"]
myPSM_filter_df["E3"]=myPSM_filter_df["S1_1"]
myPSM_filter_df["E4"]=myPSM_filter_df["S2_1"]
myPSM_filter_df["E5"]=myPSM_filter_df["C3_1"]
myPSM_filter_df["E6"]=myPSM_filter_df["C4_1"]
myPSM_filter_df["E7"]=myPSM_filter_df["C5_1"]
myPSM_filter_df["E8"]=myPSM_filter_df["C6_1"]


myPSM_filter_df["M1"]=(myPSM_filter_df["M1_2"]+myPSM_filter_df["M1_3"])/2
myPSM_filter_df["M2"]=myPSM_filter_df["M2_1"]
myPSM_filter_df["M3"]=myPSM_filter_df["M3_1"]
myPSM_filter_df["M4"]=myPSM_filter_df["M4_2"]
myPSM_filter_df["M5"]=(myPSM_filter_df["S1_2"]+myPSM_filter_df["S1_3"])/2 # S1
myPSM_filter_df["M6"]=(myPSM_filter_df["S2_2"]+myPSM_filter_df["S2_3"])/2# S2
myPSM_filter_df["M7"]=(myPSM_filter_df["S3_1"]+myPSM_filter_df["S3_2"])/2# S3
myPSM_filter_df["M8"]=(myPSM_filter_df["C1_1"]+myPSM_filter_df["C1_2"])/2# C1
myPSM_filter_df["M9"]=(myPSM_filter_df["C2_1"]+myPSM_filter_df["C2_2"]+myPSM_filter_df["C2_3"])/3# C2
myPSM_filter_df["M10"]=myPSM_filter_df["C3_2"]# C3
myPSM_filter_df["M11"]=myPSM_filter_df["C4_2"]# C4
myPSM_filter_df["M12"]=myPSM_filter_df["C5_2"]# C5
myPSM_filter_df["M13"]=myPSM_filter_df["C6_2"]# C6

myPSM_filter_df["R1"]=myPSM_filter_df["C1RR"]
myPSM_filter_df["R2"]=myPSM_filter_df["C2RR"]
myPSM_filter_df["R3"]=myPSM_filter_df["C3RR"]
myPSM_filter_df["R4"]=myPSM_filter_df["C4RR"]
myPSM_filter_df["R5"]=myPSM_filter_df["C5RR"]

myPSM_filter_df["R6"]=myPSM_filter_df["C1R"]
myPSM_filter_df["R7"]=myPSM_filter_df["C2R"]
myPSM_filter_df["R8"]=myPSM_filter_df["C3R"]
myPSM_filter_df["R9"]=myPSM_filter_df["C4R"]
myPSM_filter_df["R10"]=myPSM_filter_df["C5R"]


myPSM_filter_df=myPSM_filter_df[["Accession","E1","E2","E3","E4","E5","E6","E7","E8","M1","M2","M3","M4","M5","M6"
                    ,"M7","M8","M9","M10","M11","M12","M13","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10"]]

zo= myPSM_filter_df[["Accession","E1","E2","E3","E4","E5","E6","E7","E8","M1","M2","M3","M4","M5","M6"
                    ,"M7","M8","M9","M10","M11","M12","M13","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10"]]

s=zo[zo.columns[1:]].std(axis=1)
m=zo[zo.columns[1:]].mean(axis=1)
#
c=s/m

z= myPSM_filter_df.dropna(thresh=len(z.T)*(0.3))

z.fillna(0.01,inplace=True)
z.reset_index(inplace=True)
myPSM_filter_df=z[z.columns[1:]]


myPSM_filter_df.to_excel("/home/ali/Desktop/heatmap_after_80percent_filtration.xlsx")


myPSM_filter_df=myPSM_filter_df[myPSM_filter_df.columns[1:]]




mydf=myPSM_filter_df.T



X_std=StandardScaler().fit_transform(mydf)

#
pca=PCA(n_components=31)


principalComponents=pca.fit_transform(X_std)
print(principalComponents)

print(X_std)
features=range(pca.n_components)



PCA_components=pd.DataFrame(principalComponents)

component0= PCA_components[0].tolist()

component1= PCA_components[1].tolist()

fig,ax=plt.subplots()

ax.annotate("M1",(component0[0]+0.15,component1[0]+0.2),fontsize=10)
ax.annotate("M4",(component0[1]-0.4,component1[1]+0.45),fontsize=10)
ax.annotate("S1",(component0[2]-0.2,component1[2]+0.3),fontsize=10)
ax.annotate("S2",(component0[3]+0.15,component1[3]+0.2),fontsize=10)
ax.annotate("C3",(component0[4]+0.15,component1[4]+0.2),fontsize=10)
ax.annotate("C4",(component0[5]+0.15,component1[5]),fontsize=10)
ax.annotate("C5",(component0[6]+0.1,component1[6]+0.1),fontsize=10)
ax.annotate("C6",(component0[7]+0.1,component1[7]+0.2),fontsize=10)
ax.annotate("M1",(component0[8]+0.15,component1[8]+0.2),fontsize=10)
ax.annotate("M2",(component0[9]-0.4,component1[9]+0.45),fontsize=10)
ax.annotate("M3",(component0[10]-0.2,component1[10]+0.3),fontsize=10)
ax.annotate("M4",(component0[11]+0.15,component1[11]+0.2),fontsize=10)
ax.annotate("S1",(component0[12]+0.15,component1[12]+0.2),fontsize=10)
ax.annotate("S2",(component0[13]+0.15,component1[13]),fontsize=10)
ax.annotate("S3",(component0[14]+0.15,component1[14]+0.2),fontsize=10)
ax.annotate("C1",(component0[15]-0.7,component1[15]+0.2),fontsize=10)
ax.annotate("C2",(component0[16]+0.1,component1[16]+0.2),fontsize=10)
ax.annotate("C3",(component0[17]+0.3,component1[17]-0.3),fontsize=10)
ax.annotate("C4",(component0[18]+0.1,component1[18]+0.2),fontsize=10)
ax.annotate("C5",(component0[19]-0.6,component1[19]+0.2),fontsize=10)
ax.annotate("C6",(component0[20]+0.2,component1[20]-0.2),fontsize=10)

ax.annotate("M1",(component0[21]+0.1,component1[21]+0.2),fontsize=10)
ax.annotate("M4",(component0[22]+0.1,component1[22]+0.2),fontsize=10)
ax.annotate("S1",(component0[23]+0.1,component1[23]+0.2),fontsize=10)
ax.annotate("S2",(component0[24]+0.1,component1[24]+0.2),fontsize=10)
ax.annotate("S3",(component0[25]+0.1,component1[25]+0.2),fontsize=10)
ax.annotate("C1",(component0[26]+0.1,component1[26]+0.2),fontsize=10)
ax.annotate("C2",(component0[27]+0.1,component1[27]+0.2),fontsize=10)
ax.annotate("C4",(component0[28]+0.1,component1[28]+0.2),fontsize=10)
ax.annotate("C5",(component0[29]+0.1,component1[29]+0.2),fontsize=10)
ax.annotate("C6",(component0[30]+0.1,component1[30]+0.2),fontsize=10)

from matplotlib.patches import Patch

handles2=[Patch(facecolor="darkturquoise"),Patch(facecolor="orange"),Patch(facecolor="brown")]

plt.scatter([component0[0]],[component1[0]],marker='o',color="#DAC4F7",label="Early Infection")
plt.scatter([component0[1]],[component1[1]],marker='o',color="#DAC4F7")
plt.scatter([component0[2]],[component1[2]],marker='o',color="#DAC4F7")
plt.scatter([component0[3]],[component1[3]],marker='o',color="#DAC4F7")
plt.scatter([component0[4]],[component1[4]],marker='o',color="#DAC4F7")
plt.scatter([component0[5]],[component1[5]],marker='o',color="#DAC4F7")
plt.scatter([component0[6]],[component1[6]],marker='o',color="#DAC4F7")
plt.scatter([component0[7]],[component1[7]],marker='o',color="#DAC4F7")
#
#
plt.scatter([component0[8]],[component1[8]],marker='o',color="#F4989C",label="Maximum Infection")
plt.scatter([component0[9]],[component1[9]],marker='o',color="#F4989C")
plt.scatter([component0[10]],[component1[10]],marker='o',color="#F4989C")
plt.scatter([component0[11]],[component1[11]],marker='o',color="#F4989C")
plt.scatter([component0[12]],[component1[12]],marker='o',color="#F4989C")
plt.scatter([component0[13]],[component1[13]],marker='o',color="#F4989C")
plt.scatter([component0[14]],[component1[14]],marker='o',color="#F4989C")
plt.scatter([component0[15]],[component1[15]],marker='o',color="#F4989C")
plt.scatter([component0[16]],[component1[16]],marker='o',color="#F4989C")
plt.scatter([component0[17]],[component1[17]],marker='o',color="#F4989C")
plt.scatter([component0[18]],[component1[18]],marker='o',color="#F4989C")
plt.scatter([component0[19]],[component1[19]],marker='o',color="#F4989C")
plt.scatter([component0[20]],[component1[20]],marker='o',color="#F4989C")

plt.scatter([component0[21]],[component1[21]],marker='o',color="#D6F6DD",label="Post-Infection")
plt.scatter([component0[22]],[component1[22]],marker='o',color="#D6F6DD")
plt.scatter([component0[23]],[component1[23]],marker='o',color="#D6F6DD")
plt.scatter([component0[24]],[component1[24]],marker='o',color="#D6F6DD")
plt.scatter([component0[25]],[component1[25]],marker='o',color="#D6F6DD")
plt.scatter([component0[26]],[component1[26]],marker='o',color="#D6F6DD")
plt.scatter([component0[27]],[component1[27]],marker='o',color="#D6F6DD")
plt.scatter([component0[28]],[component1[28]],marker='o',color="#D6F6DD")
plt.scatter([component0[29]],[component1[29]],marker='o',color="#D6F6DD")
plt.scatter([component0[30]],[component1[30]],marker='o',color="#D6F6DD")


plt.xlabel("Principal Component 1", fontsize=20)
plt.ylabel("Principal Component 2",fontsize=20)

a= plt.legend(bbox_to_anchor=(0., 0.32, 1.51, 0.3),fontsize=15,title="Time",prop={'size': 15})
a.set_title('Time',prop={'size':15})


plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.subplots_adjust(left=0.12,bottom=0.13,right=0.6)
plt.show()
