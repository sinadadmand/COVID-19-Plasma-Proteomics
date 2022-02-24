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


healthy=pd.read_excel(path_healthy)
severe= pd.read_excel(path_severe_i)
moderate= pd.read_excel(path_mild_i)
critical= pd.read_excel(path_critical_i)

#
# zo= healthy[["H1","H2","H3","H4","H5","H6","H7","H8","H9"]]
#
# s=zo.std(axis=1)
# m=zo.mean(axis=1)
#
# c=s/m
# print(len(healthy))
# healthy["std"]=c*100
# healthy=healthy[healthy["std"]<=50]
# healthy.reset_index(inplace=True)
# print(len(healthy))
# healthy=healthy[healthy.columns[1:]]
# print(healthy)
# print(healthy)
#
# zo= moderate[["M1_1","M1_2","M1_3","M2_1","M3_1"]]
#
# s=zo.std(axis=1)
# m=zo.mean(axis=1)
#
# c=s/m
# print(len(moderate))
# moderate["std"]=c*100
# moderate=moderate[moderate["std"]<=50]
# moderate.reset_index(inplace=True)
# print(len(moderate))
# moderate=moderate[moderate.columns[1:]]
#
# zo= severe[["S1_1","S1_2","S1_3","S2_1","S2_2","S2_3","S3_1","S3_2"]]
#
# s=zo.std(axis=1)
# m=zo.mean(axis=1)
#
# c=s/m
# print(len(severe))
# severe["std"]=c*100
# severe=severe[severe["std"]<=50]
# severe.reset_index(inplace=True)
# print(len(severe))
# severe=severe[severe.columns[1:]]
#
#
# zo= critical[["C1_1","C1_2","C1_3","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1","C5_2","C6_1","C6_2","C8_1","C8_2"]]
#
# s=zo.std(axis=1)
# m=zo.mean(axis=1)
#
# c=s/m
# print(len(critical))
# critical["std"]=c*100
# critical=critical[critical["std"]<=50]
# critical.reset_index(inplace=True)
# print(len(critical))
# critical=critical[critical.columns[1:]]
#



merged=pd.merge(healthy,moderate, on="Accession", how="outer")
merged=pd.merge(merged,severe, on="Accession", how="outer")
z=pd.merge(merged,critical, on="Accession", how="outer")

print(len(z))

z=z[["Accession","genes_1","genes_2","genes_3","genes_4","H1","H2","H3","H4","H5","M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2","S1_1","S1_2","S1_3","S2_1"
     ,"S2_2","S2_3","S3_1","S3_2","C1_1","C1_2","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1","C5_2",
     "C6_1","C6_2"]]

# z.columns=["Accession","H1","H2","H3","H4","H5","H6","M1_1","M1_2","M1_3","M2_1","M3_1","S1_1","S1_2","S1_3","S2_1"
#      ,"S2_2","S2_3","S3_1","S3_2","C1_1","C1_2","C1_3","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1","C5_2","C6_1","C6_2"
#      ,"C7_1","C7_2"]

z["genes_1"].fillna(z["genes_2"],inplace=True)
z["genes_1"].fillna(z["genes_3"],inplace=True)
z["genes_1"].fillna(z["genes_4"],inplace=True)

# z.to_excel("/home/a/Desktop/res.xlsx")
z= z.dropna(thresh=len(z.T)*(0.0))
z.fillna(0,inplace=True)
z.reset_index(inplace=True)
myPSM_filter_df=z[z.columns[1:]]
print(myPSM_filter_df)

# myPSM_filter_df["M1"]=(myPSM_filter_df["M1_1"]+myPSM_filter_df["M1_2"]+myPSM_filter_df["M1_3"])/3
myPSM_filter_df["M1"]=myPSM_filter_df[["M1_1","M1_2","M1_3"]].median(axis=1)


myPSM_filter_df["M2"]=(myPSM_filter_df["M2_1"])/1
myPSM_filter_df["M3"]=(myPSM_filter_df["M3_1"])/1
myPSM_filter_df["M4"]=myPSM_filter_df[["M4_1","M4_2"]].median(axis=1)


myPSM_filter_df["S1"]=myPSM_filter_df[["S1_1","S1_2","S1_3"]].median(axis=1)
myPSM_filter_df["S2"]=myPSM_filter_df[["S2_1","S2_2","S2_3"]].median(axis=1)
myPSM_filter_df["S3"]=myPSM_filter_df[["S3_1","S3_2"]].median(axis=1)


# myPSM_filter_df["S1"]=(myPSM_filter_df["S1_3"])/1
# myPSM_filter_df["S2"]=(myPSM_filter_df["S2_3"])/1
# myPSM_filter_df["S3"]=(myPSM_filter_df["S3_2"])/1


myPSM_filter_df["C1"]=myPSM_filter_df[["C1_1","C1_2"]].median(axis=1)
myPSM_filter_df["C2"]=myPSM_filter_df[["C2_1","C2_2","C2_3"]].median(axis=1)
myPSM_filter_df["C3"]=myPSM_filter_df[["C3_1","C3_2"]].median(axis=1)
myPSM_filter_df["C4"]=myPSM_filter_df[["C4_1","C4_2"]].median(axis=1)
myPSM_filter_df["C5"]=myPSM_filter_df[["C5_1","C5_2"]].median(axis=1)
myPSM_filter_df["C6"]=myPSM_filter_df[["C6_1","C6_2"]].median(axis=1)


myPSM_filter_df=myPSM_filter_df[["Accession","genes_1","H1","H2","H3","H4","H5","M1","M2","M3","M4","S1","S2","S3","C1","C2","C3","C4","C5",
                  "C6"]]


zo= myPSM_filter_df[["H1","H2","H3","H4","H5"]]

s=zo.std(axis=1)
m=zo.mean(axis=1)
#
c=s/m

myPSM_filter_df["std"]=c*100
myPSM_filter_df=myPSM_filter_df[(myPSM_filter_df["std"]<70)]
myPSM_filter_df.reset_index(inplace=True)
print(myPSM_filter_df)
myPSM_filter_df=myPSM_filter_df[myPSM_filter_df.columns[1:-1]]

# myPSM_filter_df.to_excel("/home/a/Desktop/pca_vars.xlsx")


print(len(myPSM_filter_df))
#
myPSM_filter_df=myPSM_filter_df[myPSM_filter_df.columns[2:]]

print(myPSM_filter_df)

mydf=myPSM_filter_df.T



X_std=StandardScaler().fit_transform(mydf)

#
pca=PCA(n_components=18)


principalComponents=pca.fit_transform(X_std)
print(principalComponents)

print(X_std)
features=range(pca.n_components)
#
#
#
PCA_components=pd.DataFrame(principalComponents)
#
component0= PCA_components[0].tolist()

component1= PCA_components[1].tolist()

fig,ax=plt.subplots()
#
ax.annotate("H1",(component0[0]+0.15,component1[0]+0.2),fontsize=20)
ax.annotate("H2",(component0[1]-0.45,component1[1]+0.45),fontsize=20)
ax.annotate("H3",(component0[2]-0.6,component1[2]+0.3),fontsize=20)
ax.annotate("H4",(component0[3]+0.15,component1[3]+0.2),fontsize=20)
ax.annotate("H5",(component0[4]+0.15,component1[4]+0.2),fontsize=20)
# ax.annotate("H6",(component0[5]-0.55,component1[5]-1.2),fontsize=20)
ax.annotate("P1",(component0[5]+0.15,component1[5]+0.2),fontsize=20)
ax.annotate("P2",(component0[6]+0.1,component1[6]+0.2),fontsize=20)
ax.annotate("P3",(component0[7]-2,component1[7]+0.5),fontsize=20)
ax.annotate("P4",(component0[8]-0.2,component1[8]+0.3),fontsize=20)
ax.annotate("P5",(component0[9]+0.15,component1[9]+0.1),fontsize=20)
ax.annotate("P6",(component0[10]-0.55,component1[10]+0.2),fontsize=20)
ax.annotate("P7",(component0[11]+0.45,component1[11]-0.5),fontsize=20)
ax.annotate("P8",(component0[12]+0.35,component1[12]+0.1),fontsize=20)
ax.annotate("P9",(component0[13]+0.1,component1[13]+0.2),fontsize=20)
ax.annotate("P10",(component0[14]+0.15,component1[14]+0.2),fontsize=20)
ax.annotate("P11",(component0[15]+0.15,component1[15]+0.2),fontsize=20)
ax.annotate("P12",(component0[16]+0.15,component1[16]+0.2),fontsize=20)
ax.annotate("P13",(component0[17]+0.25,component1[17]-0.1),fontsize=20)



plt.scatter([component0[0]],[component1[0]],marker='*',color="blue",label="Healthy")
plt.scatter([component0[1]],[component1[1]],marker='*',color="blue")
plt.scatter([component0[2]],[component1[2]],marker='*',color="blue")
plt.scatter([component0[3]],[component1[3]],marker='*',color="blue")
plt.scatter([component0[4]],[component1[4]],marker='*',color="blue")
# plt.scatter([component0[5]],[component1[5]],marker='*',color="blue")
# plt.scatter([component0[6]],[component1[6]],marker='*',color="blue")

# plt.scatter([component0[6]],[component1[6]],marker='s',color="darkturquoise",label='moderate')
# plt.scatter([component0[7]],[component1[7]],marker='s',color="darkturquoise")
# plt.scatter([component0[8]],[component1[8]],marker='s',color="darkturquoise")
# plt.scatter([component0[9]],[component1[9]],marker='s',color="darkturquoise")
# # plt.scatter([component0[10]],[component1[10]],marker='s',color="darkturquoise")
# # plt.scatter([component0[11]],[component1[11]],marker='s',color="darkturquoise")
# # plt.scatter([component0[12]],[component1[12]],marker='s',color="darkturquoise")
#
# plt.scatter([component0[10]],[component1[10]],marker='X',color="orange",label="severe")
# plt.scatter([component0[11]],[component1[11]],marker='X',color="orange")
# plt.scatter([component0[12]],[component1[12]],marker='X',color="orange")
# plt.scatter([component0[16]],[component1[16]],marker='X',color="orange")
# plt.scatter([component0[17]],[component1[17]],marker='X',color="orange")
# plt.scatter([component0[18]],[component1[18]],marker='X',color="orange")
# plt.scatter([component0[19]],[component1[19]],marker='X',color="orange")
# plt.scatter([component0[20]],[component1[20]],marker='X',color="orange")

# plt.scatter([component0[12]],[component1[12]],marker='^',color="brown",label="critical")
plt.scatter([component0[5]],[component1[5]],marker='^',color="brown",label="Patients")
plt.scatter([component0[6]],[component1[6]],marker='^',color="brown")
plt.scatter([component0[7]],[component1[7]],marker='^',color="brown")
plt.scatter([component0[8]],[component1[8]],marker='^',color="brown")
plt.scatter([component0[9]],[component1[9]],marker='^',color="brown")
plt.scatter([component0[10]],[component1[10]],marker='^',color="brown")
plt.scatter([component0[11]],[component1[11]],marker='^',color="brown")
plt.scatter([component0[12]],[component1[12]],marker='^',color="brown")
plt.scatter([component0[13]],[component1[13]],marker='^',color="brown")
plt.scatter([component0[14]],[component1[14]],marker='^',color="brown")
plt.scatter([component0[15]],[component1[15]],marker='^',color="brown")
plt.scatter([component0[16]],[component1[16]],marker='^',color="brown")
plt.scatter([component0[17]],[component1[17]],marker='^',color="brown")



# plt.scatter([component0[27]],[component1[27]],marker='^',color="brown")
# plt.scatter([component0[28]],[component1[28]],marker='^',color="brown")
# plt.scatter([component0[29]],[component1[29]],marker='^',color="brown")
# plt.scatter([component0[30]],[component1[30]],marker='^',color="brown")
# plt.scatter([component0[31]],[component1[31]],marker='^',color="brown")
# plt.scatter([component0[32]],[component1[32]],marker='^',color="brown")
# plt.scatter([component0[33]],[component1[33]],marker='^',color="brown")
# plt.scatter([component0[34]],[component1[34]],marker='^',color="brown")
# plt.scatter([component0[18]],[component1[18]],marker='^',color="brown")

plt.xlabel("Principal Component 1", fontsize=25)
plt.ylabel("Principal Component 2",fontsize=25)

plt.legend(prop={'size': 20})
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

plt.subplots_adjust(left=0.09,bottom=0.13,right=0.6)
plt.show()

# print(len(z))
print(len(myPSM_filter_df))
