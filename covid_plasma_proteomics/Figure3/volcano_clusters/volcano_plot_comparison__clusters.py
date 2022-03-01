import pandas as pd
from matplotlib.patches import Patch
import pandas as pd
import numpy as np
import pandas as pd
import numpy as np
from functools import reduce
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import mpl_toolkits
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar




path_critical_i="/home/ali/Desktop/combined_data/combined_filtered_critical.xlsx"
path_healthy="/home/ali/Desktop/combined_data/combined_filtered_healthy.xlsx"
path_severe_i="/home/ali/Desktop/combined_data/combined_filtered_severe.xlsx"
path_mild_i="/home/ali/Desktop/combined_data/combined_filtered_moderate.xlsx"
path_ms_rec="/home/ali/Documents/projects/COVID_projects/Combined_data/combined_filtered_moderate_severe_recovery.xlsx"
path_cr_rec="/home/ali/Documents/projects/COVID_projects/Combined_data/combined_filtered_critical_recovery.xlsx"


healty=pd.read_excel(path_healthy)
mild_i=pd.read_excel(path_mild_i)
severe_i=pd.read_excel(path_severe_i)
critical_i=pd.read_excel(path_critical_i)
ms=pd.read_excel(path_ms_rec)
cr=pd.read_excel(path_cr_rec)



merged=pd.merge(healty,mild_i, on="Accession", how="outer")
merged=pd.merge(merged,severe_i, on="Accession", how="outer")
merged=pd.merge(merged,critical_i, on="Accession", how="outer")
merged=pd.merge(merged,ms, on="Accession", how="outer")
merged=pd.merge(merged,cr, on="Accession", how="outer")

z=merged[["Accession","genes_1","genes_2","genes_3","genes_4","genes_1R","genes1_RR","H1","H2","H3","H4","H5","M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2","S1_1","S1_2","S1_3","S2_1"
     ,"S2_2","S2_3","S3_1","S3_2","C1_1","C1_2","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1","C5_2",
     "C6_1","C6_2","C1RR","C2RR","C3RR","C4RR","C5RR","C1R","C2R","C3R","C4R","C5R"]]

z["genes_1"].fillna(z["genes_2"],inplace=True)
z["genes_1"].fillna(z["genes_3"],inplace=True)
z["genes_1"].fillna(z["genes_4"],inplace=True)
z["genes_1"].fillna(z["genes_1R"],inplace=True)
z["genes_1"].fillna(z["genes1_RR"],inplace=True)



z=z[["Accession","genes_1","H1","H2","H3","H4","H5","M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2","S1_1","S1_2","S1_3","S2_1"
     ,"S2_2","S2_3","S3_1","S3_2","C1_1","C1_2","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1","C5_2",
         "C6_1","C6_2","C1RR","C2RR","C3RR","C4RR","C5RR","C1R","C2R","C3R","C4R","C5R"]]

myPSM_filter_df= z.dropna(thresh=0.40)

myPSM_filter_df.fillna(0.001,inplace=True)
myPSM_filter_df.reset_index(inplace=True)
myPSM_filter_df=myPSM_filter_df[myPSM_filter_df.columns[1:]]
print(myPSM_filter_df)

#
searcher_p="/home/ali/Documents/projects/COVID_projects/results/last_version-20211018T081343Z-001/last_version/critical_healthy/critical_healthy.xlsx"
searcher_m="/home/ali/Documents/projects/COVID_projects/results/last_version-20211018T081343Z-001/last_version/moderate_healthy/moderate_healthy.xlsx"
searcher_s="/home/ali/Documents/projects/COVID_projects/results/last_version-20211018T081343Z-001/last_version/severe_healthy/severe_healthy.xlsx"




myPSM_filter_df.reset_index(inplace=True)
myPSM_filter_df=myPSM_filter_df[myPSM_filter_df.columns[1:]]

myPSM_filter_df.to_excel("/home/ali/Desktop/pca_data.xlsx")


myPSM_filter_df["M1_1"]=myPSM_filter_df["M1_1"]
myPSM_filter_df["M1_2"]=myPSM_filter_df["M1_2"]
myPSM_filter_df["M1_3"]=myPSM_filter_df["M1_3"]
myPSM_filter_df["M2_1"]=myPSM_filter_df["M2_1"]
myPSM_filter_df["M3_1"]=myPSM_filter_df["M3_1"]

myPSM_filter_df["M4_1"]=myPSM_filter_df["M4_1"]
myPSM_filter_df["M4_2"]=myPSM_filter_df["M4_2"]


myPSM_filter_df["S1_1"]=myPSM_filter_df["S1_1"]
myPSM_filter_df["S1_2"]=myPSM_filter_df["S1_2"]
myPSM_filter_df["S1_3"]=myPSM_filter_df["S1_3"]
myPSM_filter_df["S2_1"]=myPSM_filter_df["S2_1"]
myPSM_filter_df["S2_2"]=myPSM_filter_df["S2_2"]
myPSM_filter_df["S2_3"]=myPSM_filter_df["S2_3"]
myPSM_filter_df["S3_1"]=myPSM_filter_df["S3_1"]
myPSM_filter_df["S3_2"]=myPSM_filter_df["S3_2"]

myPSM_filter_df["C1_1"]=myPSM_filter_df["C1_1"]
myPSM_filter_df["C1_2"]=myPSM_filter_df["C1_2"]
myPSM_filter_df["C2_1"]=myPSM_filter_df["C2_1"]
myPSM_filter_df["C2_2"]=myPSM_filter_df["C2_2"]
myPSM_filter_df["C2_3"]=myPSM_filter_df["C2_3"]
myPSM_filter_df["C3_1"]=myPSM_filter_df["C3_1"]
myPSM_filter_df["C3_2"]=myPSM_filter_df["C3_2"]
myPSM_filter_df["C4_1"]=myPSM_filter_df["C4_1"]
myPSM_filter_df["C4_2"]=myPSM_filter_df["C4_2"]
myPSM_filter_df["C5_1"]=myPSM_filter_df["C5_1"]
myPSM_filter_df["C5_2"]=myPSM_filter_df["C5_2"]
myPSM_filter_df["C6_1"]=myPSM_filter_df["C6_1"]
myPSM_filter_df["C6_2"]=myPSM_filter_df["C6_2"]

myPSM_filter_df["C1RR"]=myPSM_filter_df["C1RR"]
myPSM_filter_df["C2RR"]=myPSM_filter_df["C2RR"]
myPSM_filter_df["C3RR"]=myPSM_filter_df["C3RR"]
myPSM_filter_df["C4RR"]=myPSM_filter_df["C4RR"]
myPSM_filter_df["C5RR"]=myPSM_filter_df["C5RR"]

myPSM_filter_df["C1R"]=myPSM_filter_df["C1R"]
myPSM_filter_df["C2R"]=myPSM_filter_df["C2R"]
myPSM_filter_df["C3R"]=myPSM_filter_df["C3R"]
myPSM_filter_df["C4R"]=myPSM_filter_df["C4R"]
myPSM_filter_df["C5R"]=myPSM_filter_df["C5R"]


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


myPSM_filter_df=myPSM_filter_df[["genes_1","H1","H2","H3","H4","H5","E1","E2","E3","E4","E5","E6","E7","E8","M1","M2","M3","M4","M5","M6","M7","M8",
                                 "M9","M10","M11","M12","M13","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10"]]


healthy_samples= myPSM_filter_df[["H1","H2","H3","H4","H5"]]


zz= myPSM_filter_df


cluster1=zz[["R10","R9","E5","M10","R8","E7","M12","E6","M11","E8","M13"]]

cluster2=zz[["M9","R7","M8","R6","M2","M3","E4","R4","M7","R5","E2","M4","R2","R3","E3","E1","M1","M5","M6","R1"]]

a1_=cluster1.mean(axis=1)
a2_=cluster2.mean(axis=1)


tra_1=cluster1.T
tra_2=cluster2.T
import scipy

zz["p_value"]= scipy.stats.ttest_ind(tra_1,tra_2,equal_var=True)[1]




cl_mean1= cluster1[cluster1.columns].mean(axis=1)
cl_mean2= cluster2[cluster2.columns].mean(axis=1)

average= np.log2(cl_mean2/cl_mean1)

zz["average"]=pd.Series(average)




zz["t_test"]= scipy.stats.ttest_ind(tra_1,tra_2,equal_var=True,nan_policy="omit")[0]

zz["log 10 pvalue"]=-np.log10(zz["p_value"])



zz.loc[(zz['p_value']<0.05)&(zz['average']>0.3),'Sample downregulated (-) and upregulated (+) proteins']='+'
zz.loc[(zz['p_value']<0.05)&(zz['average']<-0.3),'Sample downregulated (-) and upregulated (+) proteins']='-'
zz['Sample downregulated (-) and upregulated (+) proteins']=zz['Sample downregulated (-) and upregulated (+) proteins'].fillna('Not Valid')

df1_Int_valid=zz[(zz['Sample downregulated (-) and upregulated (+) proteins']=='+')]
df1_Int_valid2=zz[zz['Sample downregulated (-) and upregulated (+) proteins']=='-']
df1_Int_notvalid=zz[zz['Sample downregulated (-) and upregulated (+) proteins']=='Not Valid']

df1_Int_valid.to_excel("/home/ali/Desktop/cluster1_specific.xlsx")
df1_Int_valid.reset_index(inplace=True)
df1_Int_valid2.reset_index(inplace=True)
df1_Int_notvalid.reset_index(inplace=True)

df1_Int_valid2.to_excel("/home/ali/Desktop/cluster2_specific.xlsx")


highest=df1_Int_valid.sort_values(by="average",ascending=False)

highest.reset_index(inplace=True)
highest=highest[highest.columns[1:]]
print(highest["average"])

highest.to_excel("/home/ali/Desktop/cluster1_specific.xlsx")

print(highest)


lowest=df1_Int_valid2.sort_values(by="average",ascending=True)

first_10={}
for j in range(1,100):
    x_axis=highest.iloc[j][39]
    y_axis=highest.iloc[j][41]
    first_10[highest.iloc[j][1]]=[x_axis,y_axis]

print(first_10)

least7={}
for j in range(0,25):
    x_axis=lowest.iloc[j][39]
    y_axis=lowest.iloc[j][41]
    least7[lowest.iloc[j][1]]=[x_axis,y_axis]

fig,ax=plt.subplots()

plt.plot('average', 'log 10 pvalue',data=df1_Int_valid,linestyle='',marker='o', markersize=5.5,alpha=0.50,color='slateblue',label="Significantly higher abundant in Cluster 2 (n=106)")
plt.plot('average', 'log 10 pvalue',data=df1_Int_valid2,linestyle='',marker='o', markersize=5.5,alpha=0.50,color='green',label="Significantly higher abundant in Cluster 1 (n=28)")
plt.plot('average', 'log 10 pvalue',data=df1_Int_notvalid,linestyle='',marker='o', markersize=3.5,alpha=0.40,color='gray',fillstyle= 'bottom',label="Insignificant between clusters")



plt.ylabel("-log 10 pvalue",fontsize=20)
plt.xlabel("log2(cluster2/cluster1)",fontsize=20)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.xlim(-4,5)

ax.annotate("FN1",(3.96417885048086+0.1,3.25511087446691),fontsize=10)
ax.annotate(list(first_10.keys())[0],(list(first_10.values())[0][0]+0.1,list(first_10.values())[0][1]),fontsize=10)
ax.annotate(list(first_10.keys())[1],(list(first_10.values())[1][0]+0.1,list(first_10.values())[1][1]+0.1),fontsize=10)
ax.annotate(list(first_10.keys())[2],(list(first_10.values())[2][0]+0.1,list(first_10.values())[2][1]),fontsize=10)
ax.annotate(list(first_10.keys())[3],(list(first_10.values())[3][0]+0.1,list(first_10.values())[3][1]+0.1),fontsize=10)
ax.annotate(list(first_10.keys())[4],(list(first_10.values())[4][0]+0.25,list(first_10.values())[4][1]+0.1),fontsize=10)
ax.annotate(list(first_10.keys())[5],(list(first_10.values())[5][0]-0.3,list(first_10.values())[5][1]+0.35),fontsize=10)
ax.annotate(list(first_10.keys())[6],(list(first_10.values())[6][0]+0.1,list(first_10.values())[6][1]),fontsize=10)
ax.annotate(list(first_10.keys())[7],(list(first_10.values())[7][0]+0.15,list(first_10.values())[7][1]),fontsize=10)
ax.annotate(list(first_10.keys())[9],(list(first_10.values())[9][0]+0.15,list(first_10.values())[9][1]),fontsize=10)
ax.annotate(list(first_10.keys())[10],(list(first_10.values())[10][0],list(first_10.values())[10][1]-0.35),fontsize=10)
ax.annotate(list(first_10.keys())[11],(list(first_10.values())[11][0]+0.05,list(first_10.values())[11][1]),fontsize=10)
ax.annotate(list(first_10.keys())[15],(list(first_10.values())[15][0]+0.15,list(first_10.values())[15][1]+1.1),fontsize=10)
ax.annotate(list(first_10.keys())[17],(list(first_10.values())[17][0]-0.15,list(first_10.values())[17][1]+0.17),fontsize=10)

ax.annotate(list(first_10.keys())[18],(list(first_10.values())[18][0]-0.02,list(first_10.values())[18][1]+0.15),fontsize=10)

ax.annotate(list(first_10.keys())[21],(list(first_10.values())[21][0]+0.1,list(first_10.values())[21][1]-0.2),fontsize=10)
ax.annotate(list(first_10.keys())[22],(list(first_10.values())[22][0]+0.1,list(first_10.values())[22][1]),fontsize=10)
ax.annotate(list(first_10.keys())[23],(list(first_10.values())[23][0]-0.15,list(first_10.values())[23][1]+0.23),fontsize=10)
ax.annotate(list(first_10.keys())[24],(list(first_10.values())[24][0],list(first_10.values())[24][1]-0.7),fontsize=10)
ax.annotate(list(first_10.keys())[27],(list(first_10.values())[27][0],list(first_10.values())[27][1]+0.05),fontsize=10)
ax.annotate(list(first_10.keys())[32],(list(first_10.values())[32][0]+0.05,list(first_10.values())[32][1]+0.05),fontsize=10)
ax.annotate(list(first_10.keys())[35],(list(first_10.values())[35][0]-0.1,list(first_10.values())[35][1]+0.1),fontsize=10)
ax.annotate(list(first_10.keys())[39],(list(first_10.values())[39][0]-0.1,list(first_10.values())[39][1]+0.1),fontsize=10)
ax.annotate(list(first_10.keys())[41],(list(first_10.values())[41][0],list(first_10.values())[41][1]-0.3),fontsize=10)
ax.annotate(list(first_10.keys())[48],(list(first_10.values())[48][0],list(first_10.values())[48][1]-1.3),fontsize=10)
ax.annotate(list(first_10.keys())[50],(list(first_10.values())[50][0]-0.1,list(first_10.values())[50][1]+0.15),fontsize=10)
ax.annotate(list(first_10.keys())[56],(list(first_10.values())[56][0]+0.1,list(first_10.values())[56][1]-0.02),fontsize=10)
ax.annotate(list(first_10.keys())[58],(list(first_10.values())[58][0]+0.1,list(first_10.values())[58][1]-0.05),fontsize=10)
ax.annotate(list(first_10.keys())[69],(list(first_10.values())[69][0]+0.1,list(first_10.values())[69][1]-0.05),fontsize=10)
ax.annotate(list(first_10.keys())[70],(list(first_10.values())[70][0]-0.3,list(first_10.values())[70][1]+0.1),fontsize=10)
ax.annotate(list(first_10.keys())[84],(list(first_10.values())[84][0]-0.1,list(first_10.values())[84][1]+0.1),fontsize=10)
ax.annotate(list(first_10.keys())[84],(list(first_10.values())[84][0]-0.1,list(first_10.values())[84][1]+0.1),fontsize=10)
ax.annotate(list(first_10.keys())[91],(list(first_10.values())[91][0],list(first_10.values())[91][1]),fontsize=10)


ax.annotate(list(least7.keys())[0],(list(least7.values())[0][0]-0.35,list(least7.values())[0][1]+0.09),fontsize=10)
ax.annotate(list(least7.keys())[1],(list(least7.values())[1][0]-0.25,list(least7.values())[1][1]+0.1),fontsize=10)
ax.annotate(list(least7.keys())[2],(list(least7.values())[2][0]-0.25,list(least7.values())[2][1]+0.07),fontsize=10)
ax.annotate(list(least7.keys())[3],(list(least7.values())[3][0]-0.25,list(least7.values())[3][1]+0.07),fontsize=10)
ax.annotate(list(least7.keys())[4],(list(least7.values())[4][0]-0.15,list(least7.values())[4][1]+0.32),fontsize=10)
ax.annotate(list(least7.keys())[5],(list(least7.values())[5][0]-0.15,list(least7.values())[5][1]+0.07),fontsize=10)
ax.annotate(list(least7.keys())[6],(list(least7.values())[6][0]-0.2,list(least7.values())[6][1]+0.07),fontsize=10)
ax.annotate(list(least7.keys())[7],(list(least7.values())[7][0]-0.8,list(least7.values())[7][1]),fontsize=10)
ax.annotate(list(least7.keys())[8],(list(least7.values())[8][0]+0.05,list(least7.values())[8][1]),fontsize=10)
ax.annotate(list(least7.keys())[9],(list(least7.values())[9][0],list(least7.values())[9][1]+0.15),fontsize=10)
ax.annotate(list(least7.keys())[10],(list(least7.values())[10][0]-0.1,list(least7.values())[10][1]-0.4),fontsize=10)
# ax.annotate(list(least7.keys())[11],(list(least7.values())[11][0],list(least7.values())[11][1]+0.1),fontsize=10)
ax.annotate(list(least7.keys())[12],(list(least7.values())[12][0]+0.1,list(least7.values())[12][1]),fontsize=10)
ax.annotate(list(least7.keys())[15],(list(least7.values())[15][0]-0.01,list(least7.values())[15][1]+0.1),fontsize=10)
# ax.annotate(list(least7.keys())[15],(list(least7.values())[15][0],list(least7.values())[15][1]+0.1),fontsize=10)
ax.annotate(list(least7.keys())[19],(list(least7.values())[19][0],list(least7.values())[19][1]+0.18),fontsize=10)
# ax.annotate(list(least7.keys())[21],(list(least7.values())[21][0],list(least7.values())[21][1]-0.15),fontsize=10)
ax.annotate(list(least7.keys())[23],(list(least7.values())[23][0]-0.4,list(least7.values())[23][1]+0.15),fontsize=10)


plt.subplots_adjust(top=0.95,right=0.607,left=0.064)

plt.legend(bbox_to_anchor=(0., 0.32, 1.71, 0.3),fontsize=12)

plt.show()
