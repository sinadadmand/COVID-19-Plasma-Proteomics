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
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar



path_critical_i="/home/ali/Desktop/combined_data/combined_filtered_critical.xlsx"
path_healthy="/home/ali/Desktop/combined_data/combined_filtered_healthy.xlsx"
path_severe_i="/home/ali/Desktop/combined_data/combined_filtered_severe.xlsx"
path_mild_i="/home/ali/Desktop/combined_data/combined_filtered_moderate.xlsx"

healty=pd.read_excel(path_healthy)
mild_i=pd.read_excel(path_mild_i)
severe_i=pd.read_excel(path_severe_i)
critical_i=pd.read_excel(path_critical_i)

merged=pd.merge(healty,mild_i, on="Accession", how="outer")
merged=pd.merge(merged,severe_i, on="Accession", how="outer")
merged=pd.merge(merged,critical_i, on="Accession", how="outer")

z=merged[["Accession","genes_1","genes_2","genes_3","genes_4","H1","H2","H3","H4","H5",
                         "M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2","S1_1","S1_2","S1_3",
                         "S2_1","S2_2","S2_3","S3_1","S3_2","C1_1","C1_2","C2_1","C2_2","C2_3",
                         "C3_1","C3_2","C4_1","C4_2","C5_1","C5_2","C6_1","C6_2"]]

z["genes_1"].fillna(z["genes_2"],inplace=True)
z["genes_1"].fillna(z["genes_3"],inplace=True)
z["genes_1"].fillna(z["genes_4"],inplace=True)

z=z[["Accession","genes_1","H1","H2","H3","H4","H5",
                         "M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2","S1_1","S1_2","S1_3",
                         "S2_1","S2_2","S2_3","S3_1","S3_2","C4_1","C4_2","C5_1","C5_2","C6_1","C6_2","C1_1","C1_2","C2_1","C2_2","C2_3",
                         "C3_1","C3_2"]]

myPSM_filter_df= z.dropna(thresh=0.40)

myPSM_filter_df.fillna(0.0001,inplace=True)
myPSM_filter_df.reset_index(inplace=True)
myPSM_filter_df=myPSM_filter_df[myPSM_filter_df.columns[1:]]
print(myPSM_filter_df)


searcher_p="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2C/critical_vs_healthy/critical_healthy.xlsx"
# second_p="/home/a/Desktop/covid_project_3_July/critical_severe/critical_severe.xlsx"
# searche="/home/a/Desktop/regulation/critical/upregulated_ones.xlsx"
# searcher_p="/home/a/Desktop/critical.xlsx"
search=pd.read_excel(searcher_p)
# second=pd.read_excel(second_p)

# search=search[(search["p_values"]<=0.01)]


# second=second[(second["log2(Mean(First/Second))"]>0.75)|(second["log2(Mean(First/Second))"]<-0.75)]

# second=second[second["p_values"]<0.04]
# search=search[(search["Accession"].isin(second["Accession"]))]
#
# search=search[(search["Accession"]!="P00738")&(search["Accession"]!="P0C0L5")]

# sec=pd.read_excel(second_p)
#print(search)
# cuttoff5=["VWF","S100A9","QSOX1","IHH","IGKV4-1","IGHV5-10-1","IGHV3-20","IGHV2-70D","IGFBP2","F9","CD163"]
# cuttoff4=['ITIH3','IGHD','IGKV1-27','IGKV4-1','IGFBP2','IGHV3-20','VWF','CST3','S100A9','IGHV2-70D','PKM',
#           'MB','IHH','AZGP1','F9','IGHV1-24']
#
# print(myPSM_filter_df)
# sec=sec[(sec["log2(Mean(First/Second))"]>0.75)|(sec["log2(Mean(First/Second))"]<-0.75)]

myPSM_filter_df=myPSM_filter_df[(myPSM_filter_df["Accession"].isin(search["Accession"]))]
myPSM_filter_df.reset_index(inplace=True)
myPSM_filter_df=myPSM_filter_df[myPSM_filter_df.columns[1:]]
print(myPSM_filter_df)
# sec.to_excel("/home/a/Desktop/nene.xlsx")
# print(len(myPSM_filter_df))
#
avg=(myPSM_filter_df["H1"]+myPSM_filter_df["H2"]+myPSM_filter_df["H3"]+myPSM_filter_df["H4"]+myPSM_filter_df["H5"])/5

myPSM_filter_df["M1_1"]=myPSM_filter_df["M1_1"]/avg
myPSM_filter_df["M1_2"]=myPSM_filter_df["M1_2"]/avg
myPSM_filter_df["M1_3"]=myPSM_filter_df["M1_3"]/avg
myPSM_filter_df["M2_1"]=myPSM_filter_df["M2_1"]/avg
myPSM_filter_df["M3_1"]=myPSM_filter_df["M3_1"]/avg
myPSM_filter_df["M4_1"]=myPSM_filter_df["M4_1"]/avg
myPSM_filter_df["M4_2"]=myPSM_filter_df["M4_2"]/avg


myPSM_filter_df["S1_1"]=myPSM_filter_df["S1_1"]/avg
myPSM_filter_df["S1_2"]=myPSM_filter_df["S1_2"]/avg
myPSM_filter_df["S1_3"]=myPSM_filter_df["S1_3"]/avg
myPSM_filter_df["S2_1"]=myPSM_filter_df["S2_1"]/avg
myPSM_filter_df["S2_2"]=myPSM_filter_df["S2_2"]/avg
myPSM_filter_df["S2_3"]=myPSM_filter_df["S2_3"]/avg
myPSM_filter_df["S3_1"]=myPSM_filter_df["S3_1"]/avg
myPSM_filter_df["S3_2"]=myPSM_filter_df["S3_2"]/avg
#
myPSM_filter_df["C1_1"]=myPSM_filter_df["C1_1"]/avg
myPSM_filter_df["C1_2"]=myPSM_filter_df["C1_2"]/avg
# myPSM_filter_df["C1_3"]=myPSM_filter_df["C1_3"]/avg
myPSM_filter_df["C2_1"]=myPSM_filter_df["C2_1"]/avg
myPSM_filter_df["C2_2"]=myPSM_filter_df["C2_2"]/avg
myPSM_filter_df["C2_3"]=myPSM_filter_df["C2_3"]/avg
myPSM_filter_df["C3_1"]=myPSM_filter_df["C3_1"]/avg
myPSM_filter_df["C3_2"]=myPSM_filter_df["C3_2"]/avg
myPSM_filter_df["C4_1"]=myPSM_filter_df["C4_1"]/avg
myPSM_filter_df["C4_2"]=myPSM_filter_df["C4_2"]/avg
myPSM_filter_df["C5_1"]=myPSM_filter_df["C5_1"]/avg
myPSM_filter_df["C5_2"]=myPSM_filter_df["C5_2"]/avg
myPSM_filter_df["C6_1"]=myPSM_filter_df["C6_1"]/avg
myPSM_filter_df["C6_2"]=myPSM_filter_df["C6_2"]/avg

# myPSM_filter_df["V1"]=myPSM_filter_df["V1"]/avg
# myPSM_filter_df["V2"]=myPSM_filter_df["V2"]/avg
# myPSM_filter_df["V3"]=myPSM_filter_df["V3"]/avg
# myPSM_filter_df["V4"]=myPSM_filter_df["V4"]/avg
#
genes=myPSM_filter_df["genes_1"].tolist()

zz= myPSM_filter_df[myPSM_filter_df.columns[8:]]

print(zz)
zz=np.log2(zz)
print(search)
print(zz)

# genes= zz["genes"].tolist()

# zz=zz[zz.columns[0:-2]]
print(zz)




pallete=sns.color_palette()

color_dict={}
color_dict2={}
color_dict3={}

# for col in zz.columns:
#     if col=="SI_1" or col=="SI_2" or col=="SI_3"or col=="SI_4"or col=="SI_5"or col=="SI_6":
#         color_dict3[col]="darkturquoise"
#     elif col=="SI_1" or col=="SI_2" or col=="SI_3"or col=="SI_4" or col=="SI_5"or col=="SI_6" or col=="SI_7"or col=="SI_8":
#         color_dict3[col]="orange"
#     else:
#         color_dict3[col]="brown"

# #
for col in zz.columns:
    # if  col=="H1" or col=="H2" or col=="H3":
    #     color_dict2[col]="lightgreen"
    if col=="M1_1" or col=="M1_2" or col=="M1_3"or col=="M2_1"or col=="M3_1"or col=="M4_1"or col=="M4_2":
        color_dict2[col]="darkturquoise"
    elif col=="S1_1" or col=="S1_2" or col=="S1_3"or col=="S2_1" or col=="S2_2"or col=="S2_3" or \
            col=="S3_1"or col=="S3_2":
        color_dict2[col]="orange"
    # elif col=="V1" or col=="V2" or col=="V3" or col=="V4":
    #     color_dict2[col] = "green"
    else:
        color_dict2[col]="brown"

color_col2=pd.Series(color_dict2)
for col in zz.columns:
    if col=="M3_1"or col=="S1_1"or col=="S1_2"or col=="S1_3" or col=="S3_1"or col=="S3_2" or col=="C1_1"or col=="C1_2"or col=="C1_3"or col=="C3_1"or col=="C3_2":
        color_dict[col]="navy"
    else:
        color_dict[col] = "crimson"

color_col=pd.Series(color_dict)
print(color_col)
print(color_dict2)
color_col2= pd.Series(color_dict2)

k= sns.clustermap(zz,method="complete",metric="cityblock",cmap="coolwarm",row_cluster=True,col_cluster=False,col_colors=[color_col2,color_col],yticklabels=genes,vmin=-6,vmax=6,xticklabels=False,colors_ratio=0.015,cbar_kws={"label":"log2(FC)"})#

# k.cax.set_visible(False)



pp=["darkturquoise","orange","brown"]
counter=0
for i in ["Moderate","Severe","Critical"]:
    k.ax_col_dendrogram.bar(0,0,color=pp[counter],label=i)
    counter+=1

s=k.ax_col_dendrogram.legend(title='Type', loc="upper right", ncol=1, bbox_to_anchor=(0.4,0.75), bbox_transform=plt.gcf().transFigure,prop={'size':11},title_fontsize=11)


plt.setp(k.ax_heatmap.yaxis.get_majorticklabels(),rotation=0,fontsize=6,fontstyle="normal")
plt.subplots_adjust(top=1.24,right=0.5,wspace=0,bottom=0.02)

k.cax.set_position([0.63,.22,0.02,0.15])

my_color1=[pallete[0],(52, 205, 235)]
# my_color2=[pallete[3],pallete[4],pallete[5],pallete[6]]

handles=[Patch(facecolor="navy"),Patch(facecolor="crimson")]
handles2=[Patch(facecolor="darkturquoise"),Patch(facecolor="orange"),Patch(facecolor="rosybrown")]

a=plt.legend(handles,["Female","Male"],title="Sex",bbox_to_anchor=(0.7,0.55),bbox_transform=plt.gcf().transFigure,loc='upper right',prop={'size':11},title_fontsize=11)
a2=plt.legend(handles2,["Moderate","Severe","Critical"],title="Type",bbox_to_anchor=(0.7,0.75),bbox_transform=plt.gcf().transFigure,loc='upper right',prop={'size':11},title_fontsize=11)

plt.gca().add_artist(a)

plt.show()


