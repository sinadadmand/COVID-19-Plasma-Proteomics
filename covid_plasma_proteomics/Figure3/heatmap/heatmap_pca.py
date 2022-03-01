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

myPSM_filter_df= z.dropna(thresh=0)

myPSM_filter_df.fillna(0.01,inplace=True)
myPSM_filter_df.reset_index(inplace=True)
myPSM_filter_df=myPSM_filter_df[myPSM_filter_df.columns[1:]]
print(myPSM_filter_df)


searcher_p="/home/ali/Documents/projects/COVID_projects/results/last_version-20211018T081343Z-001/last_version/critical_healthy/critical_healthy.xlsx"
searcher_m="/home/ali/Documents/projects/COVID_projects/results/last_version-20211018T081343Z-001/last_version/moderate_healthy/moderate_healthy.xlsx"
searcher_s="/home/ali/Documents/projects/COVID_projects/results/last_version-20211018T081343Z-001/last_version/severe_healthy/severe_healthy.xlsx"


path_="/home/ali/Desktop/heatmap_after_80percent_filtration.xlsx"
search=pd.read_excel(path_)

myPSM_filter_df=myPSM_filter_df[(myPSM_filter_df["Accession"].isin(search["Accession"]))]



myPSM_filter_df.reset_index(inplace=True)
myPSM_filter_df=myPSM_filter_df[myPSM_filter_df.columns[1:]]
myPSM_filter_df.to_excel("/home/ali/Desktop/pca_data.xlsx")


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

myPSM_filter_df["C1_1"]=myPSM_filter_df["C1_1"]/avg
myPSM_filter_df["C1_2"]=myPSM_filter_df["C1_2"]/avg
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

myPSM_filter_df["C1RR"]=myPSM_filter_df["C1RR"]/avg
myPSM_filter_df["C2RR"]=myPSM_filter_df["C2RR"]/avg
myPSM_filter_df["C3RR"]=myPSM_filter_df["C3RR"]/avg
myPSM_filter_df["C4RR"]=myPSM_filter_df["C4RR"]/avg
myPSM_filter_df["C5RR"]=myPSM_filter_df["C5RR"]/avg

myPSM_filter_df["C1R"]=myPSM_filter_df["C1R"]/avg
myPSM_filter_df["C2R"]=myPSM_filter_df["C2R"]/avg
myPSM_filter_df["C3R"]=myPSM_filter_df["C3R"]/avg
myPSM_filter_df["C4R"]=myPSM_filter_df["C4R"]/avg
myPSM_filter_df["C5R"]=myPSM_filter_df["C5R"]/avg


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




myPSM_filter_df=myPSM_filter_df[["Accession","genes_1","E1","E2","E3","E4","E5","E6","E7","E8","M1","M2","M3","M4","M5","M6","M7","M8",
                                 "M9","M10","M11","M12","M13","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10"]]

genes=myPSM_filter_df["genes_1"].tolist()
zz= myPSM_filter_df[myPSM_filter_df.columns[2:]]


zz=np.log2(zz)





cluster2=[]

colu=zz.columns

pallete=sns.color_palette()

color_dict={}
color_dict2={}
color_dict3={}

for col in zz.columns:
    # if  col=="H1" or col=="H2" or col=="H3":
    #     color_dict2[col]="lightgreen"
    if col=="E1" or col=="E2" or col=="E3"or col=="E4"or col=="E5"or col=="E6":
        color_dict2[col]="#DAC4F7"
    elif col=="M1" or col=="M2" or col=="M3"or col=="M4" or col=="M5"or col=="M6" or \
            col=="M7"or col=="M8"or col=="M9"or col=="M10"or col=="M11"or col=="M12"or col=="M13":
        color_dict2[col]="#F4989C"
    else:
        color_dict2[col]="#D6F6DD"

color_col2=pd.Series(color_dict2)
for col in zz.columns:
    if col=="E1"or col=="E2"or col=="M1"or col=="M2" or col=="M3"or col=="M4" or col=="R1"or col=="R2":
        color_dict[col]="darkturquoise"
    elif col=="E3"or col=="E4"or col=="M5"or col=="M6" or col=="M7"or col=="R3"or col=="R4"or col=="R5":
        color_dict[col] = "orange"

    else:
        color_dict[col] = "brown"

color_col=pd.Series(color_dict)

k= sns.clustermap(zz,method="complete",metric="canberra",cmap="coolwarm",row_cluster=False,col_cluster=True,col_colors=[color_col,color_col2],yticklabels=genes,vmin=-6,vmax=6,xticklabels=False,colors_ratio=0.015,cbar_kws={"label":"log2(FC)"},tree_kws={'colors':["slateblue"]*19+["green"]*10+["black"]*1,})


for a in k.ax_col_dendrogram.collections:
    a.set_linewidth(1.5)



pp=["darkturquoise","orange","brown"]
counter=0
for i in ["Moderate","Severe","Critical"]:
    k.ax_col_dendrogram.bar(0,0,color=pp[counter],label=i)
    counter+=1

s=k.ax_col_dendrogram.legend(title='Type', loc="upper right", ncol=1, bbox_to_anchor=(0.37,0.75), bbox_transform=plt.gcf().transFigure,prop={'size':11},title_fontsize=11)



plt.setp(k.ax_heatmap.yaxis.get_majorticklabels(),rotation=0,fontsize=5,fontstyle="normal")
plt.subplots_adjust(top=0.95,right=0.5,wspace=0,bottom=0.1)

k.cax.set_position([0.6,.22,0.02,0.15])

my_color1=[pallete[0],(52, 205, 235)]

handles=[Patch(facecolor="#DAC4F7"),Patch(facecolor="#F4989C"),Patch(facecolor="#D6F6DD")]
handles2=[Patch(facecolor="darkturquoise"),Patch(facecolor="orange"),Patch(facecolor="brown")]

a=plt.legend(handles,["Early Infection","Maximum Infection","Post-Infection"],title="Time",bbox_to_anchor=(0.7,0.57),bbox_transform=plt.gcf().transFigure,loc='upper right',prop={'size':11},title_fontsize=11)
a2=plt.legend(handles2,["Moderate","Severe","Critical"],title="Type",bbox_to_anchor=(0.67,0.75),bbox_transform=plt.gcf().transFigure,loc='upper right',prop={'size':11},title_fontsize=11)

plt.gca().add_artist(a)

plt.show()
