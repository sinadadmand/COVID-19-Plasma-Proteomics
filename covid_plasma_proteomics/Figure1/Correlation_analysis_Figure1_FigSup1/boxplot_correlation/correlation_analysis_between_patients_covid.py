import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

path_healthy="/home/ali/Desktop/combined_data/combined_filtered_healthy.xlsx"
path_severe_i="/home/ali/Desktop/combined_data/combined_filtered_severe.xlsx"
path_mild_i="/home/ali/Desktop/combined_data/combined_filtered_moderate.xlsx"
path_critical_i="/home/ali/Desktop/combined_data/combined_filtered_critical.xlsx"


healthy=pd.read_excel(path_healthy)
severe=pd.read_excel(path_severe_i)
moderate= pd.read_excel(path_mild_i)
critical= pd.read_excel(path_critical_i)

# Let me drop nan values from healthy subjects and patients, and draw the correlation graphs for them!

healthy=healthy.dropna()
healthy.reset_index(inplace=True)
# healthy= healthy[healthy.columns[2:]]
healthy=healthy[["Accession","genes_2","H1","H2","H3","H4","H5"]]
print(healthy)

moderate=moderate.dropna()
moderate.reset_index(inplace=True)
moderate= moderate[moderate.columns[2:]]
print(moderate)

severe=severe.dropna()
severe.reset_index(inplace=True)
severe= severe[severe.columns[2:]]
print(severe)

critical=critical.dropna()
critical.reset_index(inplace=True)
critical= critical[critical.columns[2:]]
print(critical)



healthy=healthy[healthy["genes_2"]!="ALB"]
moderate=moderate[moderate["genes_3"]!="ALB"]
severe=severe[severe["genes_4"]!="ALB"]
critical=critical[critical["genes_1"]!="ALB"]

healthy=healthy[healthy.columns[2:]]
moderate=moderate[moderate.columns[2:]]
severe=severe[severe.columns[2:]]
critical=critical[critical.columns[2:]]


critical=critical[["C1_1","C1_2","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1",
                    "C5_2","C6_1","C6_2"]]




print(critical)
import numpy as np

fig,ax=plt.subplots(figsize=(20,20),dpi=140)

correlate1=healthy.corr(method='pearson')
correlate1=correlate1.where(np.tril(np.ones(correlate1.shape)).astype(np.bool))

a= sns.heatmap(correlate1,cmap="YlGnBu",annot=False,annot_kws={"size":20})

plt.xticks(fontsize=20,rotation=90)
plt.yticks(fontsize=20,va="center",rotation=0)
cbar= a.collections[0].colorbar
cbar.ax.tick_params(labelsize=18)

overall4=[]
for j in range(0,len(correlate1.columns)):
    mylist=correlate1[correlate1.columns[j]].to_list()
    for k in mylist:
        overall4.append(k)

healthy_score = [x for x in overall4 if str(x) != 'nan' and str(x)!=str(1.0)]



plt.show()


fig,ax=plt.subplots(figsize=(20,20),dpi=140)


correlate2=moderate.corr(method='pearson')
correlate2=correlate2.where(np.tril(np.ones(correlate2.shape)).astype(np.bool))
correlate2.to_excel("/home/ali/Desktop/moderate_check_corr.xlsx")

a2= sns.heatmap(correlate2,cmap="YlGnBu",annot=False,annot_kws={"size":20})
plt.xticks(fontsize=20,rotation=90)
plt.yticks(fontsize=20,va="center",rotation=0)
cbar= a2.collections[0].colorbar
cbar.ax.tick_params(labelsize=18)


overall3=[]
for j in range(0,len(correlate2.columns)):
    mylist=correlate2[correlate2.columns[j]].to_list()
    for k in mylist:
        overall3.append(k)

moderate_score = [x for x in overall3 if str(x) != 'nan' and str(x)!=str(1.0)]

print(moderate_score)


plt.show()

#
#
fig,ax=plt.subplots(figsize=(20,20),dpi=140)
#
correlate3=severe.corr(method='pearson')
correlate3=correlate3.where(np.tril(np.ones(correlate3.shape)).astype(np.bool))

a3= sns.heatmap(correlate3,cmap="YlGnBu",annot=False,annot_kws={"size":20})
plt.xticks(fontsize=20,rotation=90)
plt.yticks(fontsize=20,va="center",rotation=0)
cbar= a3.collections[0].colorbar
cbar.ax.tick_params(labelsize=18)

overall2=[]
for j in range(0,len(correlate3.columns)):
    mylist=correlate3[correlate3.columns[j]].to_list()
    for k in mylist:
        overall2.append(k)

severe_score = [x for x in overall2 if str(x) != 'nan' and str(x)!=str(1.0)]

print(severe_score)

plt.show()
#

fig,ax=plt.subplots(figsize=(20,20),dpi=140)

correlate4=critical.corr(method='pearson')
correlate4=correlate4.where(np.tril(np.ones(correlate4.shape)).astype(np.bool))
print(correlate4)

# I am goinf to form a box plot for correlation analysis!

a4= sns.heatmap(correlate4,cmap="YlGnBu",annot=False,annot_kws={"size":17})
plt.xticks(fontsize=20,rotation=90)
plt.yticks(fontsize=20,va="center",rotation=0)
cbar= a4.collections[0].colorbar
cbar.ax.tick_params(labelsize=18)

overall=[]
for j in range(0,len(correlate4.columns)):
    mylist=correlate4[correlate4.columns[j]].to_list()
    for k in mylist:
        overall.append(k)

critical_score = [x for x in overall if str(x) != 'nan' and str(x)!=str(1.0)]

plt.show()


s4=["Critical"]*len(critical_score)
s1=["Healthy"]*len(healthy_score)
s2=["Moderate"]*len(moderate_score)
s3=["Severe"]*len(severe_score)

A=healthy_score+moderate_score+severe_score+critical_score
S=s1+s2+s3+s4

box_plot_dataframe=pd.DataFrame({"Pearson Correlation Coefficient (R)":A,"sample":S})


print(box_plot_dataframe)

# box_plot_dataframe.to_excel("/home/a/Desktop/res.xlsx")

sns.set(style="whitegrid")

print(box_plot_dataframe)


from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statannot import add_stat_annotation


tukey=pairwise_tukeyhsd(endog=box_plot_dataframe["Pearson Correlation Coefficient (R)"],groups=box_plot_dataframe["sample"],alpha=0.5)

print(tukey)
mysample=["Healthy","Moderate","Severe","Critical"]


palette={"Healthy":"gray","Moderate":"deepskyblue","Severe":"green","Critical":"brown"}

flierprops=dict(marker='o',markersize=3)
ax=sns.boxplot(x="sample",y="Pearson Correlation Coefficient (R)",data=box_plot_dataframe, flierprops=flierprops,palette=palette)
# ax=sns.swarmplot(x="sample",y="Pearson Correlation Coefficient (R)",data=box_plot_dataframe,color=".25")
ax.set(xlabel=None)
add_stat_annotation(ax,data=box_plot_dataframe ,x="sample", y="Pearson Correlation Coefficient (R)", order=mysample,
                    box_pairs=[(mysample[0], mysample[1]),(mysample[0], mysample[2]),(mysample[0], mysample[3]),

                              ],
                    test='Kruskal', text_format='star', loc='outside', verbose=0,fontsize=15)


plt.ylim(0.8,1)
plt.ylabel("Pearson Correlation Coefficient (R)",fontsize=20)
# plt.xlabel(fontsize=14)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.show()
#
