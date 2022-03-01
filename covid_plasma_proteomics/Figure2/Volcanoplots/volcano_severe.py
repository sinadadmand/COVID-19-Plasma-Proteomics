import pandas as pd
import matplotlib.pyplot as plt
import seaborn
import numpy as np
from scipy import stats
import statsmodels.stats as sts
import permutation_test as p
from scipy.stats import mannwhitneyu,fisher_exact
from scipy import stats
import seaborn as sns
import matplotlib


path1="/home/ali/Desktop/combined_data/combined_filtered_severe.xlsx"
path2="/home/ali/Desktop/combined_data/combined_filtered_healthy.xlsx"

# Let me read the data
data1=pd.read_excel(path1)
data2=pd.read_excel(path2)
# data3=pd.read_excel(path3)
# Let me merge the combined data
z= pd.merge(data1,data2, on="Accession", how="outer")
# z= pd.merge(z,data3, on="Accession", how="outer")

data_filter= z[["Accession","genes_2","genes_4","S1_1","S1_2","S1_3","S2_1","S2_2","S2_3","S3_1","S3_2",
                   "H1","H2","H3","H4","H5"]]


data_filter["genes_2"].fillna(data_filter["genes_4"],inplace=True)

data_filter=data_filter[["Accession","genes_2","S1_1","S1_2","S1_3","S2_1","S2_2","S2_3","S3_1","S3_2",
                   "H1","H2","H3","H4","H5"]]



a=data_filter.min(axis=1,skipna=True)
print(a.min())


# # Let me filter out 0.2 % threshold of the nan values and fill the remaining nans with small value of 0.001
data_filter= data_filter.dropna(thresh=len(data_filter.T)*(0.4))
data_filter.fillna(0.0001,inplace=True)
data_filter.reset_index(inplace=True)
data_filter=data_filter[data_filter.columns[1:]]
print(data_filter)


# Let me construct the critical data by taking mean values of each sample points between the patients!!

data_filter["Severe1"]=data_filter[["S1_1","S1_2","S1_3"]].mean(axis=1)
data_filter["Severe2"]=data_filter[["S2_1","S2_2","S2_3"]].mean(axis=1)
data_filter["Severe3"]=data_filter[["S3_1","S3_2"]].mean(axis=1)

data_filter=data_filter[["Accession","genes_2","H1","H2","H3","H4","H5","Severe1","Severe2","Severe3"]]

df_filter_only_values=data_filter[["Severe1","Severe2","Severe3"]]
df_filter2_only_values=data_filter[["H1","H2","H3","H4","H5"]]

print(df_filter_only_values)
print(df_filter2_only_values)

transposed_data1=df_filter_only_values.T
transposed_data2=df_filter2_only_values.T


my_overall=pd.DataFrame({"Accession":data_filter["Accession"].tolist(),"Protein Name":data_filter["genes_2"].tolist(),"First Rep1":data_filter["Severe1"],"First Rep2":data_filter["Severe2"],
                         "First Rep3":data_filter["Severe3"],
                            "Second Rep1":data_filter["H1"],"Second Rep2":data_filter["H2"],"Second Rep3":data_filter["H3"],"Second Rep4":data_filter["H4"],
                            "Second Rep5":data_filter["H5"]
                         })



my_overall['p_values']= pd.DataFrame(stats.ttest_ind(transposed_data1,transposed_data2,axis=0,equal_var=False)[1])


my_overall['-log10 p value']=-np.log10(my_overall['p_values'])
my_overall.loc[my_overall['p_values']<=0.05,'significance']='+'
my_overall.loc[my_overall['p_values']>0.05,'significance']='-'

ppp=((my_overall["First Rep1"]+my_overall["First Rep2"]+my_overall["First Rep3"])/3)/((my_overall["Second Rep1"]+my_overall["Second Rep2"]+my_overall["Second Rep3"]+my_overall["Second Rep4"]+my_overall["Second Rep5"])/5)

my_overall['log2(Mean(First/Second))']= np.log2(ppp)

my_overall=my_overall[my_overall["p_values"]!=0]

S= my_overall[my_overall["p_values"]<0.05]

my_overall.loc[(my_overall['significance']=='+')&(my_overall['log2(Mean(First/Second))']>0),'Sample downregulated (-) and upregulated (+) proteins']='+'
my_overall.loc[(my_overall['significance']=='+')&(my_overall['log2(Mean(First/Second))']<-0),'Sample downregulated (-) and upregulated (+) proteins']='-'
my_overall['Sample downregulated (-) and upregulated (+) proteins']=my_overall['Sample downregulated (-) and upregulated (+) proteins'].fillna('Not Valid')

df1_Int_valid=my_overall[(my_overall['Sample downregulated (-) and upregulated (+) proteins']=='+')]
df1_Int_valid2=my_overall[my_overall['Sample downregulated (-) and upregulated (+) proteins']=='-']
df1_Int_notvalid=my_overall[my_overall['Sample downregulated (-) and upregulated (+) proteins']=='Not Valid']


path_mild_up="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2C/critical_vs_healthy/upregulated_ones.xlsx"
path_mild_down="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2C/critical_vs_healthy/downregulated_ones.xlsx"
path_severe_up="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2A/moderate_vs_healthy/upregulated_ones.xlsx"
path_severe_down="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2A/moderate_vs_healthy/downregulated_ones.xlsx"

severe_u=pd.read_excel(path_severe_up)
severe_d=pd.read_excel(path_severe_down)

mild_u=pd.read_excel(path_mild_up)
mild_d=pd.read_excel(path_mild_down)
mild=pd.concat([mild_u,mild_d],axis=0,ignore_index=True)
severe=pd.concat([severe_d,severe_u],axis=0,ignore_index=True)
critical=pd.concat([df1_Int_valid,df1_Int_valid2],axis=0,ignore_index=True)

critical_mild=critical[(critical["Accession"].isin(mild["Accession"]))&(~critical["Accession"].isin(severe["Accession"]))]
critical_severe=critical[(~critical["Accession"].isin(mild["Accession"]))&(critical["Accession"].isin(severe["Accession"]))]
common=critical[(critical["Accession"].isin(mild["Accession"]))&(critical["Accession"].isin(severe["Accession"]))]
only_critical=critical[(~critical["Accession"].isin(mild["Accession"]))&(~critical["Accession"].isin(severe["Accession"]))]




df1_Int_valid.sort_values(by="log2(Mean(First/Second))",inplace=True,ignore_index=True,ascending=False)
df1_Int_valid2.sort_values(by="log2(Mean(First/Second))",inplace=True,ignore_index=True,ascending=True)

df1_Int_valid.to_excel("/home/ali/Desktop/upregulated_ones.xlsx")
df1_Int_valid2.to_excel("/home/ali/Desktop/downregulated_ones.xlsx")

count=0
my_gene_coordinates={}



highest= critical.sort_values(by="log2(Mean(First/Second))",ascending=False)
highest.to_excel("/home/ali/Desktop/neybu.xlsx")
#
lowest=critical.sort_values(by="log2(Mean(First/Second))",ascending=True)
first_10={}

for j in range(0,55):
    x_axis=highest.iloc[j][13]
    y_axis=highest.iloc[j][11]
    first_10[highest.iloc[j][1]]=[x_axis,y_axis]

print(first_10)

least7={}
for j in range(0,10):
    x_axis=lowest.iloc[j][13]
    y_axis=lowest.iloc[j][11]
    least7[lowest.iloc[j][1]]=[x_axis,y_axis]

print(first_10)





fig,ax=plt.subplots()
plt.plot('log2(Mean(First/Second))','-log10 p value',data=critical_severe,linestyle='',marker='o', markersize=5.5,alpha=1,color='blue',label="severe-moderate (n=11)")
plt.plot('log2(Mean(First/Second))','-log10 p value',data=common,linestyle='',marker='o', markersize=5.5,alpha=1,color='purple',label="common (n=26)")
plt.plot('log2(Mean(First/Second))','-log10 p value',data=only_critical,linestyle='',marker='o', markersize=5.5,alpha=1,color='orange',label="only severe (n=20)")
plt.plot('log2(Mean(First/Second))','-log10 p value',data=critical_mild,linestyle='',marker='o', markersize=5.5,alpha=1,color='green',label="severe-critical (n=8)")

marker_list=["o","v"]
string=["No previous data","Previously found"]

#
h = [plt.plot([],[], color="black", marker=marker_list[i],ls="")[0] for i in range(0,2)]


plt.plot('log2(Mean(First/Second))','-log10 p value',data=df1_Int_notvalid,linestyle='',marker='o', markersize=5.5,alpha=0.50,color='gray')

# plt.axvline(-0.3, color='black', linestyle='dashed', linewidth=1)
plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
plt.hlines(-np.log10(0.05),-10,10,colors='black',linestyles='dashed',linewidth=1)
plt.xlim(-8,11)
plt.ylabel(("-log10(P value)"))
plt.xlabel(("log2(fold change)"))


s= critical.iloc[1][0]




ax.annotate(list(first_10.keys())[0],(list(first_10.values())[0][0]-0.45,list(first_10.values())[0][1]+0.1),fontsize=12)
ax.annotate(list(first_10.keys())[1],(list(first_10.values())[1][0]+0.1,list(first_10.values())[1][1]),fontsize=12)
ax.annotate(list(first_10.keys())[2],(list(first_10.values())[2][0]+0.14,list(first_10.values())[2][1]-0.05),fontsize=12)
ax.annotate(list(first_10.keys())[3],(list(first_10.values())[3][0]+0.15,list(first_10.values())[3][1]),fontsize=12)
ax.annotate(list(first_10.keys())[4],(list(first_10.values())[4][0]+0.1,list(first_10.values())[4][1]),fontsize=12)
ax.annotate(list(first_10.keys())[5],(list(first_10.values())[5][0]+0.15,list(first_10.values())[5][1]-0.07),fontsize=12)
ax.annotate(list(first_10.keys())[6],(list(first_10.values())[6][0]+0.15,list(first_10.values())[6][1]+0.03),fontsize=12)
ax.annotate(list(first_10.keys())[7],(list(first_10.values())[7][0]+0.15,list(first_10.values())[7][1]-0.01),fontsize=12)
ax.annotate(list(first_10.keys())[8],(list(first_10.values())[8][0],list(first_10.values())[8][1]+0.09),fontsize=12)
ax.annotate(list(first_10.keys())[9],(list(first_10.values())[9][0]+0.1,list(first_10.values())[9][1]+0.03),fontsize=12)
# ax.annotate(list(first_10.keys())[10],(list(first_10.values())[10][0],list(first_10.values())[10][1]+0.26),fontsize=12)
# ax.annotate(list(first_10.keys())[11],(list(first_10.values())[11][0]+0.05,list(first_10.values())[11][1]),fontsize=12)



ax.annotate(list(first_10.keys())[12],(list(first_10.values())[12][0]-0.15,list(first_10.values())[12][1]+0.05),fontsize=12)
#
# ax.annotate(list(first_10.keys())[13],(list(first_10.values())[13][0]+0.15,list(first_10.values())[13][1]+0.05),fontsize=12)

# ax.annotate(list(first_10.keys())[14],(list(first_10.values())[14][0],list(first_10.values())[14][1]+0.05),fontsize=12)
# ax.annotate(list(first_10.keys())[15],(list(first_10.values())[15][0],list(first_10.values())[15][1]),fontsize=12)
ax.annotate(list(first_10.keys())[16],(list(first_10.values())[16][0]+0.15,list(first_10.values())[16][1]),fontsize=12)
ax.annotate(list(first_10.keys())[17],(list(first_10.values())[17][0]+0.13,list(first_10.values())[17][1]+0.05),fontsize=12)
# ax.annotate(list(first_10.keys())[20],(list(first_10.values())[20][0]+0.13,list(first_10.values())[20][1]),fontsize=12)
# ax.annotate(list(first_10.keys())[21],(list(first_10.values())[21][0]+0.13,list(first_10.values())[21][1]),fontsize=12)
ax.annotate(list(first_10.keys())[21],(list(first_10.values())[21][0]+0.13,list(first_10.values())[21][1]+0.05),fontsize=12)

ax.annotate(list(first_10.keys())[25],(list(first_10.values())[25][0]+0.13,list(first_10.values())[25][1]+0.05),fontsize=12)
ax.annotate(list(first_10.keys())[26],(list(first_10.values())[26][0]+0.13,list(first_10.values())[26][1]+0.05),fontsize=12)

ax.annotate(list(first_10.keys())[29],(list(first_10.values())[29][0]-0.15,list(first_10.values())[29][1]+0.12),fontsize=12)
ax.annotate(list(first_10.keys())[35],(list(first_10.values())[35][0]+0.18,list(first_10.values())[35][1]),fontsize=12)
ax.annotate(list(first_10.keys())[36],(list(first_10.values())[36][0]+0.18,list(first_10.values())[36][1]),fontsize=12)
ax.annotate(list(first_10.keys())[45],(list(first_10.values())[45][0]-0.48,list(first_10.values())[45][1]+0.06),fontsize=12)
ax.annotate(list(first_10.keys())[44],(list(first_10.values())[44][0]+0.1,list(first_10.values())[44][1]),fontsize=12)
ax.annotate(list(first_10.keys())[46],(list(first_10.values())[46][0]-0.6,list(first_10.values())[46][1]+0.05),fontsize=12)

ax.annotate(list(first_10.keys())[47],(list(first_10.values())[47][0]-0.5,list(first_10.values())[47][1]+0.1),fontsize=12)
ax.annotate(list(first_10.keys())[48],(list(first_10.values())[48][0]-0.12,list(first_10.values())[48][1]-0.17),fontsize=12)
ax.annotate(list(first_10.keys())[49],(list(first_10.values())[49][0]-0.1,list(first_10.values())[49][1]+0.16),fontsize=12)
ax.annotate(list(first_10.keys())[50],(list(first_10.values())[50][0]-0.85,list(first_10.values())[50][1]+0.19),fontsize=12)
ax.annotate(list(first_10.keys())[51],(list(first_10.values())[51][0]-0.3,list(first_10.values())[51][1]+0.05),fontsize=12)
ax.annotate(list(first_10.keys())[54],(list(first_10.values())[54][0]-0.35,list(first_10.values())[54][1]+0.07),fontsize=12)
ax.annotate(list(first_10.keys())[53],(list(first_10.values())[53][0]-0.65,list(first_10.values())[53][1]+0.07),fontsize=12)
ax.annotate(list(first_10.keys())[52],(list(first_10.values())[52][0]-0.65,list(first_10.values())[52][1]+0.2),fontsize=12)

plt.xlabel('log2(Severe/Healthy)', fontsize=25)
plt.ylabel('-log10(P value)', fontsize=25)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

plt.xlim(-6.7,9)

plt.subplots_adjust(right=0.6,bottom=0.15,left=0.06,top=0.955)
plt.show()

