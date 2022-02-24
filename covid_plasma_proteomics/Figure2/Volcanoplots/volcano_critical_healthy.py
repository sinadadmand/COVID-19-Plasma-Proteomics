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
path1="/home/ali/Desktop/combined_data/combined_filtered_critical.xlsx"
path2="/home/ali/Desktop/combined_data/combined_filtered_healthy.xlsx"

# Let me read the data
data1=pd.read_excel(path1)
data2=pd.read_excel(path2)

# Let me merge the combined data
z= pd.merge(data1,data2, on="Accession", how="outer")
data_filter= z[["Accession","genes_1","genes_2","H1","H2","H3","H4","H5","C1_1","C1_2","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1","C5_2","C6_1",
                   "C6_2"]]

data_filter["genes_1"].fillna(data_filter["genes_2"],inplace=True)
data_filter["genes_2"].fillna(data_filter["genes_1"],inplace=True)

data_filter=data_filter[["Accession","genes_1","H1","H2","H3","H4","H5","C1_1","C1_2","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1","C5_2","C6_1",
                   "C6_2"]]

# Let me filter out 0.2 % threshold of the nan values and fill the remaining nans with small value of 0.001

# a=data_filter.min(axis=1,skipna=True)
# print(a.min())

# data_filter.to_excel("/home/a/Desktop/res.xlsx")

data_filter= data_filter.dropna(thresh=len(data_filter.T)*(0.40))
data_filter.fillna(0.0001,inplace=True)
data_filter.reset_index(inplace=True)
data_filter=data_filter[data_filter.columns[1:]]
print(data_filter)

#
# data_filter.to_excel("/home/a/Desktop/before.xlsx")

# Let me construct the critical data by taking mean values of each sample points between the patients!!

# data_filter.to_excel("/home/a/Desktop/then.xlsx")
data_filter["Critical1"]=(data_filter["C1_1"]+data_filter["C1_2"])/2
data_filter["Critical2"]=(data_filter["C2_1"]+data_filter["C2_3"]+data_filter["C2_2"])/3
data_filter["Critical3"]=(data_filter["C3_1"]+data_filter["C3_2"])/2
data_filter["Critical4"]=(data_filter["C4_1"]+data_filter["C4_2"])/2
data_filter["Critical5"]=(data_filter["C5_1"]+data_filter["C5_2"])/2
data_filter["Critical6"]=(data_filter["C6_1"]+data_filter["C6_2"])/2
# data_filter["Critical8"]=(data_filter["C8_1"]+data_filter["C8_2"])/2

# data_filter=data1[data1["Accession"].isin(data2["Accession"])]

# z=data1[~data1["Accession"].isin(data2["Accession"])]
#
# z.to_excel("/home/a/Desktop/notfound.xlsx")
#
# data2_filter=data2[data2["Accession"].isin(data1["Accession"])]

# data_filter=data_filter[data_filter.columns[1:]]
# data2_filter=data2_filter[data2_filter.columns[1:]]
#
# print(data2_filter)
# print(data_filter)
 # #
# data_filter=data_filter.sort_values(by="Accession",ignore_index=True)
# data2_filter=data2_filter.sort_values(by="Accession",ignore_index=True)
#
#
# # data2_filter=data2_filter[["Accession","genes","first","second","third"]]
#
# # print(data2_filter)
#
# print("number = {}".format(len(data_filter)))
#
# data2_filter["fourth"]=data2_filter["first"]
# data2_filter["fifth"]=data2_filter["second"]
# data2_filter["sixth"]=data2_filter["third"]
# data2_filter["seventh"]=(data2_filter["first"])
# data2_filter["eighth"]=data2_filter["second"]
#
# df_filter_only_values=data_filter[data_filter.columns[2:]]
# df_filter2_only_values=data2_filter[data2_filter.columns[2:]]
# #

df_filter_only_values=data_filter[["Critical1","Critical2","Critical3","Critical4","Critical5","Critical6"]]
# df_filter2_only_values=data_filter[["H1","H9","H3","H5","H6","H7","H8"]]
df_filter2_only_values=data_filter[["H1","H2","H3","H4","H5"]]
# data_filter.to_excel("/home/a/Desktop/as.xlsx")

print(df_filter_only_values)
print(df_filter2_only_values)

transposed_data1=df_filter_only_values.T
transposed_data2=df_filter2_only_values.T

#
my_overall=pd.DataFrame({"Accession":data_filter["Accession"].tolist(),"Protein Name":data_filter["genes_1"].tolist(),"First Rep1":data_filter["Critical1"],"First Rep2":data_filter["Critical2"],
                         "First Rep3":data_filter["Critical3"],"First Rep4":data_filter["Critical4"],
                        "First Rep5":data_filter["Critical5"],"First Rep6":data_filter["Critical6"],"Second Rep1":data_filter["H1"],
                         "Second Rep2":data_filter["H2"],"Second Rep3":data_filter["H3"],"Second Rep4":data_filter["H4"]
                            ,"Second Rep5":data_filter["H5"]
                         })

# my_overall.to_excel("/home/a/Desktop/critical_mutual.xlsx")
print(transposed_data1)
print(transposed_data2)
# my_overall.to_excel("/home/a/Desktop/my_excell.xlsx")
#
my_overall['p_values']= pd.DataFrame(stats.ttest_ind(transposed_data1,transposed_data2,axis=0,equal_var=False)[1])


# my_overall.to_excel("/home/a/Desktop/res.xlsx")

my_overall['-log10 p value']=-np.log10(my_overall['p_values'])

# my_overall.to_excel("/home/a/Desktop/myoverall.xlsx")
my_overall.loc[my_overall['p_values']<=0.05,'significance']='+'
my_overall.loc[my_overall['p_values']>0.05,'significance']='-'

ppp=((my_overall["First Rep1"]+my_overall["First Rep2"]+my_overall["First Rep3"]+my_overall["First Rep4"]+my_overall["First Rep5"]+my_overall["First Rep6"])/6)/((my_overall["Second Rep1"]+my_overall["Second Rep2"]+my_overall["Second Rep3"]+my_overall["Second Rep4"]+my_overall["Second Rep5"])/5)

my_overall['log2(Mean(First/Second))']= np.log2(ppp)
print(len(my_overall))
print(len(my_overall[my_overall["log2(Mean(First/Second))"]>0]))

my_overall=my_overall[my_overall["p_values"]!=0]

S= my_overall[my_overall["p_values"]<0.05]
print(len(S))

my_overall.to_excel("/home/ali/Desktop/result.xlsx")

my_overall.loc[(my_overall['significance']=='+')&(my_overall['log2(Mean(First/Second))']>0),'Sample downregulated (-) and upregulated (+) proteins']='+'
my_overall.loc[(my_overall['significance']=='+')&(my_overall['log2(Mean(First/Second))']<-0),'Sample downregulated (-) and upregulated (+) proteins']='-'
my_overall.loc[(my_overall['significance']=='+')&(my_overall['log2(Mean(First/Second))']<0)&(my_overall['log2(Mean(First/Second))']>-0),'Sample downregulated (-) and upregulated (+) proteins']='0.5 cutt-off_extracted'

my_overall['Sample downregulated (-) and upregulated (+) proteins']=my_overall['Sample downregulated (-) and upregulated (+) proteins'].fillna('Not Valid')
# my_overall.to_excel("/home/a/Desktop/resso.xlsx")

df1_Int_valid=my_overall[(my_overall['Sample downregulated (-) and upregulated (+) proteins']=='+')]
df1_Int_valid2=my_overall[my_overall['Sample downregulated (-) and upregulated (+) proteins']=='-']
df1_Int_notvalid=my_overall[my_overall['Sample downregulated (-) and upregulated (+) proteins']=='Not Valid']



path_mild_up="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2A/moderate_vs_healthy/upregulated_ones.xlsx"
path_mild_down="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2A/moderate_vs_healthy/downregulated_ones.xlsx"
path_severe_up="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2B/severe_vs_healthy/upregulated_ones.xlsx"
path_severe_down="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2B/severe_vs_healthy/downregulated_ones.xlsx"

severe_u=pd.read_excel(path_severe_up)
severe_d=pd.read_excel(path_severe_down)

mild_u=pd.read_excel(path_mild_up)
mild_d=pd.read_excel(path_mild_down)
mild=pd.concat([mild_u,mild_d],axis=0,ignore_index=True)
severe=pd.concat([severe_d,severe_u],axis=0,ignore_index=True)
print(len(df1_Int_valid))
print(len(df1_Int_valid2))
critical=pd.concat([df1_Int_valid,df1_Int_valid2],axis=0,ignore_index=True)

critical_mild=critical[(critical["Accession"].isin(mild["Accession"]))&(~critical["Accession"].isin(severe["Accession"]))]
critical_severe=critical[(~critical["Accession"].isin(mild["Accession"]))&(critical["Accession"].isin(severe["Accession"]))]
common=critical[(critical["Accession"].isin(mild["Accession"]))&(critical["Accession"].isin(severe["Accession"]))]
only_critical=critical[(~critical["Accession"].isin(mild["Accession"]))&(~critical["Accession"].isin(severe["Accession"]))]

print(len(critical_mild))
print(len(critical_severe))
print(len(common))
print(len(only_critical))
# found_in_mild=
#
#
#


print(df1_Int_valid["log2(Mean(First/Second))"])
print(df1_Int_valid2["log2(Mean(First/Second))"])

print(len(df1_Int_valid))
print(len(df1_Int_valid2))
print(len(df1_Int_notvalid)+len(df1_Int_valid2)+len(df1_Int_valid))


df1_Int_valid.sort_values(by="log2(Mean(First/Second))",inplace=True,ignore_index=True,ascending=False)
df1_Int_valid2.sort_values(by="log2(Mean(First/Second))",inplace=True,ignore_index=True,ascending=True)
df1_Int_valid.to_excel("/home/ali/Desktop/upregulated_ones.xlsx")
df1_Int_valid2.to_excel("/home/ali/Desktop/downregulated_ones.xlsx")

# my_proteome_path="/home/ali/Desktop/masaustu/COVID19_proteomics_database.txt"
#
# proteome=pd.read_csv(my_proteome_path,delimiter="\n")
#
# proteome=proteome[proteome.columns[0]].tolist()
#
#
# proteome_names=[]
# for each in proteome:
#     if each.startswith("#"):
#         continue
#     else:
#         proteome_names.append(each)
# print(proteome_names)



# z=["PI16","LCP1","FETUB","AZGP1","CRP","CFI","ORM1","ORM2","SERPINA3","S100A9","CETP"]
count=0
my_gene_coordinates={}

print(len(critical))

# critical.loc[critical["Protein Name"].isin(proteome_names),"Literature"]="+"
# critical.loc[~critical["Protein Name"].isin(proteome_names),"Literature"]="-"
# print(proteome)
# critical.sort_values(by="Literature",inplace=True)
# critical.reset_index(inplace=True)
# critical=critical[critical.columns[1:]]
# critical.to_excel("/home/a/Desktop/new_critical.xlsx")

# for j in critical["Protein Name"].tolist():
#     for k in proteome_names:
#         if j==k:
#             xaxis= critical.iloc[count][16]
#             yaxis= critical.iloc[count][14]
#             my_gene_coordinates[j]=[xaxis,yaxis]
#     count+=1
#
# print(critical.columns)

print(len(critical))
# print(my_gene_coordinates)
print(len(my_gene_coordinates))


highest= critical.sort_values(by="log2(Mean(First/Second))",ascending=False)
highest.to_excel("/home/ali/Desktop/neybu.xlsx")




lowest=critical.sort_values(by="log2(Mean(First/Second))",ascending=True)
first_10={}

for j in range(0,61):
    x_axis=highest.iloc[j][16]
    y_axis=highest.iloc[j][14]
    first_10[highest.iloc[j][1]]=[x_axis,y_axis]

print(first_10)
#
# least7={}
# for j in range(0,25):
#     x_axis=lowest.iloc[j][15]
#     y_axis=lowest.iloc[j][13]
#     least7[lowest.iloc[j][1]]=[x_axis,y_axis]
#

# print(first_10)
#

#
# print(len(critical_severe))

# critical_severe_v= critical_severe[critical_severe["Protein Name"].isin(proteome_names)]
# critical_severe_n= critical_severe[~critical_severe["Protein Name"].isin(proteome_names)]
#
# common_v= common[common["Protein Name"].isin(proteome_names)]
# common_n= common[~common["Protein Name"].isin(proteome_names)]
#
# only_critical_v=only_critical[only_critical["Protein Name"].isin(proteome_names)]
# only_critical_n= only_critical[~only_critical["Protein Name"].isin(proteome_names)]


# critical_mild_v=critical_mild[critical_mild["Protein Name"].isin(proteome_names)]
# critical_mild_n= critical_mild[~critical_mild["Protein Name"].isin(proteome_names)]
# print(len(critical_mild_n)+len(critical_mild_v)+len(only_critical_n)+len(only_critical_v)+len(common_n)+len(common_v)+len(critical_severe_v)+len(critical_severe_n))
#
fig,ax=plt.subplots()
# plt.plot('log2(Mean(First/Second))','-log10 p value',data=critical_severe_n,linestyle='',marker='o', markersize=5.5,alpha=5,color='orange',label="critical-severe")
#
# plt.plot('log2(Mean(First/Second))','-log10 p value',data=common_n,linestyle='',marker='o', markersize=5.5,alpha=5,color='purple',label="common")
# plt.plot('log2(Mean(First/Second))','-log10 p value',data=only_critical_n,linestyle='',marker='o', markersize=5.5,alpha=5,color='brown',label="only critical")
# plt.plot('log2(Mean(First/Second))','-log10 p value',data=critical_mild_n,linestyle='',marker='o', markersize=5.5,alpha=5,color='darkturquoise',label="critical-moderate")
#
plt.plot('log2(Mean(First/Second))','-log10 p value',data=critical_severe,linestyle='',marker='o', markersize=5.5,alpha=1,color='green',label="critical-severe (n=8)")
plt.plot('log2(Mean(First/Second))','-log10 p value',data=common,linestyle='',marker='o', markersize=5.5,alpha=1,color='purple',label="common (n=26)")
plt.plot('log2(Mean(First/Second))','-log10 p value',data=only_critical,linestyle='',marker='o', markersize=5.5,alpha=1,color='brown',label="only critical (n=34)")
plt.plot('log2(Mean(First/Second))','-log10 p value',data=critical_mild,linestyle='',marker='o', markersize=5.5,alpha=1,color='pink',label="critical-moderate (n=28)")

# leg = plt.legend(loc=(1.05,0.35),title="Type",prop={'size':15},title_fontsize=15)
# ax.add_artist(leg)
#
marker_list=["o","v"]
string=["No previous data","Previously found"]
#
#
h = [plt.plot([],[], color="black", marker=marker_list[i],ls="")[0] for i in range(0,2)]
#
# print(h)
# leg2= plt.legend(handles=h, labels=string,loc=(1.03,0.3),title="Change in Expression level",prop={'size':11},title_fontsize=11)
# ax.add_artist(leg2)

# plt.legend()



plt.plot('log2(Mean(First/Second))','-log10 p value',data=df1_Int_notvalid,linestyle='',marker='o', markersize=5.5,alpha=0.50,color='gray')
# plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
# plt.axvline(0.3, color='black', linestyle='dashed', linewidth=1)
plt.hlines(-np.log10(0.05),-10,10,colors='black',linestyles='dashed',linewidth=1)
# plt.xlim(-8,11)
plt.ylabel(("-log10(P value)"))
plt.xlabel(("log2(fold change)"))
#

# print(critical[critical.columns[:]])
s= critical.iloc[1][0]


ax.annotate(list(first_10.keys())[0],(list(first_10.values())[0][0]-0.15,list(first_10.values())[0][1]+0.12),fontsize=12)
ax.annotate(list(first_10.keys())[1],(list(first_10.values())[1][0]+0.15,list(first_10.values())[1][1]+0.12),fontsize=12)
ax.annotate(list(first_10.keys())[2],(list(first_10.values())[2][0]+0.15,list(first_10.values())[2][1]+0.05),fontsize=12)
ax.annotate(list(first_10.keys())[3],(list(first_10.values())[3][0]-0.15,list(first_10.values())[3][1]+0.12),fontsize=12)
ax.annotate(list(first_10.keys())[4],(list(first_10.values())[4][0]+0.15,list(first_10.values())[4][1]+0.12),fontsize=12)
ax.annotate(list(first_10.keys())[5],(list(first_10.values())[5][0]+0.15,list(first_10.values())[5][1]),fontsize=12)
ax.annotate(list(first_10.keys())[6],(list(first_10.values())[6][0]-0.15,list(first_10.values())[6][1]+0.1),fontsize=12)
ax.annotate(list(first_10.keys())[7],(list(first_10.values())[7][0]+0.1,list(first_10.values())[7][1]+0.1),fontsize=12)
ax.annotate(list(first_10.keys())[8],(list(first_10.values())[8][0]+0.15,list(first_10.values())[8][1]-0.15),fontsize=12)
ax.annotate(list(first_10.keys())[9],(list(first_10.values())[9][0],list(first_10.values())[9][1]-0.32),fontsize=12)
ax.annotate(list(first_10.keys())[10],(list(first_10.values())[10][0]+0.1,list(first_10.values())[10][1]+0.09),fontsize=12)
ax.annotate(list(first_10.keys())[11],(list(first_10.values())[11][0]+0.1,list(first_10.values())[11][1]+0.09),fontsize=12)
ax.annotate(list(first_10.keys())[12],(list(first_10.values())[12][0]+0.1,list(first_10.values())[12][1]+0.05),fontsize=12)
ax.annotate(list(first_10.keys())[13],(list(first_10.values())[13][0]+0.13,list(first_10.values())[13][1]-0.05),fontsize=12)
ax.annotate(list(first_10.keys())[14],(list(first_10.values())[14][0]+0.1,list(first_10.values())[14][1]+0.09),fontsize=12)
ax.annotate(list(first_10.keys())[15],(list(first_10.values())[15][0]-0.3,list(first_10.values())[15][1]+0.09),fontsize=12)

# ax.annotate(list(first_10.keys())[16],(list(first_10.values())[16][0]-0.35,list(first_10.values())[16][1]+0.15),fontsize=12)
ax.annotate(list(first_10.keys())[17],(list(first_10.values())[17][0]+0.15,list(first_10.values())[17][1]-0.05),fontsize=12)

ax.annotate(list(first_10.keys())[18],(list(first_10.values())[18][0],list(first_10.values())[18][1]+0.34),fontsize=12)
ax.annotate(list(first_10.keys())[19],(list(first_10.values())[19][0]+0.1,list(first_10.values())[19][1]-0.1),fontsize=12)
# ax.annotate(list(first_10.keys())[20],(list(first_10.values())[20][0]+0.1,list(first_10.values())[20][1]-0.1),fontsize=12)
ax.annotate(list(first_10.keys())[21],(list(first_10.values())[21][0]+0.15,list(first_10.values())[21][1]),fontsize=12)
# ax.annotate(list(first_10.keys())[22],(list(first_10.values())[22][0]-0.15,list(first_10.values())[22][1]+0.21),fontsize=12)
# ax.annotate(list(first_10.keys())[23],(list(first_10.values())[23][0]+0.1,list(first_10.values())[23][1]),fontsize=12)
ax.annotate(list(first_10.keys())[24],(list(first_10.values())[24][0]+0.15,list(first_10.values())[24][1]),fontsize=12)
ax.annotate(list(first_10.keys())[25],(list(first_10.values())[25][0]+0.1,list(first_10.values())[25][1]),fontsize=12)
# ax.annotate(list(first_10.keys())[26],(list(first_10.values())[26][0]+0.1,list(first_10.values())[26][1]-0.03),fontsize=12)
# ax.annotate(list(first_10.keys())[27],(list(first_10.values())[27][0]-0.1,list(first_10.values())[27][1]+0.1),fontsize=12)
# # ax.annotate(list(first_10.keys())[28],(list(first_10.values())[28][0]+0.15,list(first_10.values())[28][1]),fontsize=12)
ax.annotate(list(first_10.keys())[29],(list(first_10.values())[29][0]+0.15,list(first_10.values())[29][1]),fontsize=12)
ax.annotate(list(first_10.keys())[30],(list(first_10.values())[30][0]+0.15,list(first_10.values())[30][1]),fontsize=12)
ax.annotate(list(first_10.keys())[31],(list(first_10.values())[31][0]+0.15,list(first_10.values())[31][1]-0.07),fontsize=12)
# ax.annotate(list(first_10.keys())[32],(list(first_10.values())[32][0]+0.1,list(first_10.values())[32][1]+0.04),fontsize=12)
ax.annotate(list(first_10.keys())[33],(list(first_10.values())[33][0]+0.12,list(first_10.values())[33][1]-0.05),fontsize=12)
# ax.annotate(list(first_10.keys())[34],(list(first_10.values())[34][0]-0.19,list(first_10.values())[34][1]+0.07),fontsize=12)
# ax.annotate(list(first_10.keys())[35],(list(first_10.values())[35][0]+0.1,list(first_10.values())[35][1]+0.04),fontsize=12)
# ax.annotate(list(first_10.keys())[36],(list(first_10.values())[36][0]+0.1,list(first_10.values())[36][1]+0.05),fontsize=12)
# # ax.annotate(list(first_10.keys())[37],(list(first_10.values())[37][0]+0.05,list(first_10.values())[37][1]-0.15),fontsize=12)
# # ax.annotate(list(first_10.keys())[38],(list(first_10.values())[38][0]+0.25,list(first_10.values())[38][1]),fontsize=12)
# # ax.annotate(list(first_10.keys())[39],(list(first_10.values())[39][0]+0.15,list(first_10.values())[39][1]),fontsize=12)
# # ax.annotate(list(first_10.keys())[40],(list(first_10.values())[40][0]+0.15,list(first_10.values())[40][1]),fontsize=12)
# ax.annotate(list(first_10.keys())[41],(list(first_10.values())[41][0]+0.1,list(first_10.values())[41][1]+0.2),fontsize=12)
# ax.annotate(list(first_10.keys())[42],(list(first_10.values())[42][0]+0.1,list(first_10.values())[42][1]),fontsize=12)
# # ax.annotate(list(first_10.keys())[43],(list(first_10.values())[43][0],list(first_10.values())[43][1]+0.3),fontsize=12)
ax.annotate(list(first_10.keys())[44],(list(first_10.values())[44][0]+0.1,list(first_10.values())[44][1]),fontsize=12)
# # ax.annotate(list(first_10.keys())[45],(list(first_10.values())[45][0]-0.28,list(first_10.values())[45][1]+0.01),fontsize=12)
ax.annotate(list(first_10.keys())[46],(list(first_10.values())[46][0]+0.15,list(first_10.values())[46][1]-0.01),fontsize=12)
ax.annotate(list(first_10.keys())[50],(list(first_10.values())[50][0]-0.34,list(first_10.values())[50][1]+0.08),fontsize=12)


ax.annotate(list(first_10.keys())[51],(list(first_10.values())[51][0]-0.02,list(first_10.values())[51][1]+0.08),fontsize=12)
ax.annotate(list(first_10.keys())[52],(list(first_10.values())[52][0],list(first_10.values())[52][1]+0.08),fontsize=12)
ax.annotate(list(first_10.keys())[53],(list(first_10.values())[53][0]-0.15,list(first_10.values())[53][1]+0.08),fontsize=12)

ax.annotate(list(first_10.keys())[54],(list(first_10.values())[54][0]-0.35,list(first_10.values())[54][1]+0.07),fontsize=12)
ax.annotate(list(first_10.keys())[55],(list(first_10.values())[55][0]-0.55,list(first_10.values())[55][1]+0.28),fontsize=12)
ax.annotate(list(first_10.keys())[56],(list(first_10.values())[56][0]-0.5,list(first_10.values())[56][1]+0.12),fontsize=12)
ax.annotate(list(first_10.keys())[57],(list(first_10.values())[57][0]-0.7,list(first_10.values())[57][1]+0.21),fontsize=12)

#
# ax.annotate(list(first_10.keys())[58],(list(first_10.values())[58][0]+0.15,list(first_10.values())[58][1]),fontsize=12)
ax.annotate(list(first_10.keys())[59],(list(first_10.values())[59][0]+0.1,list(first_10.values())[59][1]),fontsize=12)
ax.annotate(list(first_10.keys())[60],(list(first_10.values())[60][0]-0.6,list(first_10.values())[60][1]+0.1),fontsize=12)
# # ax.annotate(list(first_10.keys())[61],(list(first_10.values())[61][0],list(first_10.values())[61][1]-0.05),fontsize=8)
# #
# ax.annotate(list(first_10.keys())[61],(list(first_10.values())[61][0]+0.1,list(first_10.values())[61][1]),fontsize=12)
# ax.annotate(list(first_10.keys())[63],(list(first_10.values())[63][0]+0.15,list(first_10.values())[63][1]),fontsize=12)
# # ax.annotate(list(first_10.keys())[64],(list(first_10.values())[64][0]-0.18,list(first_10.values())[64][1]+0.05),fontsize=12)
# ax.annotate(list(first_10.keys())[65],(list(first_10.values())[65][0]+0.1,list(first_10.values())[65][1]+0.07),fontsize=12)
# # ax.annotate(list(first_10.keys())[66],(list(first_10.values())[66][0],list(first_10.values())[66][1]+0.1),fontsize=12)
# ax.annotate(list(first_10.keys())[67],(list(first_10.values())[67][0]+0.01,list(first_10.values())[67][1]),fontsize=12)
# # ax.annotate(list(first_10.keys())[68],(list(first_10.values())[68][0]+0.05,list(first_10.values())[68][1]),fontsize=12)
#
#
# ax.annotate(list(first_10.keys())[71],(list(first_10.values())[71][0]+0.15,list(first_10.values())[71][1]+0.03),fontsize=12)
# ax.annotate(list(first_10.keys())[72],(list(first_10.values())[72][0]+0.1,list(first_10.values())[72][1]-0.01),fontsize=12)
# ax.annotate(list(first_10.keys())[81],(list(first_10.values())[81][0]-0.1,list(first_10.values())[81][1]+0.27),fontsize=12)



#
# ax.annotate(list(least7.keys())[0],(list(least7.values())[0][0]-0.55,list(least7.values())[0][1]+0.07),fontsize=12)
# ax.annotate(list(least7.keys())[1],(list(least7.values())[1][0]-0.35,list(least7.values())[1][1]+0.25),fontsize=12)
# ax.annotate(list(least7.keys())[2],(list(least7.values())[2][0]-0.55,list(least7.values())[2][1]+0.1),fontsize=12)
# ax.annotate(list(least7.keys())[3],(list(least7.values())[3][0]-0.25,list(least7.values())[3][1]+0.07),fontsize=12)
# ax.annotate(list(least7.keys())[4],(list(least7.values())[4][0]-1.1,list(least7.values())[4][1]+0.1),fontsize=12)
# ax.annotate(list(least7.keys())[5],(list(least7.values())[5][0]-0.4,list(least7.values())[5][1]+0.19),fontsize=12)
# ax.annotate(list(least7.keys())[6],(list(least7.values())[6][0],list(least7.values())[6][1]+0.09),fontsize=12)
# ax.annotate(list(least7.keys())[7],(list(least7.values())[7][0]+0.1,list(least7.values())[7][1]+0.04),fontsize=12)
# ax.annotate(list(least7.keys())[8],(list(least7.values())[8][0]-0.28,list(least7.values())[8][1]+0.08),fontsize=12)






plt.xlabel('log2(Critical/Healthy)', fontsize=25)
plt.ylabel('-log10(P value)', fontsize=25)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

plt.xlim(-3.5,10)

plt.subplots_adjust(right=0.6,bottom=0.15,left=0.06,top=0.955)
plt.show()

