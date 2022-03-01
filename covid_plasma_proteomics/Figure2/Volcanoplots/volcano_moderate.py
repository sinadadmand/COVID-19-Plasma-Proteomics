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


path1="/home/ali/Desktop/combined_data/combined_filtered_moderate.xlsx"
path2="/home/ali/Desktop/combined_data/combined_filtered_healthy.xlsx"

# Let me read the data
data1=pd.read_excel(path1)
data2=pd.read_excel(path2)
# data3=pd.read_excel(path3)
# Let me merge the combined data
z= pd.merge(data1,data2, on="Accession", how="outer")
# z= pd.merge(z,data3, on="Accession", how="outer")

data_filter= z[["Accession","genes_2","genes_3","M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2",
                   "H1","H2","H3","H4","H5"]]



data_filter["genes_2"].fillna(data_filter["genes_3"],inplace=True)

data_filter=data_filter[["Accession","genes_2","M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2",
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


data_filter["Moderate1"]=(data_filter["M1_1"]+data_filter["M1_2"]+data_filter["M1_3"])/3
data_filter["Moderate2"]=(data_filter["M2_1"])/1
data_filter["Moderate3"]=(data_filter["M3_1"])/1
data_filter["Moderate4"]=(data_filter["M4_1"]+data_filter["M4_2"])/2

data_filter=data_filter[["Accession","genes_2","H1","H2","H3","H4","H5","Moderate1","Moderate2","Moderate3","Moderate4"]]

df_filter_only_values=data_filter[["Moderate1","Moderate2","Moderate3","Moderate4"]]
df_filter2_only_values=data_filter[["H1","H2","H3","H4","H5"]]

print(df_filter_only_values)
print(df_filter2_only_values)

transposed_data1=df_filter_only_values.T
transposed_data2=df_filter2_only_values.T


my_overall=pd.DataFrame({"Accession":data_filter["Accession"].tolist(),"Protein Name":data_filter["genes_2"].tolist(),"First Rep1":data_filter["Moderate1"],"First Rep2":data_filter["Moderate2"],
                         "First Rep3":data_filter["Moderate3"],"First Rep4":data_filter["Moderate4"],
                            "Second Rep1":data_filter["H1"],"Second Rep2":data_filter["H2"],"Second Rep3":data_filter["H3"],"Second Rep4":data_filter["H4"],
                            "Second Rep5":data_filter["H5"]
                         })


my_overall['p_values']= pd.DataFrame(stats.ttest_ind(transposed_data1,transposed_data2,axis=0,equal_var=False)[1])


my_overall['-log10 p value']=-np.log10(my_overall['p_values'])
my_overall.loc[my_overall['p_values']<=0.05,'significance']='+'
my_overall.loc[my_overall['p_values']>0.05,'significance']='-'

ppp=((my_overall["First Rep1"]+my_overall["First Rep2"]+my_overall["First Rep3"]+my_overall["First Rep4"])/4)/((my_overall["Second Rep1"]+my_overall["Second Rep2"]+my_overall["Second Rep3"]+my_overall["Second Rep4"]+my_overall["Second Rep5"])/5)

my_overall['log2(Mean(First/Second))']= np.log2(ppp)
print(len(my_overall))
print(len(my_overall[my_overall["log2(Mean(First/Second))"]>0]))

my_overall=my_overall[my_overall["p_values"]!=0]

S= my_overall[my_overall["p_values"]<0.05]
print(len(S))


my_overall.loc[(my_overall['significance']=='+')&(my_overall['log2(Mean(First/Second))']>0),'Sample downregulated (-) and upregulated (+) proteins']='+'
my_overall.loc[(my_overall['significance']=='+')&(my_overall['log2(Mean(First/Second))']<-0),'Sample downregulated (-) and upregulated (+) proteins']='-'
my_overall['Sample downregulated (-) and upregulated (+) proteins']=my_overall['Sample downregulated (-) and upregulated (+) proteins'].fillna('Not Valid')

df1_Int_valid=my_overall[(my_overall['Sample downregulated (-) and upregulated (+) proteins']=='+')]
df1_Int_valid2=my_overall[my_overall['Sample downregulated (-) and upregulated (+) proteins']=='-']
df1_Int_notvalid=my_overall[my_overall['Sample downregulated (-) and upregulated (+) proteins']=='Not Valid']


path_mild_up="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2C/critical_vs_healthy/upregulated_ones.xlsx"
path_mild_down="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2C/critical_vs_healthy/downregulated_ones.xlsx"
path_severe_up="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2B/severe_vs_healthy/upregulated_ones.xlsx"
path_severe_down="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2B/severe_vs_healthy/downregulated_ones.xlsx"

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

for j in range(0,60):
    x_axis=highest.iloc[j][14]
    y_axis=highest.iloc[j][12]
    first_10[highest.iloc[j][1]]=[x_axis,y_axis]

print(first_10)

least7={}
fig,ax=plt.subplots()

plt.plot('log2(Mean(First/Second))','-log10 p value',data=critical_severe,linestyle='',marker='o', markersize=5.5,alpha=1,color='blue',label="moderate-severe (n=11)")
plt.plot('log2(Mean(First/Second))','-log10 p value',data=common,linestyle='',marker='o', markersize=5.5,alpha=1,color='purple',label="common (n=26)")
plt.plot('log2(Mean(First/Second))','-log10 p value',data=only_critical,linestyle='',marker='o', markersize=5.5,alpha=1,color='darkturquoise',label="only moderate (n=23)")
plt.plot('log2(Mean(First/Second))','-log10 p value',data=critical_mild,linestyle='',marker='o', markersize=5.5,alpha=1,color='pink',label="moderate-critical (n=28)")

marker_list=["o","v"]
string=["No previous data","Previously found"]

h = [plt.plot([],[], color="black", marker=marker_list[i],ls="")[0] for i in range(0,2)]




plt.plot('log2(Mean(First/Second))','-log10 p value',data=df1_Int_notvalid,linestyle='',marker='o', markersize=5.5,alpha=0.50,color='gray')

plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
plt.hlines(-np.log10(0.05),-10,10,colors='black',linestyles='dashed',linewidth=1)
plt.xlim(-8,11)
plt.ylabel(("-log10(P value)"))
plt.xlabel(("log2(fold change)"))


s= critical.iloc[1][0]


#
#
ax.annotate(list(first_10.keys())[0],(list(first_10.values())[0][0]+0.15,list(first_10.values())[0][1]),fontsize=12)
ax.annotate(list(first_10.keys())[1],(list(first_10.values())[1][0]+0.15,list(first_10.values())[1][1]-0.05),fontsize=12)
ax.annotate(list(first_10.keys())[2],(list(first_10.values())[2][0]-0.25,list(first_10.values())[2][1]+0.1),fontsize=12)
ax.annotate(list(first_10.keys())[3],(list(first_10.values())[3][0],list(first_10.values())[3][1]+0.5),fontsize=12)
ax.annotate(list(first_10.keys())[4],(list(first_10.values())[4][0]-0.07,list(first_10.values())[4][1]-0.2),fontsize=12)
ax.annotate(list(first_10.keys())[5],(list(first_10.values())[5][0],list(first_10.values())[5][1]+0.3),fontsize=12)
ax.annotate(list(first_10.keys())[6],(list(first_10.values())[6][0]+0.1,list(first_10.values())[6][1]),fontsize=12)
ax.annotate(list(first_10.keys())[7],(list(first_10.values())[7][0]-0.05,list(first_10.values())[7][1]-0.28),fontsize=12)
# ax.annotate(list(first_10.keys())[8],(list(first_10.values())[8][0]+0.1,list(first_10.values())[8][1]),fontsize=12)
ax.annotate(list(first_10.keys())[9],(list(first_10.values())[9][0]+0.15,list(first_10.values())[9][1]),fontsize=12)
# ax.annotate(list(first_10.keys())[10],(list(first_10.values())[10][0]+0.1,list(first_10.values())[10][1]),fontsize=12)
ax.annotate(list(first_10.keys())[11],(list(first_10.values())[11][0]+0.1,list(first_10.values())[11][1]),fontsize=12)
#
# # #
# # #
ax.annotate(list(first_10.keys())[12],(list(first_10.values())[12][0]-0.1,list(first_10.values())[12][1]+0.1),fontsize=12)
# #
# #
ax.annotate(list(first_10.keys())[13],(list(first_10.values())[13][0]+0.1,list(first_10.values())[13][1]+0.2),fontsize=12)
# #
# ax.annotate(list(first_10.keys())[14],(list(first_10.values())[14][0]+0.1,list(first_10.values())[14][1]),fontsize=12)
#
ax.annotate(list(first_10.keys())[15],(list(first_10.values())[15][0],list(first_10.values())[15][1]+0.15),fontsize=12)
# ax.annotate(list(first_10.keys())[16],(list(first_10.values())[16][0],list(first_10.values())[16][1]+0.15),fontsize=12)

ax.annotate(list(first_10.keys())[17],(list(first_10.values())[17][0]+0.1,list(first_10.values())[17][1]+0.05),fontsize=12)
ax.annotate(list(first_10.keys())[19],(list(first_10.values())[19][0]+0.1,list(first_10.values())[19][1]),fontsize=12)
#

# ax.annotate(list(first_10.keys())[22],(list(first_10.values())[22][0]+0.15,list(first_10.values())[22][1]-0.1),fontsize=12)
ax.annotate(list(first_10.keys())[23],(list(first_10.values())[23][0]+0.15,list(first_10.values())[23][1]),fontsize=12)
ax.annotate(list(first_10.keys())[24],(list(first_10.values())[24][0]+0.15,list(first_10.values())[24][1]-0.12),fontsize=12)
# ax.annotate(list(first_10.keys())[26],(list(first_10.values())[26][0]+0.15,list(first_10.values())[26][1]),fontsize=12)

ax.annotate(list(first_10.keys())[31],(list(first_10.values())[31][0]+0.1,list(first_10.values())[31][1]+0.04),fontsize=12)

# # ax.annotate(list(first_10.keys())[47],(list(first_10.values())[47][0]+0.1,list(first_10.values())[47][1]-0.05),fontsize=12)
ax.annotate(list(first_10.keys())[48],(list(first_10.values())[48][0]-0.25,list(first_10.values())[48][1]+0.07),fontsize=12)
ax.annotate(list(first_10.keys())[49],(list(first_10.values())[49][0]-0.5,list(first_10.values())[49][1]+0.08),fontsize=12)
ax.annotate(list(first_10.keys())[50],(list(first_10.values())[50][0]-0.5,list(first_10.values())[50][1]+0.09),fontsize=12)
ax.annotate(list(first_10.keys())[51],(list(first_10.values())[51][0]+0.1,list(first_10.values())[51][1]+0.05),fontsize=12)
ax.annotate(list(first_10.keys())[52],(list(first_10.values())[52][0]-0.95,list(first_10.values())[52][1]-0.1),fontsize=12)
ax.annotate(list(first_10.keys())[53],(list(first_10.values())[53][0]-0.8,list(first_10.values())[53][1]+0.07),fontsize=12)
ax.annotate(list(first_10.keys())[54],(list(first_10.values())[54][0]+0.14,list(first_10.values())[54][1]),fontsize=12)
ax.annotate(list(first_10.keys())[55],(list(first_10.values())[55][0]+0.14,list(first_10.values())[55][1]),fontsize=12)
ax.annotate(list(first_10.keys())[56],(list(first_10.values())[56][0]+0.1,list(first_10.values())[56][1]+0.04),fontsize=12)
ax.annotate(list(first_10.keys())[57],(list(first_10.values())[57][0]+0.1,list(first_10.values())[57][1]+0.04),fontsize=12)
ax.annotate(list(first_10.keys())[58],(list(first_10.values())[58][0]-0.25,list(first_10.values())[58][1]+0.08),fontsize=12)
ax.annotate(list(first_10.keys())[59],(list(first_10.values())[59][0]-0.55,list(first_10.values())[59][1]+0.08),fontsize=12)



plt.xlabel('log2(Moderate/Healthy)', fontsize=25)
plt.ylabel('-log10(P value)', fontsize=25)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

plt.xlim(-6.7,6.7)

plt.subplots_adjust(right=0.6,bottom=0.15,left=0.06,top=0.955)
plt.show()

