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
path2="/home/ali/Desktop/combined_data/combined_filtered_moderate.xlsx"
path3="/home/ali/Desktop/combined_data/combined_filtered_severe.xlsx"

# Let me read the data
data1=pd.read_excel(path1)
data2=pd.read_excel(path2)
data3=pd.read_excel(path3)
# Let me merge the combined data
z= pd.merge(data1,data2, on="Accession", how="outer")
z= pd.merge(z,data3, on="Accession", how="outer")

data_filter= z[["Accession","genes_1","genes_3","genes_4","M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2"
                    ,"S1_1","S1_2","S1_3","S2_1","S2_2","S2_3","S3_1","S3_2",
                   "C1_1","C1_2","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1","C5_2",
                "C6_1","C6_2"]]


# "M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2"

#"S1_1","S1_2","S1_3","S2_1","S2_2","S2_3","S3_1","S3_2",
# data_filter= z[["Accession","genes_1","genes_3","C1_1","C1_2","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1","C5_2","C6_1",
#                    "C6_2","M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2"
#                 ]]

data_filter["genes_1"].fillna(data_filter["genes_3"],inplace=True)
data_filter["genes_1"].fillna(data_filter["genes_4"],inplace=True)
# data_filter["genes_"].fillna(data_filter["genes_1"],inplace=True)
# data_filter["genes_3"].fillna(data_filter["genes_4"],inplace=True)

data_filter=data_filter[["Accession","genes_1","M1_1","M1_2","M1_3","M2_1","M3_1","M4_1","M4_2"
                    ,"S1_1","S1_2","S1_3","S2_1","S2_2","S2_3","S3_1","S3_2",
                   "C1_1","C1_2","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1","C5_2",
                "C6_1","C6_2"]]



a=data_filter.min(axis=1,skipna=True)
print(a.min())


# # Let me filter out 0.2 % threshold of the nan values and fill the remaining nans with small value of 0.001
data_filter= data_filter.dropna(thresh=len(data_filter.T)*(0.4))
data_filter.fillna(0.0001,inplace=True)
data_filter.reset_index(inplace=True)
data_filter=data_filter[data_filter.columns[1:]]
print(data_filter)

data_filter.to_excel("/home/ali/Desktop/ress.xlsx")


# Let me construct the critical data by taking mean values of each sample points between the patients!!

# data_filter["H3"]=(data_filter["H3"]+data_filter["H4"]+data_filter["H5"])/3
# data_filter["H4"]=(data_filter["H6"]+data_filter["H7"]+data_filter["H8"])/3
# data_filter["H5"]=(data_filter["H8"]+data_filter["H9"]+data_filter["H10"])/3

data_filter["Moderate1"]=(data_filter["M1_1"]+data_filter["M1_2"]+data_filter["M1_3"])/3
data_filter["Moderate2"]=(data_filter["M2_1"])/1
data_filter["Moderate3"]=(data_filter["M3_1"])/1
data_filter["Moderate4"]=(data_filter["M4_1"]+data_filter["M4_2"])/2

# data_filter["Moderate1"]=data_filter[["M1_1","M1_2","M1_3"]].mean(axis=1)
# data_filter["Moderate2"]=(data_filter["M2_1"])/1
# data_filter["Moderate3"]=(data_filter["M3_1"])/1

data_filter["Severe1"]=data_filter[["S1_1","S1_2","S1_3"]].mean(axis=1)
data_filter["Severe2"]=data_filter[["S2_1","S2_2","S2_3"]].mean(axis=1)
data_filter["Severe3"]=data_filter[["S3_1","S3_2"]].mean(axis=1)

data_filter["Critical1"]=data_filter[["C1_1","C1_2"]].mean(axis=1)
data_filter["Critical2"]=data_filter[["C2_1","C2_2","C2_3"]].mean(axis=1)
data_filter["Critical3"]=data_filter[["C3_1","C3_2"]].mean(axis=1)
data_filter["Critical4"]=data_filter[["C4_1","C4_2"]].mean(axis=1)
data_filter["Critical5"]=data_filter[["C5_1","C5_2"]].mean(axis=1)
data_filter["Critical6"]=data_filter[["C6_1","C6_2"]].mean(axis=1)
# data_filter["Critical7"]=data_filter[["C8_1","C8_2"]].mean(axis=1)


# data_filter["Severe1"]=(data_filter["S1_1"]+data_filter["S1_2"]+data_filter["S1_3"])/3
# data_filter["Severe2"]=(data_filter["S2_1"]+data_filter["S2_2"]+data_filter["S2_3"])/3
# data_filter["Severe3"]=(data_filter["S3_1"]+data_filter["S3_2"])/2
#
#
#
# data_filter["Critical1"]=(data_filter["C1_1"]+data_filter["C1_2"])/2
# data_filter["Critical2"]=(data_filter["C2_1"]+data_filter["C2_3"]+data_filter["C2_2"])/3
# data_filter["Critical3"]=(data_filter["C3_1"]+data_filter["C3_2"])/2
# data_filter["Critical4"]=(data_filter["C4_1"]+data_filter["C4_2"])/2
# data_filter["Critical5"]=(data_filter["C5_1"]+data_filter["C5_2"])/2
# data_filter["Critical6"]=(data_filter["C6_1"]+data_filter["C6_2"])/2
# data_filter["Critical7"]=(data_filter["C8_1"]+data_filter["C8_2"])/2
#
#
# # print(data_filter)
data_filter=data_filter[["Accession","genes_1","Critical1","Critical2","Critical3","Critical4","Critical5","Critical6","Moderate1","Moderate2","Moderate3","Moderate4","Severe1","Severe2","Severe3"]]

df_filter_only_values=data_filter[["Moderate1","Moderate2","Moderate3","Moderate4","Severe1","Severe2","Severe3"]]
df_filter2_only_values=data_filter[["Critical1","Critical2","Critical3","Critical4","Critical5","Critical6"]]

print(df_filter_only_values)
print(df_filter2_only_values)

transposed_data1=df_filter_only_values.T
transposed_data2=df_filter2_only_values.T


my_overall=pd.DataFrame({"Accession":data_filter["Accession"].tolist(),"Protein Name":data_filter["genes_1"].tolist(),"First Rep1":data_filter["Moderate1"],"First Rep2":data_filter["Moderate2"],
                         "First Rep3":data_filter["Moderate3"],"First Rep4":data_filter["Moderate4"],"First Rep5":data_filter["Severe1"],"First Rep6":data_filter["Severe2"],"First Rep7":data_filter["Severe3"],
                            "Second Rep1":data_filter["Critical1"],"Second Rep2":data_filter["Critical2"],"Second Rep3":data_filter["Critical3"],"Second Rep4":data_filter["Critical4"],
                            "Second Rep5":data_filter["Critical5"],"Second Rep6":data_filter["Critical6"]
                         })


print(my_overall)
# Not differentially expressed in moderate and severe, look at them!!
criticalp='/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2C/critical_vs_healthy/critical_healthy.xlsx'
severe_p="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2B/severe_vs_healthy/severe_healthy.xlsx"
moderate_p="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2A/moderate_vs_healthy/moderate_healthy.xlsx"
moderate=pd.read_excel(moderate_p)
severe=pd.read_excel(severe_p)
critical=pd.read_excel(criticalp)
my_overall=my_overall[(~my_overall["Accession"].isin(moderate["Accession"]))&(~my_overall["Accession"].isin(severe["Accession"]))]

# my_overall.to_excel("/home/a/Desktop/moderate_mutual.xlsx")

# # my_overall.to_excel("/home/a/Desktop/my_excell.xlsx")

my_overall['p_values']= pd.DataFrame(stats.ttest_ind(transposed_data1,transposed_data2,axis=0,equal_var=False)[1])

# my_overall.to_excel("/home/a/Desktop/res.xlsx")

my_overall['-log10 p value']=-np.log10(my_overall['p_values'])
# my_overall.to_excel("/home/a/Desktop/myoverall.xlsx")
my_overall.loc[my_overall['p_values']<=0.05,'significance']='+'
my_overall.loc[my_overall['p_values']>0.05,'significance']='-'

ppp=((my_overall["Second Rep1"]+my_overall["Second Rep2"]+my_overall["Second Rep3"]+my_overall["Second Rep4"]+my_overall["Second Rep5"]+my_overall["Second Rep6"])/6)/((my_overall["First Rep1"]+my_overall["First Rep2"]+my_overall["First Rep3"]+my_overall["First Rep4"]+my_overall["First Rep5"]+my_overall["First Rep6"]+my_overall["First Rep7"])/7)

my_overall['log2(Mean(First/Second))']= np.log2(ppp)
print(len(my_overall))
print(len(my_overall[my_overall["log2(Mean(First/Second))"]>0]))

my_overall=my_overall[my_overall["p_values"]!=0]

S= my_overall[my_overall["p_values"]<0.05]
print(len(S))

# my_overall.to_excel("/home/a/Desktop/result.xlsx")

my_overall.loc[(my_overall['significance']=='+')&(my_overall['log2(Mean(First/Second))']>0),'Sample downregulated (-) and upregulated (+) proteins']='+'
my_overall.loc[(my_overall['significance']=='+')&(my_overall['log2(Mean(First/Second))']<0),'Sample downregulated (-) and upregulated (+) proteins']='-'
my_overall['Sample downregulated (-) and upregulated (+) proteins']=my_overall['Sample downregulated (-) and upregulated (+) proteins'].fillna('Not Valid')
# my_overall.to_excel("/home/a/Desktop/resso.xlsx")

df1_Int_valid=my_overall[(my_overall['Sample downregulated (-) and upregulated (+) proteins']=='+')]
df1_Int_valid2=my_overall[my_overall['Sample downregulated (-) and upregulated (+) proteins']=='-']
df1_Int_notvalid=my_overall[my_overall['Sample downregulated (-) and upregulated (+) proteins']=='Not Valid']



# path_mild_up="/home/a/Desktop/regulation/mild/upregulated_ones.xlsx"
# path_mild_down="/home/a/Desktop/regulation/mild/downregulated_ones.xlsx"
# path_severe_up="/home/a/Desktop/regulation/severe/upregulated_ones.xlsx"
# path_severe_down="/home/a/Desktop/regulation/severe/downregulated_ones.xlsx"
#
# severe_u=pd.read_excel(path_severe_up)
# severe_d=pd.read_excel(path_severe_down)
#
# mild_u=pd.read_excel(path_mild_up)
# mild_d=pd.read_excel(path_mild_down)
# mild=pd.concat([mild_u,mild_d],axis=0,ignore_index=True)
# severe=pd.concat([severe_d,severe_u],axis=0,ignore_index=True)
# print(len(df1_Int_valid))
# print(len(df1_Int_valid2))
critical=pd.concat([df1_Int_valid,df1_Int_valid2],axis=0,ignore_index=True)

# critical_mild=critical[(critical["Accession"].isin(mild["Accession"]))&(~critical["Accession"].isin(severe["Accession"]))]
# critical_severe=critical[(~critical["Accession"].isin(mild["Accession"]))&(critical["Accession"].isin(severe["Accession"]))]
# common=critical[(critical["Accession"].isin(mild["Accession"]))&(critical["Accession"].isin(severe["Accession"]))]
# only_critical=critical[(~critical["Accession"].isin(mild["Accession"]))&(~critical["Accession"].isin(severe["Accession"]))]


df1_Int_valid.sort_values(by="log2(Mean(First/Second))",inplace=True,ignore_index=True,ascending=False)
df1_Int_valid2.sort_values(by="log2(Mean(First/Second))",inplace=True,ignore_index=True,ascending=True)



df1_Int_valid.to_excel("/home/ali/Desktop/upregulated_ones.xlsx")
df1_Int_valid2.to_excel("/home/ali/Desktop/downregulated_ones.xlsx")





# my_proteome_path="/home/a/Desktop/masaustu/COVID19_proteomics_database.txt"
#
# proteome=pd.read_csv(my_proteome_path,delimiter="\n")
#
# proteome=proteome[proteome.columns[0]].tolist()
#

# proteome_names=[]
# for each in proteome:
#     if each.startswith("#"):
#         continue
#     else:
#         proteome_names.append(each)
# print(proteome_names)
#

# count=0
# my_gene_coordinates={}
#
#
# critical.loc[critical["Protein Name"].isin(proteome_names),"Literature"]="+"
# critical.loc[~critical["Protein Name"].isin(proteome_names),"Literature"]="-"
#
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

# print(critical.columns)
#
# print(len(critical))
#
# print(len(my_gene_coordinates))
#

highest= critical.sort_values(by="log2(Mean(First/Second))",ascending=False)
# highest.to_excel("/home/a/Desktop/neybu.xlsx")

lowest=critical.sort_values(by="log2(Mean(First/Second))",ascending=True)
first_10={}

for j in range(0,9):
    x_axis=highest.iloc[j][18]
    y_axis=highest.iloc[j][16]
    first_10[highest.iloc[j][1]]=[x_axis,y_axis]
print(first_10)
least7={}
for j in range(0,6):
    x_axis=lowest.iloc[j][18]
    y_axis=lowest.iloc[j][16]
    least7[lowest.iloc[j][1]]=[x_axis,y_axis]
#

# print(highest.columns[0:6])
# print(first_10)


fig,ax=plt.subplots()


df1_Int_valid.to_excel("/home/ali/Desktop/moderate_red.xlsx")
df1_Int_valid2.to_excel("/home/ali/Desktop/moderate_blue.xlsx")
df1_Int_notvalid.to_excel("/home/ali/Desktop/moderate_gray.xlsx")





common_="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2D/venn_scheme/only_critical.xlsx"

common=pd.read_excel(common_)





# print(common)
#
#
df1_Int_valid_o=df1_Int_valid[df1_Int_valid["Accession"].isin(common["Accession"])]
df1_Int_valid_n=df1_Int_valid[~df1_Int_valid["Accession"].isin(common["Accession"])]

df1_Int_valid2_o=df1_Int_valid2[df1_Int_valid2["Accession"].isin(common["Accession"])]
df1_Int_valid2_n=df1_Int_valid2[~df1_Int_valid2["Accession"].isin(common["Accession"])]


plt.plot('log2(Mean(First/Second))','-log10 p value',data=df1_Int_valid_o,linestyle='',marker='*', markersize=5.5,alpha=1,color='purple',label="Commonly Found (n=9)")
plt.plot('log2(Mean(First/Second))','-log10 p value',data=df1_Int_valid2_o,linestyle='',marker='*', markersize=5.5,alpha=1,color='purple',label="Downregulated (n=16)")


print(df1_Int_valid_o)
print(df1_Int_valid2_o)

df1_Int_valid_o.to_excel("/home/ali/Desktop/upregulated_biomarkers.xlsx")
df1_Int_valid2_o.to_excel("/home/ali/Desktop/downregulated_biomarkers.xlsx")



# plt.plot('log2(Mean(First/Second))','-log10 p value',data=df1_Int_valid_n,linestyle='',marker='^', markersize=5.5,alpha=1,color='purple',label="Not found")
# plt.plot('log2(Mean(First/Second))','-log10 p value',data=df1_Int_valid2_n,linestyle='',marker='^', markersize=5.5,alpha=1,color='purple',label="Downregulated (n=16)")

# plt.legend(bbox_to_anchor=(0., 0.32, 1.65, 0.3),fontsize=20)

marker_list=["*","^"]
string=["Commonly found (n=8)","Not found (n=9)"]

# leg = plt.legend(loc=(1.03,0.45),title="Type",prop={'size':20},title_fontsize=20)
# ax.add_artist(leg)
#

# h = [plt.plot([],[], color="black", marker=marker_list[i],ls="")[0] for i in range(0,2)]

# print(h)
# leg2= plt.legend(handles=h, labels=string,loc=(1.06,0.38),title="Overlap in Critical-Healthy",prop={'size':20},title_fontsize=20)
# ax.add_artist(leg2)



plt.plot('log2(Mean(First/Second))','-log10 p value',data=df1_Int_notvalid,linestyle='',marker='o', markersize=5.5,alpha=0.50,color='gray')
# plt.axvline(0, color='black', linestyle='dashed', linewidth=1)

# plt.axvline(-0.3, color='black', linestyle='dashed', linewidth=1)
plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
plt.hlines(-np.log10(0.05),-8,8,colors='black',linestyles='dashed',linewidth=1)
# plt.xlim(-7,7)
plt.ylabel(("-log10(P value)"))


# print(critical[critical.columns[:]])
s= critical.iloc[1][0]
# print(s)
#
ax.annotate(list(first_10.keys())[0],(list(first_10.values())[0][0]-0.1,list(first_10.values())[0][1]+0.06),fontsize=12)
ax.annotate(list(first_10.keys())[1],(list(first_10.values())[1][0]-0.25,list(first_10.values())[1][1]+0.04),fontsize=12)
ax.annotate(list(first_10.keys())[2],(list(first_10.values())[2][0]-0.15,list(first_10.values())[2][1]+0.06),fontsize=12)
ax.annotate(list(first_10.keys())[3],(list(first_10.values())[3][0]-0.18,list(first_10.values())[3][1]+0.05),fontsize=12)
# ax.annotate(list(first_10.keys())[4],(list(first_10.values())[4][0]+0.15,list(first_10.values())[4][1]),fontsize=12)
# ax.annotate(list(first_10.keys())[5],(list(first_10.values())[5][0]+0.15,list(first_10.values())[5][1]),fontsize=12)
ax.annotate(list(first_10.keys())[8],(list(first_10.values())[8][0]-0.23,list(first_10.values())[8][1]+0.06),fontsize=12)


# # ax.annotate(list(least7.keys())[0],(list(least7.values())[0][0]-0.25,list(least7.values())[0][1]-0.2),fontsize=12)
# # ax.annotate(list(least7.keys())[1],(list(least7.values())[1][0]-0.2,list(least7.values())[1][1]+0.1),fontsize=12)
# # ax.annotate(list(least7.keys())[2],(list(least7.values())[2][0]+0.1,list(least7.values())[2][1]),fontsize=12)
# # ax.annotate(list(least7.keys())[3],(list(least7.values())[3][0]-0.6,list(least7.values())[3][1]+0.1),fontsize=12)
# # ax.annotate(list(least7.keys())[4],(list(least7.values())[4][0]-0.2,list(least7.values())[4][1]+0.2),fontsize=12)
# # ax.annotate(list(least7.keys())[5],(list(least7.values())[5][0]-0.15,list(least7.values())[5][1]+0.05),fontsize=12)
# # ax.annotate(list(least7.keys())[6],(list(least7.values())[6][0]-0.7,list(least7.values())[6][1]-0.15),fontsize=12)
# # # ax.annotate(list(least7.keys())[7],(list(least7.values())[7][0]-0.65,list(least7.values())[7][1]-0.175),fontsize=12)
# # # ax.annotate(list(least7.keys())[8],(list(least7.values())[8][0]-0.1,list(least7.values())[8][1]+0.04),fontsize=12)
# # ax.annotate(list(least7.keys())[9],(list(least7.values())[9][0]+0.05,list(least7.values())[9][1]+0.05),fontsize=12)
# # # ax.annotate(list(least7.keys())[10],(list(least7.values())[10][0]+0.05,list(least7.values())[10][1]-0.1),fontsize=12)
# # ax.annotate(list(least7.keys())[11],(list(least7.values())[11][0]+0.05,list(least7.values())[11][1]+0.03),fontsize=12)
# # # ax.annotate(list(least7.keys())[12],(list(least7.values())[12][0]+0.05,list(least7.values())[12][1]-0.1),fontsize=12)
#
#
plt.ylabel('-log10(P value)', fontsize=25)
plt.xlabel("log2(Critical/Severe-Moderate)",fontsize=25)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)


plt.xlim(-2,3)
plt.ylim(0,2.7)


plt.subplots_adjust(right=0.6,bottom=0.15,left=0.1,top=0.955)
plt.show()

# print(len(df1_Int_valid_o))
# print(len(df1_Int_valid2_o))
#
# print(len(df1_Int_valid_n))
# print(len(df1_Int_valid2_n))
# #