import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# path_c="/home/ali/Documents/projects/COVID_projects/Combined_data/combined_filtered_critical.xlsx"
# path_m="/home/ali/Documents/projects/COVID_projects/Combined_data/combined_filtered_moderate.xlsx"
# path_s="/home/ali/Documents/projects/COVID_projects/Combined_data/combined_filtered_severe.xlsx"
# path_h="/home/ali/Documents/projects/COVID_projects/Combined_data/combined_filtered_healthy.xlsx"
#
# c=pd.read_excel(path_c)
# m=pd.read_excel(path_m)
# s=pd.read_excel(path_s)
# h=pd.read_excel(path_h)
#
#
# merged_=pd.merge(c,m,on="Accession",how="outer")
# merged_=pd.merge(merged_,s,on="Accession",how="outer")
# merged_=pd.merge(merged_,h,on="Accession",how="outer")
#
# merged_.dropna(inplace=True)
#
# merged_=merged_[merged_.columns[3:]]
#
# merged_.to_excel("/home/ali/Desktop/bak.xlsx")
#

pdata= "/home/ali/Documents/correlatiion_Feb17.xlsx"

data=pd.read_excel(pdata)

print(len(data))

# data= data[data.columns[::-1]]

# correlate=data.corr(method='pearson')



correlate=pd.read_excel("/home/ali/Documents/correlation.xlsx")





# correlate=correlate.where(np.tril(np.ones(correlate.shape)).astype(np.bool))

# print(correlate)
sns.set(font_scale=1.7)

# correlate.to_excel("/home/ali/Desktop/correlation_result.xlsx")


print(correlate)


correlate=correlate[correlate.columns[1:]]

# a=correlate.columns

correlate= correlate[correlate.columns[::-1]]
correlate= correlate.iloc[::-1]

print(correlate)
#cbar_kws={'shrink': 0.3,"location":"bottom","label":"Pearson corr. coef."}

ax= sns.heatmap(correlate,cmap="coolwarm",annot=False,yticklabels=False,xticklabels=False, cbar_kws = dict(use_gridspec=False,location="bottom",shrink=0.3,label="Pearson corr. coef."),vmin=0.84)


ax.hlines([5], *ax.get_xlim(),color="black")
ax.hlines([14], *ax.get_xlim(),color="black")
ax.hlines([25], *ax.get_xlim(),color="black")


ax.vlines([5], *ax.get_xlim(),color="black")
ax.vlines([14], *ax.get_xlim(),color="black")
ax.vlines([25], *ax.get_xlim(),color="black")


# plt.text(0, -2,"Healthy", fontsize=15)
# plt.text(-2, 4,"Healthy", fontsize=15,rotation=90)
#
# plt.text(7, -2,"Moderate", fontsize=14)
# plt.text(-2, 13,"Moderate", fontsize=15,rotation=90)
#
#
# plt.text(15, -2,"Severe", fontsize=15)
#
#
# plt.text(25, -2,"Critical", fontsize=15)


plt.subplots_adjust(left=0.06,right=0.44,bottom=0)

plt.show()
