import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


pdata= "/home/ali/Documents/correlatiion_Feb17.xlsx"

data=pd.read_excel(pdata)

print(len(data))



correlate=pd.read_excel("/home/ali/Documents/correlation.xlsx")





sns.set(font_scale=1.7)



print(correlate)


correlate=correlate[correlate.columns[1:]]


correlate= correlate[correlate.columns[::-1]]
correlate= correlate.iloc[::-1]

print(correlate)

ax= sns.heatmap(correlate,cmap="coolwarm",annot=False,yticklabels=False,xticklabels=False, cbar_kws = dict(use_gridspec=False,location="bottom",shrink=0.3,label="Pearson corr. coef."),vmin=0.84)


ax.hlines([5], *ax.get_xlim(),color="black")
ax.hlines([14], *ax.get_xlim(),color="black")
ax.hlines([25], *ax.get_xlim(),color="black")


ax.vlines([5], *ax.get_xlim(),color="black")
ax.vlines([14], *ax.get_xlim(),color="black")
ax.vlines([25], *ax.get_xlim(),color="black")


plt.subplots_adjust(left=0.06,right=0.44,bottom=0)

plt.show()
