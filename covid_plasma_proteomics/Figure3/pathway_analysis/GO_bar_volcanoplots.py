import pandas
import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
from math import log
import matplotlib.pyplot as mplot
import pandas as pd

label_size = 8
mpl.rcParams['ytick.labelsize'] = label_size


data=pd.read_excel("/home/ali/Desktop/clusters/go_volcanoplots/bp_volcanoplots.xlsx")
print(data)


data=data[["Term Name Common","P value Common"]]

data["P value Common"]=-np.log10(data["P value Common"])
data.sort_values(by="P value Common",ascending=False,inplace=True)
data= data.head(n=10)
data.sort_values(by="P value Common",ascending=True,inplace=True)



# terms_=["Regulation of Insulin-like Growth Factor (IGF) transport"]+terms_
# print(terms_)
df = pandas.DataFrame(dict(graph=data["Term Name Common"].to_list(),
                           n=data["P value Common"].to_list()))
ind = np.arange(len(df))
width = 0.45

fig, ax = plt.subplots()

ax.barh(ind, df.n, width, color='blueviolet')
#ax.barh(ind + width, df.n, width, color='green', label='48 hours')
hfont={'fontname':'Arial'}

plt.ylabel('ylabel', **hfont)

plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
ax.set(yticks=ind + 1*width/8, yticklabels=df.graph, ylim=[2*width-1, len(df)])
# plt.xlim(0,13,1)
ax.set_title("Pathway Analysis",fontsize=20)
plt.xlabel('-log10(Adjusted p Value) ',fontsize=20)
# ax.legend()

from math import log
plt.axvline(1, color='lightgray', linestyle='dashed', linewidth=1)
plt.axvline(2, color='lightgray', linestyle='dashed', linewidth=1)
plt.axvline(3, color='lightgray', linestyle='dashed', linewidth=1)
plt.axvline(4, color='lightgray', linestyle='dashed', linewidth=1)
plt.axvline(5, color='lightgray', linestyle='dashed', linewidth=1)
plt.axvline(6, color='lightgray', linestyle='dashed', linewidth=1)
# plt.axvline(35, color='lightgray', linestyle='dashed', linewidth=1)
# plt.axvline(40, color='lightgray', linestyle='dashed', linewidth=1)
# plt.grid(color="lightgray")
# plt.text(0.35,6.5,"Only Critical/Healthy",fontsize=20)
plt.text(-0.5,11,"Critical-Severe-Moderate/Healthy",fontsize=20)

plt.show()

