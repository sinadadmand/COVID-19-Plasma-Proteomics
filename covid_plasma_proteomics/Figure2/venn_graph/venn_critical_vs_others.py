import pandas as pd
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

critical_path= "/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2C/critical_vs_healthy/critical_healthy.xlsx"
moderate_path="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2A/moderate_vs_healthy/moderate_healthy.xlsx"
severe_path="//home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2B/severe_vs_healthy/severe_healthy.xlsx"

# shu_path="/home/a/Documents/shu.xlsx"
# shen_path= "/home/a/Downloads/shen.xlsx"

# Let me write the path of the proteins found in

critical= pd.read_excel(critical_path,engine="openpyxl")
moderate=pd.read_excel(moderate_path)
severe=pd.read_excel(severe_path)

print(len(critical))
print(len(severe))

only_critical=len(critical[(~critical["Accession"].isin(severe["Accession"]))&(~critical["Accession"].isin(moderate["Accession"]))])
critical_moderate=len(critical[(~critical["Accession"].isin(severe["Accession"]))&(critical["Accession"].isin(moderate["Accession"]))])
critical_severe=len(critical[(critical["Accession"].isin(severe["Accession"]))&(~critical["Accession"].isin(moderate["Accession"]))])
severe_moderate=len(severe[(~severe["Accession"].isin(critical["Accession"]))&(severe["Accession"].isin(moderate["Accession"]))])
only_severe=len(severe[(~severe["Accession"].isin(critical["Accession"]))&(~severe["Accession"].isin(moderate["Accession"]))])
only_moderate=len(moderate[(~moderate["Accession"].isin(critical["Accession"]))&(~moderate["Accession"].isin(severe["Accession"]))])
common=len(critical[(critical["Accession"].isin(severe["Accession"]))&(critical["Accession"].isin(moderate["Accession"]))])

7

only_critical=only_critical
critical_moderate=critical_moderate
critical_severe=critical_severe



severe_moderate=severe_moderate
only_severe=only_severe
only_moderate=only_moderate
common=common




# critical_only=len(critical)-len(common)-len(critical_shu)-len(critical_shen)
# shen_only=len(shen)-len(critical_shen)-len(common)-len(shen_shu)
# shu_only=len(shu)-len(critical_shu)-len(common)-len(shen_shu)

out= venn3(subsets= (only_critical,only_severe,critical_severe,only_moderate,critical_moderate,severe_moderate,common), set_labels=("Critical","Severe","Moderate"),alpha=0.5)

for text in out.set_labels:
    text.set_fontsize(20)
for text in out.subset_labels:
    text.set_fontsize(20)
plt.show()
#
# print(common["Protein Name"])