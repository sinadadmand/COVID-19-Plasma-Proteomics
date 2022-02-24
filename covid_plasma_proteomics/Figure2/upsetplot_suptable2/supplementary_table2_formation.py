import pandas as pd

# Let me gather all the data frames from the literature I found!!!

# Shen at al
path_shen="/home/ali/Desktop/literature_proteins/literature_proteins/shen.xlsx"
shen=pd.read_excel(path_shen)

# Shu at al
path_shu="/home/ali/Desktop/literature_proteins/literature_proteins/shu.xlsx"
shu=pd.read_excel(path_shu)

# Park at al
path_park= "/home/ali/Desktop/literature_proteins/literature_proteins/park_groups.xlsx"
park=pd.read_excel(path_park)

# Demichev at al
path_demichev= "/home/ali/Desktop/literature_proteins/literature_proteins/demichev_groups.xlsx"
demichev=pd.read_excel(path_demichev)

# Geyer et al 1

path_geyer1="/home/ali/Desktop/literature_proteins/literature_proteins/geyers_paper_significantly_difference_first_day_of_sampling_COVIDvs_Healthy.xlsx"
path_geyer2="/home/ali/Desktop/literature_proteins/literature_proteins/geyers_paper_significantly_difference_highest_signal_COVID_vs_Healthy.xlsx"
path_geyer3="/home/ali/Desktop/literature_proteins/literature_proteins/geyers_paper_significantly_difference_highest_signal_vs_first_timepoint_COVID_vs_COVID.xlsx"


# This experiment
gene_list_messner=["A1BG","ACTB","C1R","C1S","C8A","CD14","CFB","CFH","CFI","CRP","ACTG1","FGA","FGB","FGG","HP","ITIH3",
           "ITIH4","LBP","LGALS3BP","LRG1","SAA1","SAA2","SERPINA10","ALB","APOA1","APOC1","GSN","TF"]



# Let me rearrange my experimental data!

critical_ones= "/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2C/critical_vs_healthy/critical_healthy.xlsx"
severe_ones="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2B/severe_vs_healthy/severe_healthy.xlsx"
moderate_ones="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure2A/moderate_vs_healthy/moderate_healthy.xlsx"

crit= pd.read_excel(critical_ones)
sev=pd.read_excel(severe_ones)
mod=pd.read_excel(moderate_ones)

crit=crit[["Accession","Protein Name"]]
sev=sev[["Accession","Protein Name"]]
mod=mod[["Accession","Protein Name"]]

merged1=list(pd.concat([crit,sev,mod])["Protein Name"].unique())
merged2=list(pd.concat([crit,sev,mod])["Accession"].unique())

data_exp= pd.DataFrame({"Accession":merged2,"Protein Name":merged1})

# Shen data
shen["Protein Name"]= [i.split(",")[0] for i in shen["Names"]]

# Demcihev data

'''
In the data, there are severity p values, here, I have got the proteins which have p values less than 0.05
'''

demichev=demichev[demichev["Severity.Association.P.Value"]<0.05]

demichev.reset_index(inplace=True)
demichev=demichev[demichev.columns[1:]]

demichev=demichev[["Measurement","Protein"]]
print(demichev)

demichev["Measurement"]=[i.split(";") for i in demichev["Measurement"]]

print(demichev)
demichev= demichev.explode("Measurement")
print(demichev)
demichev.reset_index(inplace=True)
demichev=demichev[demichev.columns[1:]]
demichev.columns=["Protein Name","Description"]


# Geyer et al

geyer1= pd.read_excel(path_geyer1)
geyer2= pd.read_excel(path_geyer2)
geyer3= pd.read_excel(path_geyer3)

data_geyer=pd.concat([geyer1,geyer2,geyer3])

data_geyer["Gene names"]=[str(i).split(";") for i in data_geyer["Gene names"]]
data_geyer= data_geyer.explode("Gene names")


# data_geyer=data_geyer["Gene names"].dropna()

a= list(set(data_geyer["Gene names"].to_list()))

print(a)

geyer=pd.DataFrame({"Protein Name":a,"name":a})



# Park at al

print(park["Protein Name"])

park= park.dropna(subset=["Protein Name"])
park.reset_index(inplace=True)
park=park[park.columns[1:]]
print(len(park))
park["Protein Name"]=[i.split(";") for i in park["Protein Name"]]
park= demichev.explode("Protein Name")


# Let me format the Messener data!

mess= pd.DataFrame({"Protein Name":gene_list_messner,"Protein":gene_list_messner})

# # Let me now check the data format
print(demichev)
print(data_exp)
print(park)
print(shu)
print(shen)
print(mess)

all= pd.concat([demichev,data_exp,park,shu,shen,mess,geyer])

# all.to_excel("/home/a/Desktop/all.xlsx")


names= list(set(all["Protein Name"].dropna().to_list()))

print(len(names))
print(names)



# Let me create each column now!!

# Shen data list;
#
shen_list=[]
shu_list=[]
exp_list=[]
park_list=[]
demichev_list=[]
mess_list=[]
geyer_list=[]

for each in names:
    if each in shen["Protein Name"].to_list():
        shen_list.append("+")
    else:
        shen_list.append("-")

for each in names:
    if each in shu["Protein Name"].to_list():
        shu_list.append("+")
    else:
        shu_list.append("-")

for each in names:
    if each in data_exp["Protein Name"].to_list():
        exp_list.append("+")
    else:
        exp_list.append("-")

for each in names:
    if each in park["Protein Name"].to_list():
        park_list.append("+")
    else:
        park_list.append("-")


for each in names:
    if each in demichev["Protein Name"].to_list():
        demichev_list.append("+")
    else:
        demichev_list.append("-")


for each in names:
    if each in mess["Protein Name"].to_list():
        mess_list.append("+")
    else:
        mess_list.append("-")

for each in names:
    if each in geyer["Protein Name"].to_list():
        geyer_list.append("+")
    else:
        geyer_list.append("-")



data= pd.DataFrame({"Protein Names":names,"Geyer at al":geyer_list,"Shu at al":shu_list,"Shen at al":shen_list,"Park at al":park_list,"Demichev at al.":demichev_list,"Messner at al.":mess_list,"This experiment":exp_list})

data.to_excel("/home/ali/Desktop/SupplementaryTable2.xlsx")
