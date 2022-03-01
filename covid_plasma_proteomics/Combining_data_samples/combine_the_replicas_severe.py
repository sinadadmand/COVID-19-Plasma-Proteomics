import pandas as pd
import numpy as np
from functools import reduce

my_pathlist=["/home/ali/Documents/projects/COVID_projects/COVID_data/Severe/S1/KUTTAM_20201011_ATS_S1-1.xlsx",
             "/home/ali/Documents/projects/COVID_projects/COVID_data/Severe/S1/KUTTAM_20201005_ATS_S1-2.xlsx",
             "/home/ali/Documents/projects/COVID_projects/COVID_data/Severe/S1/KUTTAM_20201029_ATS_S1-3.xlsx",
             "/home/ali/Documents/projects/COVID_projects/COVID_data/Severe/S2/KUTTAM_20201111_ATS_S2-1.xlsx",
             "/home/ali/Documents/projects/COVID_projects/COVID_data/Severe/S2/KUTTAM_20201106_ATS_S2-2.xlsx",
             "/home/ali/Documents/projects/COVID_projects/COVID_data/Severe/S2/KUTTAM_20201101_ATS_S2-3.xlsx",
             "/home/ali/Documents/projects/COVID_projects/COVID_data/Severe/S3/KUTTAM_20210420_ATS_S3-I(91).xlsx",
             "/home/ali/Documents/projects/COVID_projects/COVID_data/Severe/S3/KUTTAM_20210421_ATS_S3-I(92).xlsx"]




def gene_finder(my_pathlist):
    all_genes={}
    countt=0
    for i in my_pathlist:
        control=pd.read_excel(i)

        gene_description_list = control["Description"].tolist()
        GeneID_list = []
        index = []
        ind = 0
        for i in gene_description_list:
            counter = 0
            t = []
            z = []
            for j in i:
                t.append(j)
                if t[-1] == '=' and t[-2] == 'N' and t[-3] == 'G':
                    index.append(ind)
                    gene_name = i[counter + 1:-1]
                    for k in gene_name:
                        if k != ' ':
                            z.append(k)
                        else:
                            break
                    y = ''.join(z)
                    GeneID_list.append(y)
                counter = counter + 1
            ind = ind + 1
        # Let me find the indexes of the proteins that have the uncharacterized proteins!
        c2 = []
        for i in range(0, len(gene_description_list) - 1):
            if i not in index:
                c2.append(i)
        # 'c2' is a list of the indexes, now, let's form our gene list again
        one_in_desc2 = []
        for i in c2:
            if gene_description_list[i] == 'Description':
                one_in_desc2.append(i)

        count = 0

        for i in c2:
            GeneID_list.insert(i, 'Uncharachterized_Protein')
        gene = GeneID_list
        all_genes[countt]=gene
        countt=countt+1
    return all_genes
all_genes=gene_finder(my_pathlist)

# Let me read the data

control1=pd.read_excel(my_pathlist[0])
control2=pd.read_excel(my_pathlist[1])
control3=pd.read_excel(my_pathlist[2])
control4=pd.read_excel(my_pathlist[3])
control5=pd.read_excel(my_pathlist[4])
control6=pd.read_excel(my_pathlist[5])
control7=pd.read_excel(my_pathlist[6])
control8=pd.read_excel(my_pathlist[7])


# Let me add gene names to my data

control1["genes"]=pd.Series(all_genes[0])
control2["genes"]=pd.Series(all_genes[1])
control3["genes"]=pd.Series(all_genes[2])
control4["genes"]=pd.Series(all_genes[3])
control5["genes"]=pd.Series(all_genes[4])
control6["genes"]=pd.Series(all_genes[5])
control7["genes"]=pd.Series(all_genes[6])
control8["genes"]=pd.Series(all_genes[7])


# Let me do the protein based scaling!!
import statistics

list_control1=control1[(control1["genes"]=="ALB")|(control1["genes"]=="C3")|(control1["genes"]=="APOB")|(control1["genes"]=="A2M")|(control1["genes"]=="TF")|(control1["genes"]=="SERPINA1")|(control1["genes"]=="HP")]["# PSMs"].tolist()
list_control1=statistics.mean(list_control1)

list_control2=control2[(control2["genes"]=="ALB")|(control2["genes"]=="C3")|(control2["genes"]=="APOB")|(control2["genes"]=="A2M")|(control2["genes"]=="TF")|(control2["genes"]=="SERPINA1")|(control2["genes"]=="HP")]["# PSMs"].tolist()
list_control2=statistics.mean(list_control2)

list_control3=control3[(control3["genes"]=="ALB")|(control3["genes"]=="C3")|(control3["genes"]=="APOB")|(control3["genes"]=="A2M")|(control3["genes"]=="TF")|(control3["genes"]=="SERPINA1")|(control3["genes"]=="HP")]["# PSMs"].tolist()
list_control3=statistics.mean(list_control3)

list_control4=control4[(control4["genes"]=="ALB")|(control4["genes"]=="C3")|(control4["genes"]=="APOB")|(control4["genes"]=="A2M")|(control4["genes"]=="TF")|(control4["genes"]=="SERPINA1")|(control4["genes"]=="HP")]["# PSMs"].tolist()
list_control4=statistics.mean(list_control4)

list_control5=control5[(control5["genes"]=="ALB")|(control5["genes"]=="C3")|(control5["genes"]=="APOB")|(control5["genes"]=="A2M")|(control5["genes"]=="TF")|(control5["genes"]=="SERPINA1")|(control5["genes"]=="HP")]["# PSMs"].tolist()
list_control5=statistics.mean(list_control5)

list_control6=control6[(control6["genes"]=="ALB")|(control6["genes"]=="C3")|(control6["genes"]=="APOB")|(control6["genes"]=="A2M")|(control6["genes"]=="TF")|(control6["genes"]=="SERPINA1")|(control6["genes"]=="HP")]["# PSMs"].tolist()
list_control6=statistics.mean(list_control6)

list_control7=control7[(control7["genes"]=="ALB")|(control7["genes"]=="C3")|(control7["genes"]=="APOB")|(control7["genes"]=="A2M")|(control7["genes"]=="TF")|(control7["genes"]=="SERPINA1")|(control7["genes"]=="HP")]["# PSMs"].tolist()
list_control7=statistics.mean(list_control7)

list_control8=control8[(control8["genes"]=="ALB")|(control8["genes"]=="C3")|(control8["genes"]=="APOB")|(control8["genes"]=="A2M")|(control8["genes"]=="TF")|(control8["genes"]=="SERPINA1")|(control8["genes"]=="HP")]["# PSMs"].tolist()
list_control8=statistics.mean(list_control8)




# # Here, let me divide the PSMs of my data to median for the scaling!

control1["# PSMs"]=control1["# PSMs"].div(list_control1)
control2["# PSMs"]=control2["# PSMs"].div(list_control2)
control3["# PSMs"]=control3["# PSMs"].div(list_control3)
control4["# PSMs"]=control4["# PSMs"].div(list_control4)
control5["# PSMs"]=control5["# PSMs"].div(list_control5)
control6["# PSMs"]=control6["# PSMs"].div(list_control6)
control7["# PSMs"]=control7["# PSMs"].div(list_control7)
control8["# PSMs"]=control8["# PSMs"].div(list_control8)

# Let me merge the data!

z= pd.merge(control1,control2, on="Accession", how="outer")
z= pd.merge(z,control3, on="Accession", how="outer")
z= pd.merge(z,control4, on="Accession", how="outer")
z= pd.merge(z,control5, on="Accession", how="outer")
z= pd.merge(z,control6, on="Accession", how="outer")
z= pd.merge(z,control7, on="Accession", how="outer")
z= pd.merge(z,control8, on="Accession", how="outer")

new_frame= z[z.columns[[3,9,26,17,34]]]

new_frame.columns=["Accession","S1_1","S1_3","S2_2","S3_1","S1_2","S2_1","S2_3","S3_2","genes1_1","genes1_3","genes2_2","genes3_1","genes1_2",
                   "genes2_1","genes2_3","genes3_2"]


new_frame["genes1_1"].fillna(new_frame["genes1_2"],inplace=True)
new_frame["genes1_1"].fillna(new_frame["genes1_3"],inplace=True)
new_frame["genes1_1"].fillna(new_frame["genes2_1"],inplace=True)
new_frame["genes1_1"].fillna(new_frame["genes2_2"],inplace=True)
new_frame["genes1_1"].fillna(new_frame["genes2_3"],inplace=True)
new_frame["genes1_1"].fillna(new_frame["genes3_1"],inplace=True)
new_frame["genes1_1"].fillna(new_frame["genes3_2"],inplace=True)


new_frame=new_frame[["Accession","genes1_1","S1_1","S1_2","S1_3","S2_1","S2_2","S2_3","S3_1","S3_2"]]


new_frame.to_excel("/home/ali/Desktop/combined_filtered_severe.xlsx")
