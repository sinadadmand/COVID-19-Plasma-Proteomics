import pandas as pd
import numpy as np
from functools import reduce

# my_pathlist=["/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C1/KUTTAM_20201013_ATS_C1-1.xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C1/KUTTAM_20201029_ATS_C1-2.xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C1/KUTTAM_20201007_ATS_C1-3.xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C2/KUTTAM_20201112_ATS_C2-1.xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C2/KUTTAM_20201108_ATS_C2-2.xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C2/KUTTAM_20201114_ATS_C2-3.xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C3/KUTTAM_20201212_ATS_C3-1.xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C3/KUTTAM_20201213_ATS_C3-2.xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C4/KUTTAM_20210420_ATS_C4-I(91).xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C4/KUTTAM_20210421_ATS_C4-I(92).xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C5/KUTTAM_20210410_ATS_C5-E.xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C5/KUTTAM_20210412_ATS_C5-I.xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C6/KUTTAM_20210429_ATS_C6-E.xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C6/KUTTAM_20210430_ATS_C6-I.xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C8/KUTTAM_20210504_ATS_C8-E.xlsx",
#              "/home/a/Desktop/covid_project_June10/PDExcels/Criticals/C8/KUTTAM_20210519_ATS_C8-I.xlsx"]
#/home/ali/Documents/projects/COVID_projects/COVID_data

my_pathlist=[
    # "/home/a/Desktop/covid_project_June10/PDExcels/Healthy/KUTTAM_20201002_ATS_H1.xlsx",
    #          "/home/a/Desktop/covid_project_June10/PDExcels/Healthy/KUTTAM_20201010_ATS_H2.xlsx",
             "/home/ali/Documents/projects/COVID_projects/COVID_data/Healthy/KUTTAM_20210522_ATS_H3.xlsx",
             "/home/ali/Documents/projects/COVID_projects/COVID_data/Healthy/KUTTAM_20210524_ATS_H4.xlsx",
             "/home/ali/Documents/projects/COVID_projects/COVID_data/Healthy/KUTTAM_20210526_ATS_H5.xlsx",
             "/home/ali/Documents/projects/COVID_projects/COVID_data/Healthy/KUTTAM_20210514_ATS_H6.xlsx",
            # "/home/ali/Downloads/KUTTAM_20201010_ATS_35H_newH5.xlsx",
            "/home/ali/Documents/projects/COVID_projects/COVID_data/Healthy/KUTTAM_20200822_ATS_1H_new.xlsx"
    # "/home/ali/Documents/projects/COVID_projects/COVID_data/Healthy/KUTTAM_20210207_ATS_H9.xlsx",
             # "/home/a/Desktop/covid_project_June10/PDExcels/Healthy/KUTTAM_20201031_ATS_H7.xlsx",
             # "/home/a/Desktop/covid_project_June10/PDExcels/Healthy/KUTTAM_20200822_ATS_H8.xlsx",
             # "/home/ali/Documents/projects/COVID_projects/COVID_data/Healthy/KUTTAM_20210511_ATS_H10.xlsx"]

]




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
# control6=pd.read_excel(my_pathlist[5])
# control7=pd.read_excel(my_pathlist[6])
# control8=pd.read_excel(my_pathlist[7])
# control9=pd.read_excel(my_pathlist[8])
# control10=pd.read_excel(my_pathlist[9])
# control11=pd.read_excel(my_pathlist[10])
# control12=pd.read_excel(my_pathlist[11])
# control13=pd.read_excel(my_pathlist[12])
# control14=pd.read_excel(my_pathlist[13])
# control15=pd.read_excel(my_pathlist[14])
# control16=pd.read_excel(my_pathlist[15])


# Let me add gene names to my data

control1["genes"]=pd.Series(all_genes[0])
control2["genes"]=pd.Series(all_genes[1])
control3["genes"]=pd.Series(all_genes[2])
control4["genes"]=pd.Series(all_genes[3])
control5["genes"]=pd.Series(all_genes[4])
# control6["genes"]=pd.Series(all_genes[5])
# control7["genes"]=pd.Series(all_genes[6])
# control8["genes"]=pd.Series(all_genes[7])
# control9["genes"]=pd.Series(all_genes[8])
# control10["genes"]=pd.Series(all_genes[9])
# control11["genes"]=pd.Series(all_genes[10])
# control12["genes"]=pd.Series(all_genes[11])
# control13["genes"]=pd.Series(all_genes[12])
# control14["genes"]=pd.Series(all_genes[13])
# control15["genes"]=pd.Series(all_genes[14])
# control16["genes"]=pd.Series(all_genes[15])


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

# list_control6=control6[(control6["genes"]=="ALB")|(control6["genes"]=="C3")|(control6["genes"]=="APOB")|(control6["genes"]=="A2M")|(control6["genes"]=="TF")|(control6["genes"]=="SERPINA1")|(control6["genes"]=="HP")]["# PSMs"].tolist()
# list_control6=statistics.mean(list_control6)



# list_control7=control7[(control7["genes"]=="ALB")|(control7["genes"]=="C3")|(control7["genes"]=="APOB")|(control7["genes"]=="A2M")|(control7["genes"]=="TF")|(control7["genes"]=="SERPINA1")|(control7["genes"]=="HP")]["# PSMs"].tolist()
# list_control7=statistics.mean(list_control7)
#
# list_control8=control8[(control8["genes"]=="ALB")|(control8["genes"]=="C3")|(control8["genes"]=="APOB")|(control8["genes"]=="A2M")|(control8["genes"]=="TF")|(control8["genes"]=="SERPINA1")|(control8["genes"]=="HP")]["# PSMs"].tolist()
# list_control8=statistics.mean(list_control8)
#
# list_control9=control9[(control9["genes"]=="ALB")|(control9["genes"]=="C3")|(control9["genes"]=="APOB")|(control9["genes"]=="A2M")|(control9["genes"]=="TF")|(control9["genes"]=="SERPINA1")|(control9["genes"]=="HP")]["# PSMs"].tolist()
# list_control9=statistics.mean(list_control9)
#
# list_control10=control10[(control10["genes"]=="ALB")|(control10["genes"]=="C3")|(control10["genes"]=="APOB")|(control10["genes"]=="A2M")|(control10["genes"]=="TF")|(control10["genes"]=="SERPINA1")|(control10["genes"]=="HP")]["# PSMs"].tolist()
# list_control10=statistics.mean(list_control10)
#
# list_control11=control11[(control11["genes"]=="ALB")|(control11["genes"]=="C3")|(control11["genes"]=="APOB")|(control11["genes"]=="A2M")|(control11["genes"]=="TF")|(control11["genes"]=="SERPINA1")|(control11["genes"]=="HP")]["# PSMs"].tolist()
# list_control11=statistics.mean(list_control11)
#
# list_control12=control12[(control12["genes"]=="ALB")|(control12["genes"]=="C3")|(control12["genes"]=="APOB")|(control12["genes"]=="A2M")|(control12["genes"]=="TF")|(control12["genes"]=="SERPINA1")|(control12["genes"]=="HP")]["# PSMs"].tolist()
# list_control12=statistics.mean(list_control12)
#
# list_control13=control13[(control13["genes"]=="ALB")|(control13["genes"]=="C3")|(control13["genes"]=="APOB")|(control13["genes"]=="A2M")|(control13["genes"]=="TF")|(control13["genes"]=="SERPINA1")|(control13["genes"]=="HP")]["# PSMs"].tolist()
# list_control13=statistics.mean(list_control13)
#
# list_control14=control14[(control14["genes"]=="ALB")|(control14["genes"]=="C3")|(control14["genes"]=="APOB")|(control14["genes"]=="A2M")|(control14["genes"]=="TF")|(control14["genes"]=="SERPINA1")|(control14["genes"]=="HP")]["# PSMs"].tolist()
# list_control14=statistics.mean(list_control14)
#
# list_control15=control15[(control15["genes"]=="ALB")|(control15["genes"]=="C3")|(control15["genes"]=="APOB")|(control15["genes"]=="A2M")|(control15["genes"]=="TF")|(control15["genes"]=="SERPINA1")|(control15["genes"]=="HP")]["# PSMs"].tolist()
# list_control15=statistics.mean(list_control15)
#
# list_control16=control16[(control16["genes"]=="ALB")|(control16["genes"]=="C3")|(control16["genes"]=="APOB")|(control16["genes"]=="A2M")|(control16["genes"]=="TF")|(control16["genes"]=="SERPINA1")|(control16["genes"]=="HP")]["# PSMs"].tolist()
# list_control16=statistics.mean(list_control16)


# list_control10=control10[(control10["genes"]=="C3")|(control10["genes"]=="A2M")|(control10["genes"]=="C4B")|(control10["genes"]=="CP")|(control10["genes"]=="PLG")|(control10["genes"]=="GC")|(control10["genes"]=="CFB")|(control10["genes"]=="ITIH4")]["# PSMs"].tolist()
# list_control10=statistics.mean(list_control10)
#
#
# list_control11=control11[(control11["genes"]=="C3")|(control11["genes"]=="A2M")|(control11["genes"]=="C4B")|(control11["genes"]=="CP")|(control11["genes"]=="PLG")|(control11["genes"]=="GC")|(control11["genes"]=="CFB")|(control11["genes"]=="ITIH4")]["# PSMs"].tolist()
# list_control11=statistics.mean(list_control11)
#
#
# list_control12=control12[(control12["genes"]=="C3")|(control12["genes"]=="A2M")|(control12["genes"]=="C4B")|(control12["genes"]=="CP")|(control12["genes"]=="PLG")|(control12["genes"]=="GC")|(control12["genes"]=="CFB")|(control12["genes"]=="ITIH4")]["# PSMs"].tolist()
# list_control12=statistics.mean(list_control12)
#
#
# list_control13=control13[(control13["genes"]=="C3")|(control13["genes"]=="A2M")|(control13["genes"]=="C4B")|(control13["genes"]=="CP")|(control13["genes"]=="PLG")|(control13["genes"]=="GC")|(control13["genes"]=="CFB")|(control13["genes"]=="ITIH4")]["# PSMs"].tolist()
# list_control13=statistics.mean(list_control13)
#
#
# list_control14=control14[(control14["genes"]=="C3")|(control14["genes"]=="A2M")|(control14["genes"]=="C4B")|(control14["genes"]=="CP")|(control14["genes"]=="PLG")|(control14["genes"]=="GC")|(control14["genes"]=="CFB")|(control14["genes"]=="ITIH4")]["# PSMs"].tolist()
# list_control14=statistics.mean(list_control14)
#
#
# list_control15=control15[(control15["genes"]=="C3")|(control15["genes"]=="A2M")|(control15["genes"]=="C4B")|(control15["genes"]=="CP")|(control15["genes"]=="PLG")|(control15["genes"]=="GC")|(control15["genes"]=="CFB")|(control15["genes"]=="ITIH4")]["# PSMs"].tolist()
# list_control15=statistics.mean(list_control15)
#
#
# list_control16=control16[(control16["genes"]=="C3")|(control16["genes"]=="A2M")|(control16["genes"]=="C4B")|(control16["genes"]=="CP")|(control16["genes"]=="PLG")|(control16["genes"]=="GC")|(control16["genes"]=="CFB")|(control16["genes"]=="ITIH4")]["# PSMs"].tolist()
# list_control16=statistics.mean(list_control16)


# ONLY ALB based scaling

# list_control2=control2[control2["genes"]=="ALB"]["# PSMs"].tolist()[0]
# list_control3=control3[control3["genes"]=="ALB"]["# PSMs"].tolist()[0]
# list_control4=control4[control4["genes"]=="ALB"]["# PSMs"].tolist()[0]
# list_control5=control5[control5["genes"]=="ALB"]["# PSMs"].tolist()[0]
# list_control6=control6[control6["genes"]=="ALB"]["# PSMs"].tolist()[0]
# list_control7=control7[control7["genes"]=="ALB"]["# PSMs"].tolist()[0]
# list_control8=control8[control8["genes"]=="ALB"]["# PSMs"].tolist()[0]
# list_control9=control9[control9["genes"]=="ALB"]["# PSMs"].tolist()[0]
# list_control10=control10[control10["genes"]=="ALB"]["# PSMs"].tolist()[0]
# # list_control11=control11[control11["genes"]=="ALB"]["# PSMs"].tolist()[0]
# # list_control12=control12[control12["genes"]=="ALB"]["# PSMs"].tolist()[0]
# # list_control13=control13[control13["genes"]=="ALB"]["# PSMs"].tolist()[0]
# # list_control14=control14[control14["genes"]=="ALB"]["# PSMs"].tolist()[0]
# # list_control15=control15[control15["genes"]=="ALB"]["# PSMs"].tolist()[0]
# # list_control16=control16[control16["genes"]=="ALB"]["# PSMs"].tolist()[0]
#
#
# # Let me median scale my data!
#
# list_control1=np.median(np.array(control1['# PSMs'].tolist()))
# list_control2=np.median(np.array(control2['# PSMs'].tolist()))
# list_control3=np.median(np.array(control3['# PSMs'].tolist()))
# list_control4=np.median(np.array(control4['# PSMs'].tolist()))
# list_control5=np.median(np.array(control5['# PSMs'].tolist()))
# list_control6=np.median(np.array(control6['# PSMs'].tolist()))
# list_control7=np.median(np.array(control7['# PSMs'].tolist()))
# list_control8=np.median(np.array(control8['# PSMs'].tolist()))
# list_control9=np.median(np.array(control9['# PSMs'].tolist()))
# list_control10=np.median(np.array(control10['# PSMs'].tolist()))
# list_control11=np.median(np.array(control11['# PSMs'].tolist()))
# list_control12=np.median(np.array(control12['# PSMs'].tolist()))
# list_control13=np.median(np.array(control13['# PSMs'].tolist()))
# # list_control14=np.median(np.array(control14['# PSMs'].tolist()))
# # list_control15=np.median(np.array(control15['# PSMs'].tolist()))
# # list_control16=np.median(np.array(control16['# PSMs'].tolist()))
# # # list_control17=np.median(np.array(control17['# PSMs'].tolist()))
# # # list_control18=np.median(np.array(control18['# PSMs'].tolist()))
# #
# # Here, let me divide the PSMs of my data to median for the scaling!
#
control1["# PSMs"]=control1["# PSMs"].div(list_control1)
control2["# PSMs"]=control2["# PSMs"].div(list_control2)
control3["# PSMs"]=control3["# PSMs"].div(list_control3)
control4["# PSMs"]=control4["# PSMs"].div(list_control4)
control5["# PSMs"]=control5["# PSMs"].div(list_control5)
# control6["# PSMs"]=control6["# PSMs"].div(list_control6)
# control7["# PSMs"]=control7["# PSMs"].div(list_control7)
# control8["# PSMs"]=control8["# PSMs"].div(list_control8)
# control9["# PSMs"]=control9["# PSMs"].div(list_control9)
# control10["# PSMs"]=control10["# PSMs"].div(list_control10)
# control11["# PSMs"]=control11["# PSMs"].div(list_control11)
# control12["# PSMs"]=control12["# PSMs"].div(list_control12)
# control13["# PSMs"]=control13["# PSMs"].div(list_control13)
# control14["# PSMs"]=control14["# PSMs"].div(list_control14)
# control15["# PSMs"]=control15["# PSMs"].div(list_control15)
# control16["# PSMs"]=control16["# PSMs"].div(list_control16)
# # control17["# PSMs"]=control17["# PSMs"].div(list_control17)
# # # control18["# PSMs"]=control18["# PSMs"].div(list_control18)
# #
# #
# Let me merge the data!

z= pd.merge(control1,control2, on="Accession", how="outer")
z= pd.merge(z,control3, on="Accession", how="outer")
z= pd.merge(z,control4, on="Accession", how="outer")
z= pd.merge(z,control5, on="Accession", how="outer")
# z= pd.merge(z,control6, on="Accession", how="outer")
# z= pd.merge(z,control7, on="Accession", how="outer")
# z= pd.merge(z,control8, on="Accession", how="outer")
# z= pd.merge(z,control9, on="Accession", how="outer")
# z= pd.merge(z,control10, on="Accession", how="outer")
# z= pd.merge(z,control11, on="Accession", how="outer")
# z= pd.merge(z,control12, on="Accession", how="outer")
# z= pd.merge(z,control13, on="Accession", how="outer")
# z= pd.merge(z,control14, on="Accession", how="outer")
# z= pd.merge(z,control15, on="Accession", how="outer")
# z= pd.merge(z,control16, on="Accession", how="outer")
# z= pd.merge(z,control17, on="Accession", how="outer")
# # # z= pd.merge(z,control18, on="Accession", how="outer")
# #
new_frame= z[z.columns[[3,9,26,77,17,34,85]]]
print(new_frame)


# #
# # # # # CRITICAL PATIENTS
# # new_frame.columns=["Accession","C1_1","C1_3","C2_2","C3_1","C4_1","C5_1","C6_1","C8_1","C1_2","C2_1","C2_3","C3_2",
# #                    "C4_2","C5_2","C6_2","C8_2","genes1_1","genes1_3","genes2_2","genes3_1","genes4_1","genes5_1","genes6_1","genes8_1",
# #                    "genes1_2","genes2_1","genes2_3","genes3_2",
# #                    "genes4_2","genes5_2","genes6_2","genes8_2"]
# # #
#HEALTHY SUBJECTS
new_frame.columns=["Accession","H1","H3","H2","H4","H5","genes1","genes3"
                   ,"genes2","genes4","genes5"]

# #  Let me fill na values in genes with their corresponding names for gene1_1!!

# HEALTHY ONES
new_frame["genes1"].fillna(new_frame["genes2"],inplace=True)
new_frame["genes1"].fillna(new_frame["genes3"],inplace=True)
new_frame["genes1"].fillna(new_frame["genes4"],inplace=True)
new_frame["genes1"].fillna(new_frame["genes5"],inplace=True)
# new_frame["genes1"].fillna(new_frame["genes6"],inplace=True)
# new_frame["genes1"].fillna(new_frame["genes7"],inplace=True)
# new_frame["genes1"].fillna(new_frame["genes8"],inplace=True)
# # # new_frame["genes1"].fillna(new_frame["genes9"],inplace=True)
# # # new_frame["genes1"].fillna(new_frame["genes10"],inplace=True)


# # # # # CRITICAL ONES
# # new_frame["genes1_1"].fillna(new_frame["genes1_2"],inplace=True)
# new_frame["genes1_1"].fillna(new_frame["genes1_3"],inplace=True)
# new_frame["genes1_1"].fillna(new_frame["genes2_1"],inplace=True)
# new_frame["genes1_1"].fillna(new_frame["genes2_2"],inplace=True)
# new_frame["genes1_1"].fillna(new_frame["genes2_3"],inplace=True)
# new_frame["genes1_1"].fillna(new_frame["genes3_1"],inplace=True)
# new_frame["genes1_1"].fillna(new_frame["genes3_2"],inplace=True)
# new_frame["genes1_1"].fillna(new_frame["genes4_1"],inplace=True)
# new_frame["genes1_1"].fillna(new_frame["genes4_2"],inplace=True)
# new_frame["genes1_1"].fillna(new_frame["genes5_1"],inplace=True)
# new_frame["genes1_1"].fillna(new_frame["genes5_2"],inplace=True)
# new_frame["genes1_1"].fillna(new_frame["genes6_1"],inplace=True)
# new_frame["genes1_1"].fillna(new_frame["genes6_2"],inplace=True)
# # new_frame["genes1_1"].fillna(new_frame["genes7_1"],inplace=True)
# # new_frame["genes1_1"].fillna(new_frame["genes7_2"],inplace=True)
# new_frame["genes1_1"].fillna(new_frame["genes8_1"],inplace=True)
# new_frame["genes1_1"].fillna(new_frame["genes8_2"],inplace=True)
# # # #
# # # #
# new_frame=new_frame[["Accession","genes1_1","C1_1","C1_2","C1_3","C2_1","C2_2","C2_3","C3_1","C3_2","C4_1","C4_2","C5_1","C5_2",
#                    "C6_1","C6_2","C8_1","C8_2"]]
# #
new_frame=new_frame[["Accession","genes1","H1","H2","H3","H4","H5"]]
new_frame.columns=["Accession","genes_2","H1","H2","H3","H4","H5"]
# #

# # # #
# # myPSM_filter_df= new_frame.dropna(thresh=len(new_frame.T)*(0.1428))

# print(myPSM_filter_df["genes1_1"])

# myPSM_filter_df.fillna(0.001,inplace=True)
# myPSM_filter_df.reset_index(inplace=True)
#
# combine_filtered=pd.DataFrame({"Accession":myPSM_filter_df["Accession"],"genes":myPSM_filter_df["gene1"],"first":myPSM_filter_df["PSM_1"],
#                                "second":myPSM_filter_df["PSM_2"],"third":myPSM_filter_df["PSM_3"],"fourth":myPSM_filter_df["PSM_4"]
#                                })

# myPSM_filter_df=myPSM_filter_df[myPSM_filter_df.columns[1:]]
new_frame.to_excel("/home/ali/Desktop/combined_filtered_healthy.xlsx")
# # # "fifth":myPSM_filter_df["PSM_5"],"sixth":myPSM_filter_df["PSM_6"],"seventh":myPSM_filter_df["PSM_7"],"eighth":myPSM_filter_df["PSM_8"]
# # # # combine_filtered.to_excel("/home/a/Desktop/combined_filtered.xlsx")
