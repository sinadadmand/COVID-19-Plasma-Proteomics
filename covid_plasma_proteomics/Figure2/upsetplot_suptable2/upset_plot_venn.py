import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import generate_counts
example = generate_counts()






# print(example)
# from upsetplot import plot
# plot(example)
# from matplotlib import pyplot
# pyplot.show()

example.to_excel("/home/ali/Desktop/bibak.xlsx")

path="/home/ali/Desktop/results_covid_18_1_22_new_healthy/Figure2/Figure 2E/SupplementaryTable1.xlsx"
data= pd.read_excel(path)
print(data)


a = data[(data["Geyer at al"]=="-")&(data["Shu at al"]=="-")&(data["Shen at al"]=="-")
                    &(data["Park at al"]=="-")&(data["Demichev at al"]=="-")
                    &(data["Messner at al"]=="-")&(data["In this study"]=="+")]

a.to_excel("/home/ali/Desktop/sonuc.xlsx")

name=data.columns[2:]

# C1
symbols=[]
plus_minus=["+","-"]
counts=[]

for first in plus_minus:
    for second in plus_minus:
        for third in plus_minus:
            for fourth in plus_minus:
                for five in plus_minus:
                    for sixth in plus_minus:
                        for seventh in plus_minus:
                            number = len(data[(data["Geyer at al"] == first) & (data["Shu at al"] == second) & (data["Shen at al"] == third) & (
                                        data["Park at al"] == fourth)
                                              & (data["Demichev at al"] ==five) & (data["Messner at al"] == sixth) & (
                                                          data["In this study"] == seventh)])
                            symbols.append([first, second, third, fourth, five,sixth,seventh])
                            counts.append(number)

print(counts)
print(symbols)

data_frame=pd.DataFrame(symbols)
print(data_frame)
print(name)

data_frame.columns=name

data_frame.columns=name
data_frame[""]=pd.Series(counts)
print(data_frame)


data_frame.replace({"+":True,"-":False},inplace=True)

print(data_frame)

data_= data_frame.values.tolist()

print(data_)


my_all=[]
for each in data_:
    for rn in range(0,each[-1]):
        my_all.append(each[:-1])
print(my_all)

data_f=pd.DataFrame(my_all)
print(data_frame)

data_f.to_excel("/home/ali/Desktop/result.xlsx")

data_f.columns=name


data_f.columns=["Geyer at al","Shu at al","Shen at al","Park at al","Demichev at al","Messner at al","This study"]

print(data_frame)
# Now let me re-arrange our raws!

print()
# print(name)
#
#

# example= data_f.groupby([name[0],name[1],name[2],name[3],name[4]]).count()

# print(data_f)

print(data_f)

s = data_f.groupby(data_f.columns.to_list()).size()

print(s)


font = {'family' : 'normal',

        'size'   : 18}

matplotlib.rc('font', **font)

from upsetplot import plot

plot(s,facecolor="gray", element_size=20)



from matplotlib import pyplot
plt.ylabel("Intersection Size",fontsize=20)
plt.yticks(fontsize=20)
# plt.legend(fontsize=20)


pyplot.show()
