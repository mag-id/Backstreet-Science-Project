#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import re
import pandas as pd
import datetime
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter.filedialog import askopenfilename


# In[2]:


def getFolderPath():
    root = Tk()
    root.withdraw()
    folder_selected = filedialog.askdirectory()
    #folderPath.set(folder_selected)
    root.destroy()
    return (folder_selected)


# In[3]:


cur_dir = getFolderPath()


# In[4]:


file_vec = []
file_vec_short = []
#for root, dirs, files in os.walk(os.getcwd()):
for root, dirs, files in os.walk(cur_dir):
    for file in files:
        if file.endswith(".out"):
            current_way = os.path.join(root,file)
            current_file = os.path.join(file)
            file_vec.append(current_way)
            file_vec_short.append(current_file)


# In[5]:


#values_kcal_therm = []
#values_kJ_therm = []
values_kcal_SP = []
values_kJ_SP = []
therm_files = []
SP_files = []
key_words = ('VIB', 'TOT')
HOF = []
ENT = []
HC = []
ENR = []
T = []
for file in file_vec:
    with open(str(file), "r") as f:
        switch = 0
        lines = f.readlines()
        for line in lines:
            therm_prop = re.search('\s([A-Z]{3})\.\s', line)
            if 'FINAL HEAT OF FORMATION' in line:
                values_kcal = []
                values_kJ = []
                line_vec = line.split('=')
                value_1 = re.search('-?(\d+\.\d+)',line_vec[1])
                value_2 = re.search('-?(\d+\.\d+)',line_vec[2])
                values_kcal.append(float(value_1.group(0)))
                values_kJ.append(float(value_2.group(0)))
            if therm_prop and therm_prop.group(1) in key_words:
                switch = 1
                therm_val = line.split()
                if therm_prop.group(1) == 'VIB':
                    T.append(float(therm_val[0]))
                else:
                    HOF.append(float(therm_val[1]))
                    ENT.append(float(therm_val[2]))
                    HC.append(float(therm_val[3]))
                    ENR.append(float(therm_val[4]))
                    
    if switch == 1:
        #values_kcal_therm.extend(values_kcal)
        #values_kJ_therm.extend(values_kJ)
        for name in file_vec_short:
            if name in file:
                therm_files.append(name)
    else:
        values_kcal_SP.extend(values_kcal)
        values_kJ_SP.extend(values_kJ)
        for name in file_vec_short:
            if 'sp_'+ name in file:
                SP_files.append('sp_'+name)
    #print(file)
    #print(switch)
    #print(therm_files)
    #print(SP_files)


# In[48]:


m = len(T)/len(therm_files)
therm_files_mul = []
for i in range (len(therm_files)):
    for j in range(int(m)):
        therm_files_mul.append(therm_files[i])
tuples = list(zip(therm_files_mul,T))


# In[57]:


tuples = list(zip(therm_files_mul,T))
index = pd.MultiIndex.from_tuples(tuples, names = ['file','T,K'])
data_therm = pd.DataFrame({"H.O.F,KCAL/MOL":pd.Series(HOF, index = index),"ENTHALPY,CAL/MOLE":pd.Series(ENT, index = index),
                          "HEAT CAPACITY,CAL/K/MOL":pd.Series(HC, index = index), "ENTROPY,CAL/K/MOL":pd.Series(ENR, index = index)})
data_therm = data_therm.sort_values(['T,K','file'])


# In[58]:


print(data_therm)


# In[29]:


#data_therm_FHF = pd.DataFrame({"kcal/mol":pd.Series(values_kcal_therm,index = therm_files), "kJ/mol":pd.Series(values_kJ_therm,index = therm_files)})
data_SP_FHF = pd.DataFrame({"kcal/mol":pd.Series(values_kcal_SP, index = SP_files), "kJ/mol":pd.Series(values_kJ_SP, index = SP_files)})
data_SP_FHF.index.names = ['file']
data_SP_FHF = data_SP_FHF.sort_values(['file'])


# In[30]:


now = datetime.datetime.now()
newpath = str(cur_dir)+'/MOPAC_parsed/'
if not os.path.exists(newpath):
    os.makedirs(newpath)
data_SP_FHF.to_csv(str(newpath)+str(now.strftime("%d-%m-%Y_%H-%M"))+'_SP_FHF.csv')
#data_therm_FHF.to_csv(str(newpath)+str(now.strftime("%d-%m-%Y_%H-%M"))+'_therm_FHF.csv')
data_therm.to_csv(str(newpath)+str(now.strftime("%d-%m-%Y_%H-%M"))+'_therm.csv')


# In[15]:


print(data_SP_FHF)

