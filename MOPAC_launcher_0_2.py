#!/usr/bin/env python
# coding: utf-8

# In[1]:


import subprocess
import os
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter.filedialog import askopenfilename
import shutil


# In[2]:


def getFolderPath():
    root = Tk()
    root.withdraw()
    folder_selected = filedialog.askdirectory()
    #folderPath.set(folder_selected)
    root.destroy()
    return (folder_selected)


# In[3]:


def find_mopac_exe():
    way = input("Specify your MOPAC.exe way or search directory: ")
    mop_dir = re.search('MOPAC.*\.exe', way)
    if mop_dir:
        return way
    else:
        find_files = []
        #find_exe = []
        for root, dirs, files in os.walk(way):
            find_files += [os.path.join(root,name) for name in files]
        for d in range(len(find_files)):
            mop_dir = re.search('MOPAC.*\.exe', find_files[d])
            if mop_dir:
                return find_files[d]


# In[4]:


mopac_way = find_mopac_exe()
print(mopac_way)


# In[5]:


file_vec = []
file_vec_short = []
dirs = []
#for root, dirs, files in os.walk(os.getcwd()):
for root, dirs, files in os.walk(getFolderPath()):
    for file in files:
        if file.endswith(".mop"):
            current_way = os.path.join(root,file)
            current_file = os.path.join(file)
            file_vec.append(current_way)
            file_vec_short.append(current_file)


# In[6]:


for f in range(len(file_vec)):
    current_dir = file_vec[f].replace(file_vec_short[f],'')
    dirs.append(current_dir)


# In[7]:


#print(file_vec)


# In[8]:


for i in range(len(file_vec)):
    cur_folder = dirs[i]+file_vec_short[i][:-4]
    if not os.path.exists(cur_folder):
        os.makedirs(cur_folder)
    way = cur_folder+'/'+file_vec_short[i]
    shutil.move(file_vec[i], way) # Это лучше os.replace
    start_mopac_bash ='yes '' |' + mopac_way + ' ' + way
    try:
        print(way + ' is being calculating')
        #subprocess.call(mopac_way, way)        
        os.system(start_mopac_bash)
        print(way + ' calculated')
    except:
        print(way + " WASTED")

