#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import re
import scm.plams
from scm import plams
from scm.plams import *


# In[35]:


init() # Инициализация протоколов PLAMS


# In[36]:


atoms = [('H','0.942179','0.000000', '-0.017370'), ('H','-0.471089','0.815951', '-0.017370'), 
         ('N','0.000000','0.000000', '0.383210'), ('H','-0.471089','-0.815951', '-0.017370')]


# In[37]:


mol_name = 'ConformerX'


# In[38]:


'''Химико-PLAMSовский словарь'''
atom_names_to_num = {'H':1,'C':6,'N':7,'O':8,'F':9,'P':15,'S':16,'Cl':17,'As':33,'Se':34,'Br':35,'I':53,'Ru':44}


# In[39]:


'''Собираем молекулу из ..листа и палок, собственно лист на вход, молекула класса Molecule на выход'''
def molecule_builder (atoms:list):
    mol = Molecule()
    for i in range(len(atoms)):
        mol.add_atom(Atom(atnum = atom_names_to_num[atoms[i][0]], coords = (atoms[i][1], atoms[i][2],atoms[i][3])))
    mol.guess_bonds() # Достраиваем связи по координатм атомов
    print(mol)
    return (mol)


# In[40]:


'''Ищем-рыщем exe файл MOPAC, путь на входъ, путь на выходъ'''
def find_mopac_exe():
    way = input("Specify your MOPAC.exe way or search directory: ") # Указываем путь до exe файла MOPAC на машине,
                                                                    # Да, можно не точно)
    mop_dir = re.search('MOPAC.*\.exe', way)
    if mop_dir:
        return way
    else:
        find_files = []
        for root, dirs, files in os.walk(way): # Все равно найдет) 
            find_files += [os.path.join(root,name) for name in files]
        for d in range(len(find_files)):
            mop_dir = re.search('MOPAC.*\.exe', find_files[d])
            if mop_dir:
                return find_files[d] # Возвращает путь до exe файла


# In[41]:


mopac_way = find_mopac_exe() # Нашли путь до exe файла MOPAC на машине


# In[42]:


'''Сборка PUT IN файла MOPAC, на вход молекула, ее погоняло, расчетный метод, мультиплетность, и чет там про COSMO'''
def MOPAC_runner (molecule, name:str, calc_method = 'AM1',charge = 0, multiplet = 'SINGLET', NSPA=92, EPS=78.4):
    trial_job = MOPACJob(molecule = molecule, name = name)
    trial_job.settings.input.LARGE = True # Пишем простыню текста
    trial_job.settings.input.CHARGE=charge # Закладываем заряд
    trial_job.settings.input.COSCCH = True # Щепотку COSMO
    trial_job.settings.input.NSPA=NSPA # Число сегментов модели COSMO
    trial_job.settings.input.EPS=EPS # Диэлектрическая проницаемость COSMO
    trial_job.settings.input.NOOPT = True # Расчет без оптимизации
    
    '''Из*бский подгониан мультиплетности'''
    
    if multiplet == "QUINTET":
        trial_job.settings.input.QUINTET = True
    elif multiplet == "DOUBLET":
        trial_job.settings.input.DOUBLET = True
    elif multiplet == "TRIPLET":
        trial_job.settings.input.TRIPLET = True
    elif multiplet == "QUARTET":
        trial_job.settings.input.QUARTET = True
    else:
        trial_job.settings.input.SINGLET = True # По умолчанию

    '''Из*бский подгониан расчетного метода'''
    
    if calc_method == 'PM3':
        trial_job.settings.input.PM3 = True
    elif calc_method == 'PM6':
        trial_job.settings.input.PM6 = True
    elif calc_method == 'PM7':
        trial_job.settings.input.PM7 = True
    else:
        trial_job.settings.input.AM1 = True # По умолчанию
        
    trial_job._command = mopac_way # Передаем путь до exe файла MOPAC, важно, мы не в Голландии!
    
    '''Пляшем с бубном на job, шаманим бинарные входные файлы для расчета'''
    trial_job.get_input() # Генерация инпут ключевый слов
    trial_job.get_runscript() # Сборка бинарного mop файла

    results = trial_job.run() # Поехали...!
    return (results)


# In[43]:


'''Парсим результаты РАБоты MOPAC, объект рабочего класса MOPACResults на вход'''
def parsim_fast_fast(data):
    data = data.grep_output(pattern = 'FINAL HEAT OF FORMATION') # Нетривиальный способ парсить out файл, 
                                                                 # выбираем нужную строку
    data = data[0].split()
    return (data[5:9:3]) # Выбор нужных цифр из строки


# In[44]:


parsim_fast_fast(MOPAC_runner(molecule_builder(atoms),name = mol_name))# Код в 1 строчку 


# In[45]:


finish() # Завершение протоколов PLAMS

'''Ищем папки plams_workdir, в которых лежат папки с названием конформера.
   Собственно, там находятся out файлы
   Многопоточность завтра навернуть попробую'''

