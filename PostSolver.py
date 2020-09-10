#Подключение библиотек
import pandas as pd
#import numpy as np
import matplotlib.pyplot as plt

#Запись данных
def write_res(res):
    res = pd.DataFrame(res, columns=['t', 'V', 'm', 'X', 'Y', 'teta',
                                     'R', 'G',
                                     'Mach', 'delta', 'alpha', 'C_xa', 'C_ya', 'AD_K', 'm_z', 'alpha_max', 'alpha_min',
                                     'mu',
                                     'n_xa', 'n_ya', 'n_ya_potr', 'n_ya_rasp_max', 'n_ya_rasp_min'])
    res.loc[:,['teta', 'delta', 'alpha', 'alpha_max', 'alpha_min']] *= 57.3
    res.to_csv('RESULT/result.csv', sep = '\t', index = False)
    return res
  
#Сохранение графиков
def write_graf(DF, val):
    ###Доделать!!!
    fig, ax = plt.subplots()
    P = DF.plot(x = 't', y = 'V', style = ['r-'], lw = 2)
    plt.xlim(0, DF.loc[:,'t'].max())
    plt.ylim(0, DF.loc[:,'V'].max())
    plt.xlabel('t, c', labelpad = 3, fontsize = 12, fontweight = 1000)
    plt.ylabel('V, м/с', labelpad = 3, fontsize = 12, fontweight = 1000)
    plt.title('Профиль скорости', fontstyle = 'italic', fontsize = 14)
    P.legend(['V, м/с'], loc='best')
    plt.grid(True)
    plt.savefig('GRAPH/V.png')
    return 1