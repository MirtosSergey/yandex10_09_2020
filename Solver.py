#Подключение библиотек
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

#Метод Эйлера для np.array
def Metod_Eiler(init, sd, n, dt):
    Mas, percent = np.array([init]), 0
    for i in range(n):
        dX, key = sd(Mas)
        if key == 1:
            Mas = np.append(Mas, Mas[-1:,:] + np.array(dX + [0] * (Mas.shape[1] - len(dX))) * dt, axis = 0)
        else:
            break
        if i * 100 / n > percent:
            percent += 1
            print("Calculation", percent, "%")
    return Mas

#Метод половинного деления
def Half_func(f, res, a, b, eps, max_iter):
    i, c = 0, 0.5 * (a + b)
    C, s = f(c), np.sign(f(0.75 * a + 0.25 * b) - f(0.25 * a + 0.75 * b))
    while (abs(C - res) > eps and i < max_iter):
        if (s * C - res >= 0):
            a = c
        else:
            b = c
        c = 0.5 * (a + b)
        C = f(c)
        i += 1
    return c

#Линейная интерполяция
def linterp(x, a, b, A, B):
    if (x < a):
        x = a
    elif (x > b):
        x = b
    return A + (B - A) * (x - a) / (b - a)

#Подготовка данных DF для N-мерной интерполяции
def get_data_from_DF(DF, Mx, My):
    Mas_x = {}
    val_shape = []
    for x in Mx:
        Mas_x[x] = np.array(DF.loc[:,x].drop_duplicates()).astype(np.float64)
        val_shape = np.append(val_shape, DF.loc[:,x].nunique())
    Mas_y, res = {}, {}
    for y in My:
        Mas_y[y] = np.array(DF.loc[:,y]).reshape(val_shape.astype(np.int64))
        res[y] = RegularGridInterpolator(tuple(Mas_x.values()), Mas_y[y], method = "linear", bounds_error=False, fill_value=None)
    return res, Mas_x

#N-мерная интерполяция
def N_linterp(M, X, Y):
    res = np.array([])
    for y in Y:
        res = np.append(res, M[y](X))
    if len(res) == 1:
        res = res[0]
    return res

#Определение индекса в DataFrame по элементу
def Find_index(DF, val, max_i):
    if val <= DF.values[0]:
        return [DF.index[0], DF.index[1], DF.index[2]]
    elif val >= DF.tail(n=1).values[0]:
        return [DF.index[len(DF) - 2], DF.index[len(DF) - 1], DF.index[0] + max_i]
    else:
        for i in range(len(DF) - 1):
            if (DF.values[i] < val <= DF.values[i + 1]):
                if (i == len(DF) - 2):
                    return [DF.index[i], DF.index[i + 1], DF.index[0] + max_i]
                else:
                    return [DF.index[i], DF.index[i + 1], DF.index[i + 2]]
                
#N-мерная интерполяция по DataFrame с случайным шагом
def N_linterp_random(DF, point_dict, char):
    point, Mas_DF, Columns = pd.DataFrame([point_dict]), [DF], []
    for i in range(point.shape[1]):
        Pub_DF = []
        for j in range(len(Mas_DF)):
            New_DF = Mas_DF[j]
            Mas = New_DF[point.columns[i]].drop_duplicates()
            ind_1, ind_2, ind_3 = Find_index(Mas, point[point.columns[i]][0], len(New_DF))
            Pub_DF.append(New_DF.loc[np.arange(ind_1,ind_2)])
            Pub_DF.append(New_DF.loc[np.arange(ind_2,ind_3)])
        Mas_DF = Pub_DF
    for i in range(point.shape[1]):
        Columns = np.append(Columns, point.columns[i])
    Columns = np.append(Columns, char)
    for i in range(point.shape[1]):
        Pub_DF = []
        for j in range(0, len(Mas_DF), 2):
            col = point.columns[point.shape[1] - i - 1]
            New_mas = Mas_DF[j][Columns]
            New_mas[char].values[0] = linterp(point[point.columns[point.shape[1] - i - 1]][0], Mas_DF[j][col].values[0], Mas_DF[j + 1][col].values[0], Mas_DF[j][char].values[0], Mas_DF[j + 1][char].values[0])
            Pub_DF.append(New_mas)
        Mas_DF = Pub_DF
    return Mas_DF[0][char].values[0]

#Метод Эйлера для DataFrame
def Metod_Eiler_DF(n_max, dt, DF_init, func_stop, Sys_Def):
    Mas_res = pd.DataFrame([DF_init])
    for i in range(n_max):
        X, dX = Sys_Def(Mas_res)
        Mas_res = Mas_res[0:len(Mas_res)-1]
        Mas_res = Mas_res.append(X)
        X_new = X + dX * dt
        if (func_stop(X_new) == 0):
            break
        Mas_res = Mas_res.append(X_new)
    return Mas_res