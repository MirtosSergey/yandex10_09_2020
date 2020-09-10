#Подключение библиотек
import math
import pandas as pd
import PreSolver
import Solver

#Чтение конфиг-файла
config = PreSolver.read_config()
#Импорт данных с csv
C4401 = pd.read_csv('DATA/C4401.csv', delimiter = '\t')
HADC = pd.read_csv('DATA/HADC.csv', delimiter = '\t')
ADC = pd.read_csv('DATA/ADC.csv', delimiter = '\t')
CDU = pd.read_csv('DATA/CDU.csv', delimiter = '\t')
CRoc = pd.read_csv('DATA/CRoc.csv', delimiter = '\t')
#Подготовка данных
ADC.loc[:,['D', 'A']] /= 57.3
data_C4401 = Solver.get_data_from_DF(C4401, ['H'], ['P', 'ro', 'T'])
data_HADC = Solver.get_data_from_DF(HADC, ['H', 'M'], ['DC_x'])
data_ADC = Solver.get_data_from_DF(ADC, ['M', 'D', 'A'], ['C_x', 'C_y', 'C_xa', 'C_ya', 'm_z', 'C_x_d'])
data_CDU = Solver.get_data_from_DF(CDU, ['t'], ['G', 'R'])
data_CRoc = Solver.get_data_from_DF(CRoc, ['mu'], ['x_cm'])
#Функции вызова характеристик
def get_C4401(H, Char):
    #H
    if H < data_C4401[1]['H'][0]:
        H = data_C4401[1]['H'][0]
    elif H > data_C4401[1]['H'][-1]:
        H = data_C4401[1]['H'][-1]
    return Solver.N_linterp(data_C4401[0], [H], Char)
def get_HADC(H, M):
    #H
    if H < data_HADC[1]['H'][0]:
        H = data_HADC[1]['H'][0]
    elif H > data_HADC[1]['H'][-1]:
        H = data_HADC[1]['H'][-1]
    #M
    if M < data_HADC[1]['M'][0]:
        M = data_HADC[1]['M'][0]
    elif M > data_HADC[1]['M'][-1]:
        M = data_HADC[1]['M'][-1]
    return Solver.N_linterp(data_HADC[0], [H, M], ['DC_x'])
def get_ADC(M, D, A, Char):
    #M
    if M < data_ADC[1]['M'][0]:
        M = data_ADC[1]['M'][0]
    elif M > data_ADC[1]['M'][-1]:
        M = data_ADC[1]['M'][-1]
    #D
    if D < data_ADC[1]['D'][0]:
        D = data_ADC[1]['D'][0]
    elif D > data_ADC[1]['D'][-1]:
        D = data_ADC[1]['D'][-1]
    #A
    if A < data_ADC[1]['A'][0]:
        A = data_ADC[1]['A'][0]
    elif A > data_ADC[1]['A'][-1]:
        A = data_ADC[1]['A'][-1]
    return Solver.N_linterp(data_ADC[0], [M, D, A], Char)
def get_CDU(t, Char):
    #t
    if t < data_CDU[1]['t'][0]:
        t = data_CDU[1]['t'][0]
    elif t > data_CDU[1]['t'][-1]:
        t = data_CDU[1]['t'][-1]
    return Solver.N_linterp(data_CDU[0], [t], Char)
def get_CRoc(mu, Char):
    #mu
    if mu < data_CRoc[1]['mu'][0]:
        mu = data_CRoc[1]['mu'][0]
    elif mu > data_CRoc[1]['mu'][-1]:
        mu = data_CRoc[1]['mu'][-1]
    return Solver.N_linterp(data_CRoc[0], [mu], Char)

#Скорость звука
def func_a(k, R, T):
    return math.sqrt(k * R * T)

#Коэффициент момента тангажа с учетом изменения положения центра масс
def func_m_z(C_y, m_z, L_ref, x_cm_ref, mu):
    return C_y / L_ref * (get_CRoc(mu, ['x_cm']) - x_cm_ref) + m_z

#Оределения радиуса кривизны траектории
def R_path():
    return 1e1000

#Определение потребных перегрузок
def func_n_ya_potr(V, m, X, teta, R_du, alpha, g):
    return V ** 2 / (g * R_path()) + math.cos(teta) - R_du * math.sin(alpha) / (m * g)

#Функция для определение потребного угла атаки
def func_alpha_potr(V, m, X, teta, R_du, Mach, delta, alpha, q, S_ref, g):
    n_ya_potr = func_n_ya_potr(V, m, X, teta, R_du, alpha, g)
    C_ya = get_ADC(Mach, delta, alpha, ['C_ya'])
    return C_ya * q * S_ref - m * g * n_ya_potr

#Функция для определение потребного угла отклонения рулей
def func_delta_potr(Mach, delta, alpha, L_ref, x_cm_ref, _m):
    C_y, m_z = get_ADC(Mach, delta, alpha, ['C_y', 'm_z'])
    return func_m_z(C_y, m_z, L_ref, x_cm_ref, _m)

#Определение потребного режима работы
def Find_mode(eps, n_max, n_iter_max, V, m, X, teta, R_du, Mach, delta, alpha, q, S_ref, L_ref, g, x_cm_ref, mu, delta_max, delta_min, alpha_max, alpha_min):
    alpha_i, delta_i, i = 1, 1, 0
    while (abs(delta - delta_i) > eps or abs(alpha - alpha_i) > eps):
        delta_i, alpha_i = delta, alpha
        #Поиск варьируемых параметров
        def func_alpha_potr_pub(x):
            return func_alpha_potr(V, m, X, teta, R_du, Mach, delta, x, q, S_ref, g)
        alpha = Solver.Half_func(func_alpha_potr_pub, 0, alpha_min, alpha_max, eps, n_max)
        def func_delta_potr_pub(x):
            return func_delta_potr(Mach, x, alpha, L_ref, x_cm_ref, mu)
        delta = Solver.Half_func(func_delta_potr_pub, 0, delta_min, delta_max, eps, n_max)
        i += 1
        if (i > n_iter_max):
            break
    return delta, alpha

#Система дифференциальных уравнений внешней баллистики
def Sys_out_path(Matr):
    Mas = Matr[-1,:]
    #Исходные данные
    t, V, m, X, Y, teta = Mas[:6]
    eps = config['settings']['eps']
    n_max = config['settings']['n half max']
    n_iter_max = config['settings']['n iter max']
    k = config['physical const']['k']
    R =  config['physical const']['R']
    g = config['physical const']['g']
    S_ref = config['aerodynamic characteristics']['S ref']
    L_ref = config['aerodynamic characteristics']['L ref']
    H_ref = config['aerodynamic characteristics']['H ref']
    x_cm_ref = config['aerodynamic characteristics']['x cm ref']
    delta_max = config['aerodynamic characteristics']['delta max']
    delta_min = config['aerodynamic characteristics']['delta min']
    m_0 = config['rocket']['m 0']
    m_k = config['rocket']['m k']
    #Расчет
    delta_max, delta_min = math.radians(delta_max), math.radians(delta_min)
    A_max_all, A_min_all = data_ADC[1]['A'][-1], data_ADC[1]['A'][0]
    if len(Matr) == 1:
        delta, alpha = 0, 0
    else:
        delta, alpha = Matr[-2:-1,9], Matr[-2:-1,10]
    T, ro = get_C4401(Y, ['T', 'ro'])
    G_du = R_du = 0
    a = func_a(k, R, T)
    Mach = V / a
    DC_x = get_HADC(Y, Mach) - get_HADC(H_ref, Mach)
    q = 0.5 * ro * V ** 2
    mu = (m_0 - m) / (m_0 - m_k)
    if mu <= 1:
        G_du, R_du = get_CDU(t, ['G', 'R'])
    def func_alpha_max(x):
        return get_ADC(Mach, delta_min, x, ['m_z'])
    alpha_max = Solver.Half_func(func_alpha_max, 0, A_min_all, A_max_all, eps, n_max)
    def func_alpha_min(x):
        return get_ADC(Mach, delta_max, x, ['m_z'])
    alpha_min = Solver.Half_func(func_alpha_min, 0, A_min_all, A_max_all, eps, n_max)
    delta, alpha = Find_mode(eps, n_max, n_iter_max, V, m, X, teta, R_du, Mach, delta, alpha, q, S_ref, L_ref, g, x_cm_ref, mu, delta_max, delta_min, alpha_max, alpha_min)
    C_xa, C_ya, C_y, m_z = get_ADC(Mach, delta, alpha, ['C_xa', 'C_ya', 'C_y', 'm_z'])
    C_xa = C_xa + DC_x * math.cos(alpha)
    AD_K = C_ya / C_xa
    m_z_new = func_m_z(C_y, m_z, L_ref, x_cm_ref, mu)
    n_xa = (R_du * math.cos(alpha) - C_xa * q * S_ref) / (m * g)
    n_ya = (R_du * math.sin(alpha) + C_ya * q * S_ref) / (m * g)
    n_ya_potr = func_n_ya_potr(V, m, X, teta, R_du, alpha, g)
    n_ya_rasp_max = (R_du * math.sin(alpha_max) + get_ADC(Mach, delta_min, alpha_max, ['C_ya']) * q * S_ref) / (m * g)
    n_ya_rasp_min = (R_du * math.sin(alpha_min) + get_ADC(Mach, delta_max, alpha_min, ['C_ya']) * q * S_ref) / (m * g)
    #Основные параметры
    dV = (n_xa - math.sin(teta)) * g
    dm = -G_du
    dX = V * math.cos(teta)
    dY = V * math.sin(teta)
    dteta = (n_ya - math.cos(teta)) * g / V #dwz = M_z / J_z#dalpha = w_z - dteta
    #Вывод
    i = 6
    Mas[i] = R_du; i += 1
    Mas[i] = G_du; i += 1
    Mas[i] = Mach; i += 1
    Mas[i] = delta; i += 1
    Mas[i] = alpha; i += 1
    Mas[i] = C_xa; i += 1
    Mas[i] = C_ya; i += 1
    Mas[i] = AD_K; i += 1
    Mas[i] = m_z_new; i += 1
    Mas[i] = alpha_max; i += 1
    Mas[i] = alpha_min; i += 1
    Mas[i] = mu; i += 1
    Mas[i] = n_xa; i += 1
    Mas[i] = n_ya; i += 1
    Mas[i] = n_ya_potr; i += 1
    Mas[i] = n_ya_rasp_max; i += 1
    Mas[i] = n_ya_rasp_min; i += 1
    #Ограничитель
    key = 1
    if (Y < 0):
        key = 0
    return [1, dV, dm, dX, dY, dteta], key

def Calc_one_out_path():
    init = [0] * 23
    init[0] = config['rocket']['t 0']
    init[1] = config['rocket']['V 0']
    init[2] = config['rocket']['m 0']
    init[3] = config['rocket']['X 0']
    init[4] = config['rocket']['Y 0']
    init[5] = config['rocket']['teta 0']
    res = Solver.Metod_Eiler(init, Sys_out_path, config['settings']['n eiler max'], config['settings']['dt'])
    return res
    