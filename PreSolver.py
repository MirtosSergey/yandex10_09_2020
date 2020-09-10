#Подключение библиотек
import pandas as pd
import json

#Импорт config.json файла
def read_config():
    with open('config.json', 'r') as f:
        return json.load(f)
    
#Перезапись файла
def rewrite_res():
    res_old = pd.read_csv('RESULT/result.csv', delimiter = '\t')
    res_old.to_csv('RESULT/result_old.csv', sep = '\t', index = False)