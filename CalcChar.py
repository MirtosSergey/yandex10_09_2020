import pandas as pd

dict_init = {'ADC': ['M', 'D', 'A'], 'HADC': ['H', 'M'], 'C4401': ['H']}

def read_dict(dict_init):
    Mas_dict = []
    for name, val in dict_init.items():
        name += '.csv'
        Mas = pd.read_csv(name, delimiter = '\t', skiprows = 0)
        Mas.set_index(val, inplace = True)
        Mas_dict.append(Mas.stack())
    return Mas_dict

#print(read_dict(dict_init))