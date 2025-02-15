import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Читаем файлы
data = pd.read_csv("../result/result.txt",sep=' ', header=None)
data_stat = pd.read_csv("../result/stat_points.txt", sep=' ', header=None, skiprows=1)
data_value = pd.read_csv("../result/stat_points.txt", sep=' ', nrows=1, header=None)

# Проверяем данные на NaN и некорректные значения
data.replace(['-nan(ind)', 'nan', np.inf, -np.inf], np.nan, inplace=True)
data.dropna(inplace=True)
data = data.astype(float)  # Приводим к числовому формату

data_stat.replace(['-nan(ind)', 'nan', np.inf, -np.inf], np.nan, inplace=True)
data_stat.dropna(inplace=True)
data_stat = data_stat.astype(float)

# Проверяем размерности перед построением
print("Data shape:", data.shape)
print("Data_stat shape:", data_stat.shape)

# Создаем фигуру
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Извлекаем параметры
a = data_value.iloc[0, 0]
b = data_value.iloc[0, 1]
ax.set_title(f'a = {a}, b = {b}')

# Строим графики
ax.scatter(data_stat[0], data_stat[1], data_stat[2], color='red', label="Static Points")

# Проверяем, достаточно ли данных для построения траектории
if len(data.columns) > 3:
    ax.plot(data[1], data[2], data[3], label="Trajectory", color='green')

ax.legend()
plt.show()
