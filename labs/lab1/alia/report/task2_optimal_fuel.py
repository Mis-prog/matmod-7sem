import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('../result/task2/full_optimall_fuel_accuracity.txt', sep='|', header=None)
data[1] = data[1].str.replace(r'.*\s(\d+\.\d+)', r'\1', regex=True)
data[1] = data[1].str.replace('°', '', regex=False)
data[1] = pd.to_numeric(data[1])

data[2] = data[2].str.replace(r'.*\s(\d+\.\d+)', r'\1', regex=True)
data[2] = pd.to_numeric(data[2])

data_sorted = data.sort_values(by=1).reset_index(drop=True)

min_fuel_index = data_sorted[2].idxmin()
min_fuel_row = data_sorted.iloc[min_fuel_index]
# print(min_fuel_row)

# print(data_sorted.iloc[20:30])

plt.title('Зависимость топлива и угла')
plt.plot(data_sorted[1], data_sorted[2])
plt.scatter(min_fuel_row[1], min_fuel_row[2], label=f'Угол {min_fuel_row[1]}, топливо {min_fuel_row[2]}')
plt.xlabel('Угол')
plt.ylabel('Топливо')
plt.grid()
plt.legend()
plt.show()
