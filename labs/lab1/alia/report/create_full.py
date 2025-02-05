import pandas as pd
import matplotlib.pyplot as plt

# Читаем данные из файла
df = pd.read_csv("fuel_optimal_90.txt", sep="\s+", header=None, names=["Index", "Value"])

# Строим график
plt.plot(df["Index"], df["Value"], linestyle="-")
plt.xlabel("Угол")
plt.ylabel("Масса")
# plt.title("График данных из файла")
plt.grid()
plt.show()
