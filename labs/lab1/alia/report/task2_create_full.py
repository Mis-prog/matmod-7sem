import re
import matplotlib.pyplot as plt
import numpy as np

# Путь к файлу
file_path = "full_fuel.txt"

# Регулярное выражение
pattern = r"mt\s*=\s*([\d.]+)\s*для\s*угла\s*(\d+)"


MT=[]
ANGLE=[]
# Открываем файл и считываем построчно
with open(file_path, "r", encoding="utf-8") as file:
    for line in file:
        # Ищем совпадения в текущей строке
        match = re.search(pattern, line)
        if match:
            mt, angle = match.groups()
            MT.append(float(mt))
            ANGLE.append(float(angle))

MT=np.array(MT)
i=MT.argmin()

plt.scatter(ANGLE[i],MT[i],label=f'Минимальное кол-во топлива {MT[i]} для угла {ANGLE[i]}')
plt.plot(ANGLE, MT)
plt.xlabel("Угол")
plt.ylabel("Кол-во топлива")
plt.grid()
plt.legend()
plt.show()