import matplotlib.pyplot as plt
import pandas as pd

# Чтение данных
data = pd.read_csv('output.txt', sep=' ')

# Определяем количество данных для одного участка
end = 1660

# Убираем лишние точки из цикла
for i in range(4):
    start = i * end
    end = (i + 1) * end
    plt.plot(data.x2[start:end], data.y2[start:end], label=f"Марс, участок {i + 1}")

plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()