import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('../res_task2/init_data.txt', sep=' ', header=None)

print(data)

# plt.scatter([0],[0])

plt.plot(data[0],data[1],label='Планета')
plt.plot(data[2],data[3],label='Звезда')
plt.legend()
plt.show()
