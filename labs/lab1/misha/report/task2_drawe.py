import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('../res_task2/full_trajectory.csv', sep=' ')

# plt.scatter([0],[0])

# plt.plot(data[0],data[1],label='Планета')
plt.plot(data.x3,data.y3,label='Спутник')
plt.legend()
plt.show()
