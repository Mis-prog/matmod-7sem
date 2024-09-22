import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('output.txt', sep=' ')

plt.plot(data.x2,data.y2)
plt.xlabel('x')
plt.ylabel('y')
plt.show()