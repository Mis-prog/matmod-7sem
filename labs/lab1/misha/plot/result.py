import matplotlib.pyplot as plt
import pandas as pd

# Чтение данных
data = pd.read_csv('output.csv', sep=' ')

end = 36500000000

plt.plot(data.x2[:end],data.y2[:end],label='Марс')
plt.plot(data.x3[:end],data.y3[:end],ls='--',label='Фобос')

plt.legend()
plt.show()