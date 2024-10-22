import matplotlib.pyplot as plt
import pandas as pd

# Чтение данных
data = pd.read_csv('output.txt', sep=' ')

end=1000
plt.plot(data.x2[:end],data.y2[:end],label='Земля')
plt.plot(data.x3[:end],data.y3[:end],label='Луна')
plt.legend()
plt.show()