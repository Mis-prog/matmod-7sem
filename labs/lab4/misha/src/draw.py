import matplotlib.pyplot as plt
import  pandas as pd

data=pd.read_csv("../result/result.txt",sep=' ',header=None)
data_stat=pd.read_csv("../result/stat_points.txt",sep=' ',header=None,skiprows=1)
data_value = pd.read_csv("../result/stat_points.txt", sep=' ', nrows=1,header=None)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

a = data_value.iloc[0, 0] 
b = data_value.iloc[0, 1]

ax.set_title(f'a = {a}, b = {b}') 
ax.scatter(data_stat[0], data_stat[1], data_stat[2],color='red')
ax.plot(data[1], data[2], data[3], label="Trajectory", color='blue')