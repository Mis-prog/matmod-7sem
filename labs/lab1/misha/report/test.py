import pandas as pd
import matplotlib.pyplot as plt

R_start_planeta=696340e3+228e9+3390e3

# data_planeta_0=pd.read_csv('../result/task1/path_planeta_0.csv',sep=' ')
# data_planeta_1=pd.read_csv('../result/task1/path_planeta_1.csv',sep=' ')
# data_planeta_2=pd.read_csv('../result/task1/path_planeta_2.csv',sep=' ')
# data_planeta_3=pd.read_csv('../result/task1/path_planeta_3.csv',sep=' ')

# plt.title('Марс')
# plt.plot(data_planeta_0.x-R_start_planeta,data_planeta_0.y,label='1 оборотов')
# plt.plot(data_planeta_1.x-R_start_planeta,data_planeta_1.y,label='100 оборотов')
# plt.plot(data_planeta_2.x-R_start_planeta,data_planeta_2.y,label='500 оборотов')
# plt.plot(data_planeta_3.x-R_start_planeta,data_planeta_3.y,label='1000 оборотов')
# plt.xlim(1.423e11, 1.4234e11)
# plt.xlim(26.99e4,27.01e4)  
# plt.ylim(0.5e11,0.501e11) 
# plt.ylim(0,0.76e1) 
# plt.legend()
# plt.grid()
# plt.show()

R_start_planeta=R_start_planeta+(9.4e6+11.1e3)

data_spytnik_0=pd.read_csv('../result/task1/path_spytnik_0.csv',sep=' ')
data_spytnik_1=pd.read_csv('../result/task1/path_spytnik_1.csv',sep=' ')
data_spytnik_2=pd.read_csv('../result/task1/path_spytnik_2.csv',sep=' ')
data_spytnik_3=pd.read_csv('../result/task1/path_spytnik_3.csv',sep=' ')     

plt.title('Фобос')
plt.plot(data_spytnik_0.x-R_start_planeta,data_spytnik_0.y,label='1 оборотов')
plt.plot(data_spytnik_1.x-R_start_planeta,data_spytnik_1.y,label='100 оборотов')
plt.plot(data_spytnik_2.x-R_start_planeta,data_spytnik_2.y,label='500 оборотов')
plt.plot(data_spytnik_3.x-R_start_planeta,data_spytnik_3.y,label='1000 оборотов')
# plt.xlim(-0.01e11,0.01e11)  
# plt.ylim(0.01e11,0.06e11) 
plt.legend()
plt.grid()
plt.show()