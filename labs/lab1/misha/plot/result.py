import matplotlib.pyplot as plt
import pandas as pd

# Чтение данных
data = pd.read_csv('../out/output.txt', sep=' ')

end = 1660
# plt.plot(data.x1,data.y1)
plt.plot(data.x2[:end],data.y2[:end],label='Марс')
plt.plot(data.x3[:end],data.y3[:end],label='Спутник')
plt.legend()
plt.show()



# for i in range(4):
#     start = i * end
#     end = (i + 1) * end
#     plt.plot(data.x3[start:end], data.y3[start:end], label=f"Марс, участок {i + 1}")
#
# plt.xlabel('x')
# plt.ylabel('y')
# plt.legend()
# plt.grid(True)
# plt.show()