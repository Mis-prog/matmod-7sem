import matplotlib.pyplot as plt
import pandas as pd


def full_drawe():
    data = pd.read_csv('output.csv', sep=' ')
    plt.plot(data.x2, data.y2, label='Марс')
    # plt.plot(data.x3, data.y3, ls='--', label='Фобос')
    plt.legend()
    plt.grid(True)
    plt.show()


def part_drawe():
    data1 = pd.read_csv('out_0.csv', sep=' ')
    data2 = pd.read_csv('out_1.csv', sep=' ')
    data3 = pd.read_csv('out_2.csv', sep=' ')
    data4 = pd.read_csv('out_3.csv', sep=' ')
    plt.plot(data1.x, data1.y, label='1')
    plt.plot(data2.x, data2.y, label='100')
    plt.plot(data3.x, data3.y, label='500')
    plt.plot(data4.x, data4.y, label='1000')
    plt.legend()
    plt.grid(True)
    plt.show()


def part_drawe_detail():
    data1 = pd.read_csv('out_0.csv', sep=' ')
    data2 = pd.read_csv('out_1.csv', sep=' ')
    data3 = pd.read_csv('out_2.csv', sep=' ')
    data4 = pd.read_csv('out_3.csv', sep=' ')
    plt.plot(data1.x, data1.y, label='1')
    plt.plot(data2.x, data2.y, label='100')
    plt.plot(data3.x, data3.y, label='500')
    plt.plot(data4.x, data4.y, label='1000')
    plt.legend()
    plt.grid(True)
    plt.show()


def part_drawe_bias():
    data1 = pd.read_csv('out_0.csv', sep=' ')
    data2 = pd.read_csv('out_1.csv', sep=' ')
    data3 = pd.read_csv('out_2.csv', sep=' ')
    data4 = pd.read_csv('out_3.csv', sep=' ')

    # Находим минимальное значение x среди всех наборов данных
    min_x = min(data1.x.min(), data2.x.min(), data3.x.min(), data4.x.min())

    # Смещаем значения x, вычитая минимальное значение
    plt.plot(data1.x - min_x, data1.y, label='1')
    plt.plot(data2.x - min_x, data2.y, label='100')
    plt.plot(data3.x - min_x, data3.y, label='500')
    plt.plot(data4.x - min_x, data4.y, label='1000')

    plt.legend()
    plt.grid(True)
    plt.show()


part_drawe_bias()
