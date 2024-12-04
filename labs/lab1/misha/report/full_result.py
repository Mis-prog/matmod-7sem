import matplotlib.pyplot as plt
import pandas as pd


def full_drawe():
    data = pd.read_csv('../res_task1/path_full.csv', sep=' ')
    plt.plot(data.x2, data.y2, label='Марс')
    # plt.plot(data.x3, data.y3, ls='--', label='Фобос')
    plt.legend()
    plt.grid(True)
    plt.show()


def part_spytnik_drawe():
    data1 = pd.read_csv('path_spytnik_0.csv', sep=' ')
    data2 = pd.read_csv('path_spytnik_1.csv', sep=' ')
    data3 = pd.read_csv('path_spytnik_2.csv', sep=' ')
    data4 = pd.read_csv('path_spytnik_3.csv', sep=' ')
    plt.plot(data1.x, data1.y, ls='--', label='1')
    plt.plot(data2.x, data2.y, ls='--', label='100')
    plt.plot(data3.x, data3.y, ls='--', label='500')
    plt.plot(data4.x, data4.y, ls='--', label='1000')
    plt.legend()
    plt.grid(True)
    plt.show()

def part_planeta_drawe():
    data1 = pd.read_csv('path_planeta_0.csv', sep=' ')
    data2 = pd.read_csv('path_planeta_1.csv', sep=' ')
    data3 = pd.read_csv('path_planeta_2.csv', sep=' ')
    data4 = pd.read_csv('path_planeta_3.csv', sep=' ')
    plt.plot(data1.x, data1.y, ls='--', label='1')
    plt.plot(data2.x, data2.y, ls='--', label='100')
    plt.plot(data3.x, data3.y, ls='--', label='500')
    plt.plot(data4.x, data4.y, ls='--', label='1000')
    plt.legend()
    plt.grid(True)
    plt.show()

def part_drawe():
    data1 = pd.read_csv('path_spytnik_0.csv', sep=' ')
    data2 = pd.read_csv('path_spytnik_1.csv', sep=' ')
    data3 = pd.read_csv('path_spytnik_2.csv', sep=' ')
    data4 = pd.read_csv('path_spytnik_3.csv', sep=' ')
    plt.plot(data1.x, data1.y, ls='--', label='1')
    plt.plot(data2.x, data2.y, ls='--', label='100')
    plt.plot(data3.x, data3.y, ls='--', label='500')
    plt.plot(data4.x, data4.y, ls='--', label='1000')
    data5 = pd.read_csv('path_planeta_0.csv', sep=' ')
    data6 = pd.read_csv('path_planeta_1.csv', sep=' ')
    data7 = pd.read_csv('path_planeta_2.csv', sep=' ')
    data8 = pd.read_csv('path_planeta_3.csv', sep=' ')
    plt.plot(data5.x, data5.y, ls='--', label='1_pl')
    plt.plot(data6.x, data6.y, ls='--', label='100_pl')
    plt.plot(data7.x, data7.y, ls='--', label='500_pl')
    plt.plot(data8.x, data8.y, ls='--', label='1000_pl')
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
    data1 = pd.read_csv('path_spytnik_0.csv', sep=' ')
    data2 = pd.read_csv('path_spytnik_1.csv', sep=' ')
    data3 = pd.read_csv('path_spytnik_2.csv', sep=' ')
    data4 = pd.read_csv('path_planeta_3.csv', sep=' ')

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

def thousandth_branch_drawe():
    spytnik = pd.read_csv('path_spytnik_3.csv', sep=' ')
    planeta = pd.read_csv('path_planeta_3.csv',sep=' ')
    plt.plot(spytnik.x, spytnik.y, label='Фобос')
    plt.plot(planeta.x, planeta.y, label='Марс')
    plt.legend()
    plt.show()

full_drawe()