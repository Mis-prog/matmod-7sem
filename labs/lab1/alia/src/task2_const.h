// Физические константы
struct Constants {
    static constexpr double G = 6.67e-11;    // гравитационная постоянная
    static constexpr double M1 = 2.0e30;     // масса звезды (кг)
    static constexpr double M2 = 6.0e24;     // масса планеты (кг)
    static constexpr double M3 = 7.3e22;     // масса астероида (кг)
    static constexpr double R1 = 696340e3;   // радиус звезды (м)
    static constexpr double R2 = 6378e3;     // радиус планеты (м)
    static constexpr double R3 = 1737e3;     // радиус астероида (м)
    static constexpr double R12 = 150e9;     // начальное расстояние звезда-планета (м)
    static constexpr double R23 = 384e6;     // начальное расстояние планета-астероид (м)
    static constexpr double U2 = 30e3;       // начальная скорость планеты (м/с)
    static constexpr double U3 = 1e3;     // начальная скорость астероида (м/с)

    static constexpr double T = 2400.0; // время работы двигателя (с)
    static constexpr double H = 300e3; // высота орбиты (м)
    static constexpr double M0 = 120.0; // масса полезной назрузки (кг)
    static constexpr double U = 3060.0; // скорость истечения (м/c)
    static constexpr double koef = 0.025;
};