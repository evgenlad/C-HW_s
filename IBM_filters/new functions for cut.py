# что хочу здесь? функцию, которая даёт ответ, есть ли ветвление в точке
# то есть, должа браться эта точка и до нее реккурентно строятся две ветви
# сначала сделать заново с массивами, посмотреть что получается
# построение функции:
# нахождение h1 смещения по u_1, с помощью него и фукнции stay_on_surface_one_step
# сдвигаться от точек u4, u5 к данной
# если в процессе сдвига u_2 станут близкими, то ветвление кончилось.
import numpy as np


def branching_function(w, p, c, step_num):
    branching_check = True

    [e1, e2, e3, e4, e5, e6] = special_chars()
    [u01, u02, u03, u04, u05, u06] = special_points(p)

    point_u1 = w / c
    h1 = (point_u1 - u04[0]) / step_num     # шаг по первой координате на якобиане

    def theta_u4(u, p):
        return theta35(u04 + u, p)

    def theta_u4_der1(u, p):
        return theta35_der1(u04 + u, p)

    def theta_u4_der2(u, p):
        return theta35_der2(u04 + u, p)

    def theta_u5(u, p):
        return theta35(u05 + u, p)

    def theta_u5_der1(u, p):
        return theta35_der1(u05 + u, p)

    def theta_u5_der2(u, p):
        return theta35_der2(u05 + u, p)

    f0_4 = lambda x: theta_u4(x, p)
    f_der1_4 = lambda x: theta_u4_der1(x, p)
    f_der2_4 = lambda x: theta_u4_der2(x, p)

    f0_5 = lambda x: theta_u5(x, p)
    f_der1_5 = lambda x: theta_u5_der1(x, p)
    f_der2_5 = lambda x: theta_u5_der2(x, p)

    u_curr_4 = np.array([0, 0], complex)
    u_curr_5 = np.array([0, 0], complex)
    for i in range(step_num):
        # u_curr_4(5) - текущая точка, которая в цикле будет использоваться как стартовая
        # u_new_4(5) - новая точка, которая будет получена с помощью метода
        u_new_4 = stay_on_surface_one_step(u_curr_4, h1, f0_4, f_der1_4, f_der2_4)[1]
        u_new_5 = stay_on_surface_one_step(u_curr_5, h1, f0_5, f_der1_5, f_der2_5)[1]

        # проверим, насколько разветвление прекращено
        if np.abs(u_new_4[1] - u_new_5[1]) < 1e-10:
            # ветви сошлись!
            branching_check = False
