import func2 as f
import numpy as np
import matplotlib.pyplot as plt


t = f.TAU
c = f.C
n = f.N
m = f.M
p = f.OM * 1j
p[1, 1] *= 2

w_arr = np.linspace(-1 + m * 1j * t, 1 + m * 1j * t, 10)
for w in w_arr:
    z1 = f.take_limit_f(lambda x: f.x_polygon_to_hplane(f.stand_polyg(x, n, t), p, c),
                 w, t * 1j)
    z2 = f.take_limit_f(lambda x: f.x_polygon_to_hplane(f.stand_polyg(x, n, t), p, c),
                        w, -t * 1j)
    print(z1, z2)
