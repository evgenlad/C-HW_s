import numpy as np
import func2 as f


# здесь буду искать место раздвоения на отрезке M*TAU
t = f.TAU
c = f.C
n = f.N
m = f.M
p = 1j * f.OM
p[1, 1] *= 1.5

points_count = 5
interval = np.linspace(m * t * 1j, 1 + m * t * 1j, points_count)
diffs = np.zeros(points_count, complex)
for i in range(points_count):
    print('|', end='')
    diffs[i] = (f.take_limit_f(lambda x: f.whole_pol_to_plane(x, n, t, p, c), interval[i], 1j * t) -
                f.take_limit_f(lambda x: f.whole_pol_to_plane(x, n, t, p, c), interval[i], -1j * t))
print('\n', np.max(np.abs(diffs)))
