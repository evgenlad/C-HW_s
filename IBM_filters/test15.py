import numpy as np
import func2 as f


p = f.OM * 1j
t = f.TAU
c = f.C
n = f.N
m = f.M
p[1, 1] *= 2

[e1, e2, e3, e4, e5, e6] = f.special_chars()
[u1, u2, u3, u4, u5, u6] = f.special_points(p)
w = [w1, w2, w3, w4, w5, w6] = f.my_polig_vertices(n, t, m)

#print(f.stand_polyg(-1 + n * 1j * t, n, t))
for i in range(1, 7):
    w[i - 1] = f.stand_polyg_new(w[i - 1]) / c
    print(f"u{i}[0] = {f.special_points(p)[i - 1][0]},\t vertex{i} = {w[i - 1]}")
