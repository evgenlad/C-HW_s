import time
import func2 as f


# начальное время
t = f.TAU
n = f.N
p = 1j * f.OM
c = f.C
point = 1 - n * (t / 2) * 1j + 0.5
start_time = time.time()
#r = f.x_polygon_to_hplane(point, p, c)
# код, время выполнения которого нужно измерить\
r = f.x_polygon_to_hplane(point, p, c)
# конечное время
end_time = time.time()

# разница между конечным и начальным временем
elapsed_time = end_time - start_time
print('Elapsed time: ', elapsed_time)
