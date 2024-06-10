import func2 as f
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    points_count = 20
    border = 50
    lines_count = 100
    line_origins_real = np.linspace(-border, border, lines_count)
    line_origins_imaginary = np.linspace(-border * 1j, border * 1j, lines_count)
    lines_real = []
    lines_imaginary = []
    for origin_real in line_origins_real:
        for point in np.linspace(origin_real - 1j * border, origin_real + 1j * border, points_count):
            lines_real.append(point)
    for origin_imaginary in line_origins_imaginary:
        for point in np.linspace(origin_imaginary - border, origin_imaginary + border, points_count):
            lines_imaginary.append(point)
    lines_real = np.array(lines_real)
    lines_imaginary = np.array(lines_imaginary)

    r = f.dlo([1, 2, 3], [-3, -2, -1])
    dlo_image_real = f.apply_dlo(r, lines_real)
    dlo_image_imaginary = f.apply_dlo(r, lines_imaginary)

    plt.scatter([np.real(dlo_image_real), np.real(dlo_image_imaginary)],
                [np.imag(dlo_image_real), np.imag(dlo_image_imaginary)], alpha=0.3)
    plt.show()