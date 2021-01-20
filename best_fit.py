from vector import Vector
from matrix import Matrix
import matplotlib.pyplot as plt
import numpy as np

def linreg(points, degree):
    def prep(points):
        return Matrix(*[[i[0]] for i in points]), [i[1] for i in points]
    x_coord, y_coord = prep(points)
    return Vector(*[round(i, 10) + 0 for i in (x_coord.ls_regression(degree, *y_coord))])

def main(graph=False):
    x = [
        (0, -1),
        (1, 1),
        (2, 3),
        (3, 60),
    ]
    coefficients = linreg(x, 2)
    points = np.linspace(min([i[0] for i in x]), max(i[0] for i in x), 100)
    print(coefficients)
    if graph:
        fig, ax = plt.subplots()
        ax.plot(points, sum([coefficients[i] * (points ** i) for i in range(len(coefficients))]), label='best fit polynomial')
        ax.legend()
        for i in x:
            ax.plot(i[0], i[1], 'bo')
        plt.show()

if __name__ == '__main__':
    main(True)
