import math

class Vector:
    def __init__(self, *args):
        if len(args) == 0:
            self.values = (0,0)
        else:
            self.values = args
    def __iter__(self):
        return self.values.__iter__()
    def __repr__(self):
        return str(self.values)
    def __add__(self, other):
        if len(self) != len(other):
            raise IndexError('Vectors must have the same length')
        return Vector(*tuple(i + j for i, j in zip(self, other)))
    def __mul__(self, other):
        if type(other) == type(self):
            return self.inner(other)
        elif type(other) == type(1) or type(other) == type(1.0):
            return Vector(*tuple(i * other for i in self))
    def __rmul__(self, other):
        return self.__mul__(other)
    def __sub__(self, other):
        if len(self) != len(other):
            raise IndexError('Vectors must have the same length')
        return Vector(*tuple(i - j for i, j in zip(self, other)))
    def __truediv__(self, other):
        raise TypeError('Vector division does not exist')
    def __floordiv__(self, other):
        raise TypeError('Vector division does not exist')
    def norm(self):
        return math.sqrt(sum(i ** 2 for i in self))
    def normalize(self):
        x = self.norm()
        return Vector(*tuple(i / x for i in self))
    def inner(self, other):
        if len(self) != len(other):
            raise IndexError('Vectors must have the same length')
        return sum(i * j for i, j in zip(self, other))
    def matrix_multiply(self, matrix):
        """ takes row of matrix and multiplies by vector
            example:
            [1, 2, 3] * [[1, 2, 3], [4, 5, 6], [7, 8, 9]] = 
            [[1, 2, 3] * [1, 2, 3], [1, 2, 3] * [4, 5, 6], [1, 2, 3] * [7, 8, 9]] =
            [14, 32, 50]
        """
        if not all(len(i) == len(self) for i in matrix):
            raise ValueError('dimensions do not match')
        return Vector(*tuple(Vector(*i) * self for i in matrix))
    def __len__(self):
        return len(self.values)
    def __getitem__(self, key):
        return self.values[key]
    def rotate_2d_rev(self, theta = 0):
        #counterclockwise rotation about origin
        if type(theta) in (type(int()), type(float())):
            if len(self.values) == 2:
                return self.matrix_multiply([[math.cos(theta), math.sin(theta)], [-1 * math.sin(theta), math.cos(theta)]])
            raise TypeError('use a valid two dimensional vector')
        raise TypeError('theta must be an integer or float in radians')
    def rotate_2d(self, theta = 0):
        #clockwise rotation about origin
        return self.rotate_2d_rev(-1 * theta)
    def orth_proj(self, axis_index):
        if axis_index not in range(len(self)):
            raise IndexError('specified projection axis out of bounds')
        return Vector(*[self[i] for i in range(axis_index)] + [0] + [self[i] for i in range(axis_index + 1, len(self))])