from vector import Vector
import math

""" 
    this is why i need to plan summers better
"""
class Matrix:
    def __init__(self, *args):
        if len(args) == 0:
            self.values = []
            self.dim = (0, 0)
        else:
            for i in args:
                if type(i) in (type(tuple()), type(list()), type(Vector())):
                    if len(i) != len(args[0]):
                        raise IndexError('row vectors must have same length')
                    continue
                else:
                    raise TypeError('pass 1-D tuple, 1-D list, or Vector as argument')
            self.values = list(list(i) for i in args)
            self.dim = (len(self), len(self[0]))
    def __repr__(self):
        return str(self.values)
    def __iter__(self):
        return self.values.__iter__()
    def __getitem__(self, key):
        return self.values[key]
    def __len__(self):
        return len(self.values)
    def __mul__(self, other):
        if len(self[0]) != len(other):
            return ValueError('incorrect dimensions')
        r, m = [], []
        for i in range(len(self)):
            for j in range(len(other[0])):
                s = 0
                for k in range(len(other)):
                    s += (self[i][k] * other[k][j])
                r.append(s)
            m.append(r)
            r = []
        return Matrix(*m)
    def __rmul__(self, other):
        return self.__mul__(other)
    def __truediv__(self, other):
        raise TypeError('matrix division does not exist')
    def __floordiv__(self, other):
        raise TypeError('matrix division does not exist')
    def transpose(self):
        return Matrix(*map(list, zip(*self)))
    def __add__(self, other):
        if (len(self), len(self[0])) != (len(other), len(other[0])):
            raise IndexError('matricies must have the same dimensions')
        return Matrix(*[[self[i][j] + other[i][j] for j in range(len(self[0]))] for i in range(len(self))])
    def __sub__(self, other):
        if (len(self), len(self[0])) != (len(other), len(other[0])) and type(other) != type(Matrix()):
            raise IndexError('matricies must have the same dimensions')
        return Matrix(*[[self[i][j] - other[i][j] for j in range(len(self[0]))] for i in range(len(self))])      
    @staticmethod
    def __transpose_matrix(m):
        return list(map(list, zip(*m)))
    def inverse(self):
        determinant = self.det()
        if determinant == 0 or self.dim[0] != self.dim[1]:
            raise ZeroDivisionError('this matrix is not invertible')
        if self.dim == (1, 1):
            return 1 / determinant
        n = self.dim[0]
        copy = [i[:] for i in self.values]
        identity = [i[:] for i in [[0] * n] * n]
        for i in range(n):
            identity[i][i] = 1
        indx = list(range(n))
        for pivot in range(n):
            x = 1 / copy[pivot][pivot]
            for j in range(n):
                copy[pivot][j] *= x
                identity[pivot][j] *= x
            for i in indx[0 : pivot] + indx[pivot + 1:]:
                y = copy[i][pivot]
                for j in range(n):
                    copy[i][j] = copy[i][j] - y * copy[pivot][j]
                    identity[i][j] = identity[i][j] - y * identity[pivot][j]
        return Matrix(*identity)
    def det(self):
        if self.dim[0] != self.dim[1]:
            raise IndexError('cannot get determinant of non-square matrix')
        if self.dim == (1, 1):
            return self[0][0]
        try:
            _, l, u = self.lu()
        except:
            return 0
        det1, det2 = 1, 1
        for i in range(self.dim[0]):
            det1 *= l[i][i]
            det2 *= u[i][i]
        return det1 * det2
    def lin_sys(self, *args):
        """ solves **CONSISTENT** linear system with coefficient matrix self and solution set args
            ex:
                a11 * x1 + a12 * x2 + a13 * x3 = b1
                a21 * x1 + a22 * x2 + a23 * x3 = b2
                a31 * x1 + a32 * x2 + a33 * x3 = b3

            Matrix() -> [[a11, a12, a13]
                         [a21, a22, a23]
                         [a31, a32, a33]]

            Vector(args) -> [[b1]
                             [b1]
                             [b3]]

            returns Vector() of solution values for {x1, x2, x3}
        """
        if len(args) == len(self) and type(args[0]) in (type(int()), type(float())):
            return Vector(*[i[0] for i in Matrix.inverse(self) * Matrix(*([i] for i in args))])
        raise ValueError('impossible solution vector')
    def append_col(self, *col):
        #add column vector to the right end of the matrix
        if self.values == [] and all(type(x) in (type(int()), type(float())) for x in col):
            self.values = [[i] for i in col]
            self.dim[1] += 1
            return self
        elif len(col) == self.dim[0] and all(type(x) in (type(int()), type(float())) for x in col):
            self.values = [self.values[i] + [col[i]] for i in range(len(col))]
            self.dim[1] += 1
            return self
        else:
            raise IndexError('column vector length out of bounds')
    def append_row(self, *row):
        #add row vector to the bottom of the matrix
        if self.values == [] and all(type(x) in (type(int()), type(float())) for x in row):
            self.values += [*row]
            self.dim[0] += 1
            return self
        elif len(row) == self.dim[1] and all(type(x) in (type(int()), type(float())) for x in row):
            self.values += [[*row]]
            self.dim[0] += 1
            return self
        else:
            raise IndexError('row vector length out of bounds')
    @staticmethod
    def __pivot(m):
        #helper method for lu()
        n = len(m)
        identity = [[float(i == j) for i in range(n)] for j in range(n)]
        for i in range(n):
            r = max(range(i, n), key = lambda x : abs(m[x][i]))
            if i != r:
                identity[i], identity[r] = identity[r], identity[i]
        return Matrix(*identity)
    def lu(self):
        """ returns identity matrix, L matrix, and U matrix, respectively
            L is lower triangular and U is upper triangular
        """
        if self.dim[0] == self.dim[1]:
            L, U = [[0.0] * self.dim[0] for _ in range(self.dim[0])], [[0.0] * self.dim[0] for _ in range(self.dim[0])]
            P = Matrix.__pivot(self)
            PA = P * self
            for j in range(self.dim[0]):
                L[j][j] = 1.0
                for i in range(j + 1):
                    s = sum(U[k][j] * L[i][k] for k in range(i))
                    U[i][j] = PA[i][j] - s
                for i in range(j, self.dim[0]):
                    s = sum(U[k][j] * L[i][k] for k in range(j))
                    if U[j][j] == 0:
                        raise ZeroDivisionError('LU decomposition does not exist')
                    L[i][j] = (PA[i][j] - s) / U[j][j]
            return Matrix(*P), Matrix(*L), Matrix(*U)
        raise IndexError('square matrix required for LU-decomposition')
    def ref(self):
        try:
            x = self.lu()[2]
            return Matrix(*[[x[i][j] / x[i][i] for j in range(len(x[0]))] for i in range(len(x))])
        except:
            raise Exception('row echelon form does not exist')
    def lin_indep_col(self):
        """ returns true if column vectors of self are linearly independent
            returns false if column vectors of self are linearly dependent
            linear independence means no one column can be expressed as a linear combination of other columns
            column space spans R^m
        """
        if self.dim[0] == self.dim[1]:
            return self.det() != 0
        if self.dim[1] > self.dim[0] or [0] * self.dim[0] in zip(*[i for i in self]):
            return False
        return True
    def lin_indep_row(self):
        """ returns true if row vectors of self are linearly independent
            returns false if row vectors of self are linearly dependent
            linear independence means no one row can be expressed as a linear combination of other columns
            row space spans R^n
        """
        return self.transpose().lin_indep_col()
    def __gt__(self, other):
        if type(other) == type(Matrix()):
            return sum(sum(i for i in j) for j in self) > sum(sum(i for i in j) for j in other)
        raise TypeError('operand not supported')
    def __ge__(self, other):
        if type(other) == type(Matrix()):
            return sum(sum(i for i in j) for j in self) >= sum(sum(i for i in j) for j in other)
        raise TypeError('operand not supported')
    def __le__(self, other):
        if type(other) == type(Matrix()):
            return sum(sum(i for i in j) for j in self) <= sum(sum(i for i in j) for j in other)
        raise TypeError('operand not supported')
    def __lt__(self, other):
        if type(other) == type(Matrix()):
            return sum(sum(i for i in j) for j in self) < sum(sum(i for i in j) for j in other)
        raise TypeError('operand not supported')
    def __eq__(self, other):
        if type(other) == type(Matrix()) and self.dim == other.dim:
            return all(self[i] == other[i] for i in range(self.dim[0]))
        return False
    @staticmethod
    def __norm_squared(v):
        #helper method for gs_process()
        return sum([i ** 2 for i in v])
    @staticmethod
    def __subtract(v1, v2):
        #helper method for gs_process()
        v1 = [i - j for i, j in zip(v1, v2)]
        return v1
    @staticmethod
    def __scalar_multiplication(v1, scalar):
        #helper method for gs_process()
        return [i * scalar for i in v1]
    @staticmethod
    def __inner_product(v1, v2):
        #helper method for gs_process()
        return sum([i * j for i, j in zip(v1, v2)])
    def gs_process(self):
        """ returns new orthoganal and orthonormal basis, respectively
            frequently used in QR decomposition
        """
        a = self.transpose().values
        orthoganal, orthonormal = [a[0]], []
        for i in range(1, len(a)):
            x = a[i]
            for j in range(i):
                multiple = Matrix.__inner_product(a[i], orthoganal[i - j - 1]) / Matrix.__norm_squared(orthoganal[i - j - 1])
                x = Matrix.__subtract(x, Matrix.__scalar_multiplication(orthoganal[i - j - 1], multiple))
            orthoganal.append(x)
        orthonormal = [[j / math.sqrt(Matrix.__norm_squared(i)) for j in i] for i in orthoganal]
        return Matrix(*orthoganal), Matrix(*orthonormal)
    def qr(self):
        """ the first return matrix Q has column vectors orthonormal basis from Gram-Schmidt process
            the second return matrix R is an upper triangular matrix such that Q * R is self
        """
        a = self.transpose().values
        q = list(self.gs_process()[1])
        r = []
        for i in range(len(q)):
            x = [0] * i
            for j in range(i, len(q)):
                x.append(Matrix.__inner_product(q[i], a[j]))
            r.append(x)
        return Matrix(*q).transpose(), Matrix(*r)
    def ls_regression(self, degree, *response_args):
        """ creates n x m matrix of polynomial explanatory variables (x-coordinates)
            from matrix of x-coordinates, of given degree
            ex: [
                [1, x1, x1^2, ... x1^m],
                [1, x2, x2^2, ... x2^m],
                [.  ..  ....  ... ....],
                [1, xn, xn^2, ... xn^m]
            ]
            and Matrix() of response variables
            ex: [
                [y1],
                [y2],
                ....,
                [yn]
            ]
            returns Vector() of cooefficients for best fit polynomial of degree m
        """
        if len(response_args) != self.dim[0] or not all(type(i) in (type(int()), type(float())) for i in response_args):
            raise ValueError('invalid argument length or type')
        if self.dim[1] != 1:
            raise Exception('this method does not apply to this object')
        inpt = Matrix(*[[i[0] ** j for j in range(degree + 1)] for i in self])
        r = Matrix(*[[i] for i in response_args])
        return Vector(*[i for j in inpt.qr()[1].inverse() * inpt.qr()[0].transpose() * r for i in j])
    def trace(self):
        return sum(self[i][i] for i in range(self.dim[0]))
    def eigenvalues(self):
        if self.dim != (2, 2):
            raise NotImplementedError('only available for 2x2 matrices for now')
        x = self.trace()
        y = self[0][0] * self[1][1] - self[0][1] * self[1][0]
        return [0.5 * (x - math.sqrt(x ** 2 - 4 * y)), 0.5 * (x + math.sqrt(x ** 2 - 4 * y))]
    def eigenspace_matrices(self):
        x = self.eigenvalues()
        return [Matrix([i, 0], [0, i]) - self for i in x]
    def eigenspace_basis(self, normalize = False):
        """ returns list of Vector()s where the vectors span the eigenspace of the matrix
            does not work on anything other than 2x2 matricies
        """
        if self.dim == (2, 2):
            x = self.eigenspace_matrices()
            rel = [i[0] for i in x]
            output = []
            for i in rel:
                out = [None] * len(i)
                y = [abs(j) for j in i]
                gr = i[y.index(max(y))]
                out[y.index(max(y))] = 1.0
                for j in out:
                    if j == None:
                        try:
                            out[out.index(j)] = -1 * gr / i[out.index(j)]
                        except:
                            return [Vector(0, 0), Vector(0, 0)]
                output.append(Vector(*out))
            if normalize:
                return [i.normalize() for i  in output]
            return output
        raise NotImplementedError('haven\'t bothered to do more than 2x2 yet')
    def orth_diagonalize(self):
        P_T = Matrix(*self.eigenspace_basis(normalize=True))
        P = P_T.transpose()
        return P_T * self * P, P_T, P
