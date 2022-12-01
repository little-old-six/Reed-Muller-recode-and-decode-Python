import numpy as np
import itertools
import random


class RM:
    def __init__(self, r, m, n):
        self.r, self.m, self.n = r, m, n

    def start(self):
        return self.Recode()

    def G_k(self):
        k = 0
        for i in range(self.r + 1):
            k = k + np.math.factorial(self.m) / (np.math.factorial(self.m - i) * np.math.factorial(i))
        return k

    def Recode(self):
        error = 2**(self.m-self.r-1)-1
        print('The RM(%d,%d) can correct max %d errors' % (self.r, self.m, error))
        if self.n > error:
            print('Warning! Errors exceed the error correction capability of the algorithm!')
        k = int(self.G_k())
        U = np.random.randint(0, 2, k)
        print('Origional message symbol：', U)
        self.G(1, self.m)
        if self.r > 1:
            self.G_matrix()
        M = np.dot(U, self.RM) % 2
        r = random.sample(range(0, (2**self.m)), self.n)
        print('The number of noise is :', self.n, '   The location of noise:', r)
        for i in (r):
            M[i] = 1 - M[i]
        print('Received codeword symbol:', M)
        return U, M

    def G(self, r, m):
        if (r == 0) and (m != 0):
            self.RM = np.ones((1, (2 ** m)))
            return self.RM
        elif r == m:
            a = self.G(r - 1, m)
            l = np.size(a, 1)
            b = np.concatenate((np.zeros((1, (l - 1))), [[1]]), axis=1)
            self.RM = np.concatenate((a, b), axis=0)
            return self.RM

        elif r != m:
            a = self.G(r, (m - 1))
            b = self.G((r - 1), (m - 1))
            x, y = b.shape
            c = np.zeros((x, y))
            A = np.concatenate((a, a), axis=1)
            B = np.concatenate((c, b), axis=1)
            self.RM = np.concatenate((A, B), axis=0)
            return self.RM

    def G_matrix(self):
        H = np.array(self.RM)
        H_0 = H.copy()
        H = np.delete(H, 0, axis=0)
        H = H.tolist()
        for i in range(2, self.r + 1):
            c = list(itertools.combinations(H, i))
            l1 = len(c)
            l = len(c[0])
            for j1 in range(l1):
                x_lab = c[j1][0]
                for j2 in range(l):
                    x_lab = np.multiply(x_lab, c[j1][j2])
                H_0 = np.r_[H_0, [x_lab]]
        self.RM = H_0