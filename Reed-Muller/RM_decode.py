import numpy as np
import itertools
from RM_recode import RM


class RM_low_radix:
    def __init__(self, R, M, C):
        self.R, self.M, self.C = R, M, C

    def start(self):
        return self.Decode()

    def Decode(self):
        w = 2*self.C - np.ones(len(self.C))
        for i in range(self.M):
            H = self.Kro_p((i+1), self.M)
            w = np.dot(w, H)

        address, num = -1, abs(w[0])
        for i in range((2 ** self.M)):
            if abs(w[i]) > num:
                address, num = i, abs(w[i])
        U = self.R_message(address, w[address])
        return U

    def Kro_p(self, j, m):
        H = [[1, 1],
             [1, -1]]
        H = np.array(H)
        I1 = np.eye((2 ** (m - j)))
        I2 = np.eye((2 ** (j - 1)))
        R = np.kron(np.kron(I1, H), I2)
        x, y = R.shape
        for i in range(x):
            for j in range(y):
                if R[i, j] == -0:
                    R[i, j] = 0
        return R

    def Binary_trans(self, num):
        r = []
        while (not (num == 0)):
            ret = num % 2
            num = num // 2
            r = [ret] + r
        return r

    def G_k(self):
        k = 0
        for i in range(self.R + 1):
            k = k + np.math.factorial(self.M)/(np.math.factorial(self.M - i)*np.math.factorial(i))
        return int(k)

    def R_message(self, address, w):
        k = self.G_k()
        U = np.zeros(k)
        if w <= 0:
            U[0] = 0
        else:
            U[0] = 1
        B = self.Binary_trans(address)
        l = len(B)
        for i in range(l - 1, -1, -1):
            U[l - i] = B[i]
        return U


class Matrix:
    def __init__(self, r, m):
        self.r, self.m = r, m

    def start(self):
        self.G(1, self.m)
        self.G_matrix()
        return self.RM

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


class RM_high_radix:
    def __init__(self, R, M, C):
        self.R, self.M, self.C = R, M, C

    def start(self):
        U = self.Decode()
        return U

    def Decode(self):
        g = Matrix(self.R, self.M)
        self.G = g.start()
        self.G_sup = self.G[:(self.M + 1), :]
        self.G_sup = np.ones(self.G_sup.shape) - self.G_sup
        self.G_borrow = self.G.copy()
        self.X = []
        U_mes = []
        for i in range(1, self.M + 1):
            self.X += [i]
        for i in range(self.R, 0, -1):
            c = list(itertools.combinations(self.X, i))
            l1, l = len(c), len(c[0])
            cons_range = self.M - i
            message = []
            for j1 in range(l1):
                X_0 = self.G_rest_sup(l, c, j1)
                if len(X_0) == 2:
                    c_data = X_0
                else:
                    c_data = self.Delete_same(X_0, cons_range)
                if len(c_data) == 2:
                    symbol = self.R_message_one(c_data[0])
                    message = message + [symbol]
                else:
                    symbol = self.R_mess_more(c_data)         # processing multiple element
                    message = message + [symbol]
            U_mes = message + U_mes
            self.Renew_C(message)                             # get the newer received symbols
        U_mes = np.array(([self.First_symbol()] + U_mes))
        return U_mes

    def First_symbol(self):
        d = np.count_nonzero(self.C)
        if d > (len(self.C)/2):
            return 1
        else:
            return 0

    def Delete_same(self, X_0, cons_range):
        c_data = list(itertools.combinations(X_0, cons_range))
        l_x = 0
        while l_x < len(c_data):
            Tab = 0
            for l_y in range(len(c_data[0])):
                if (-c_data[l_x][l_y]) in c_data[l_x]:
                    Tab = 1
                    break
            if Tab == 1:
                c_data.remove(c_data[l_x])
            else:
                l_x += 1
        return c_data

    def G_rest_sup(self, l, c, j1):
        X_0 = self.X.copy()
        for j2 in range(l):
            X_0.remove(c[j1][j2])
        for j3 in range(len(X_0)):
            X_0 = X_0 + [-X_0[j3]]
        return X_0

    def Renew_C(self, message):
        x_G_borrow, y_G_borrow = self.G_borrow.shape
        l_message = len(message)
        G_0 = self.G_borrow[x_G_borrow - l_message:, :]
        s = np.dot(np.array(message), G_0) % 2
        self.C = (self.C + s) % 2
        self.G_borrow = self.G_borrow[:x_G_borrow - l_message, :]

    def R_message_one(self, a):
        one_num = zero_num = 0
        m1 = (self.C @ self.G[a]) % 2
        m2 = (self.C @ self.G_sup[a]) % 2
        if m1 == 0:
            zero_num += 1
        else:
            one_num += 1
        if m2 == 0:
            zero_num += 1
        else:
            one_num += 1
        if zero_num >= one_num:
            symbol = 0
        else:
            symbol = 1
        return symbol

    def R_mess_more(self, c_data):
        one_num = zero_num = 0
        x_c_data = len(c_data)
        y_c_data = len(c_data[0])
        for i_c_data in range(x_c_data):
            if c_data[i_c_data][0] < 0:
                c_temporary = self.G_sup[-(c_data[i_c_data][0])]
            else:
                c_temporary = self.G[c_data[i_c_data][0]]
            for j_c_data in range(y_c_data):
                if c_data[i_c_data][j_c_data] < 0:
                    c_temporary = c_temporary * self.G_sup[-(c_data[i_c_data][j_c_data])]
                else:
                    c_temporary = c_temporary * self.G[c_data[i_c_data][j_c_data]]
            m1 = (self.C @ c_temporary) % 2
            if m1 == 0:
                zero_num += 1
            else:
                one_num += 1
        if zero_num >= one_num:
            symbol = 0
        else:
            symbol = 1
        return symbol
