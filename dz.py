import numpy as np
import matplotlib.pyplot as plt

class CubicSpline:

    def __init__(self, x, y):
        indices = sorted(range(len(x)), key=lambda i: x[i])
        x = [x[i] for i in indices]
        y = [y[i] for i in indices]
        self.bounds = x

        self.cnt_spls = len(self.bounds)-1 # = n

        self.pol = [[0] * 4 for i in range(self.cnt_spls)]
        self.fin_pol = [[0] * 4 for i in range(self.cnt_spls)]

        a = []
        a.append(0)
        bt = 2*(x[2]-x[0])
        a.append(-(x[2]-x[1])/bt)
        b = []
        b.append(0)
        b.append(3*((y[2]-y[1])/(x[2]-x[1])-(y[1]-y[0])/(x[1]-x[0]))/bt)
        for i in range(2, self.cnt_spls-1):
            at = x[i]-x[i-1]
            bt = 2*(x[i+1]-x[i-1])
            a.append(-(x[i+1]-x[i])/(at*a[i-1]+bt))
            b.append((3*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1])) - at*b[i-1])/(at*a[i-1]+bt))
        
        self.pol[self.cnt_spls-1][2] = (3*((y[self.cnt_spls]-y[self.cnt_spls-1])/(x[self.cnt_spls]-x[self.cnt_spls-1])-(y[self.cnt_spls-1]-y[self.cnt_spls-2])/(x[self.cnt_spls-1]-x[self.cnt_spls-2])) - (x[self.cnt_spls-1]-x[self.cnt_spls-2])*b[self.cnt_spls-2])/((x[self.cnt_spls-1]-x[self.cnt_spls-2])*a[self.cnt_spls-2]+(2*(x[self.cnt_spls]-x[self.cnt_spls-2])))
        self.pol[self.cnt_spls-1][0] = y[self.cnt_spls-1]
        self.pol[self.cnt_spls-1][3] = -(self.pol[self.cnt_spls-1][2])/(3*(x[self.cnt_spls]-x[self.cnt_spls-1]))
        self.pol[self.cnt_spls-1][1] = (y[self.cnt_spls]-y[self.cnt_spls-1])/(x[self.cnt_spls]-x[self.cnt_spls-1]) - self.pol[self.cnt_spls-1][2]*(x[self.cnt_spls]-x[self.cnt_spls-1]) - self.pol[self.cnt_spls-1][3]*(x[self.cnt_spls]-x[self.cnt_spls-1])*(x[self.cnt_spls]-x[self.cnt_spls-1])

        self.fin_pol[self.cnt_spls-1][0] = self.pol[self.cnt_spls-1][0] - self.pol[self.cnt_spls-1][1]*x[self.cnt_spls-1] + self.pol[self.cnt_spls-1][2]*x[self.cnt_spls-1]*x[self.cnt_spls-1] - self.pol[self.cnt_spls-1][3]*x[self.cnt_spls-1]*x[self.cnt_spls-1]*x[self.cnt_spls-1]
        self.fin_pol[self.cnt_spls-1][1] = self.pol[self.cnt_spls-1][1] - self.pol[self.cnt_spls-1][2]*2*x[self.cnt_spls-1] + self.pol[self.cnt_spls-1][3]*3*x[self.cnt_spls-1]*x[self.cnt_spls-1]
        self.fin_pol[self.cnt_spls-1][2] = self.pol[self.cnt_spls-1][2] - self.pol[self.cnt_spls-1][3]*3*x[self.cnt_spls-1]
        self.fin_pol[self.cnt_spls-1][3] = self.pol[self.cnt_spls-1][3]

        for i in range(self.cnt_spls-2, 0, -1):
            self.pol[i][2] = a[i]*self.pol[i+1][2]+b[i]

        self.pol[0][2] = 0
        
        for i in range(0, self.cnt_spls-1):
            self.pol[i][0] = y[i]
            self.pol[i][1] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (2*self.pol[i][2]+self.pol[i+1][2])*(x[i+1]-x[i])/3
            self.pol[i][3] = (self.pol[i+1][2]-self.pol[i][2])/(3*(x[i+1]-x[i]))

            self.fin_pol[i][0] = self.pol[i][0] - self.pol[i][1]*x[i] + self.pol[i][2]*x[i]*x[i] - self.pol[i][3]*x[i]*x[i]*x[i]
            self.fin_pol[i][1] = self.pol[i][1] - self.pol[i][2]*2*x[i] + self.pol[i][3]*3*x[i]*x[i]
            self.fin_pol[i][2] = self.pol[i][2] - self.pol[i][3]*3*x[i]
            self.fin_pol[i][3] = self.pol[i][3]


    def __call__(self, xa):
        y = []
        for x in xa:
            flg = True
            if x <= self.bounds[0]:
                y.append(self.fin_pol[0][3]*x*x*x + self.fin_pol[0][2]*x*x + self.fin_pol[0][1]*x + self.fin_pol[0][0])
                continue
            for i in range(1, self.cnt_spls+1):
                if x <= self.bounds[i]:
                    y.append(self.fin_pol[i-1][3]*x*x*x + self.fin_pol[i-1][2]*x*x + self.fin_pol[i-1][1]*x + self.fin_pol[i-1][0])
                    flg = False
                    break
            if flg:
                y.append(self.fin_pol[self.cnt_spls-1][3]*x*x*x + self.fin_pol[self.cnt_spls-1][2]*x*x + self.fin_pol[self.cnt_spls-1][1]*x + self.fin_pol[self.cnt_spls-1][0])
        return y


x = [int(x) for x in list(input('Введите значения по x: ')) if x != ' ']
y = [int(x) for x in list(input('Введите значения по y: ')) if x != ' ']

indices = sorted(range(len(x)), key=lambda i: x[i])
x = [x[i] for i in indices]
y = [y[i] for i in indices]

cs = CubicSpline(x, y)

fig, ax = plt.subplots()

xs = np.arange(x[0], x[1], 0.01)
ax.plot(xs, cs(xs), label='y={:.3f}x\u00B3{:+.3f}x\u00B2{:+.3f}x{:+.3f}, x\u2208[{}, {}]'.format(cs.fin_pol[0][3], cs.fin_pol[0][2], cs.fin_pol[0][1], cs.fin_pol[0][0], x[0], x[1]))

for i in range(1, len(x)-1):
    xs = np.arange(x[i], x[i+1], 0.01)
    ax.plot(xs, cs(xs), label='y={:.3f}x\u00B3{:+.3f}x\u00B2{:+.3f}x{:+.3f}, x\u2208({}, {}]'.format(cs.fin_pol[i][3], cs.fin_pol[i][2], cs.fin_pol[i][1], cs.fin_pol[i][0], x[i], x[i+1]))

plt.plot(x, y, 'o', label='Контрольные точки')

plt.grid(True)
plt.legend()
plt.show()
