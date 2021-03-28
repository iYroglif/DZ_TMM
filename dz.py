import numpy as np
#from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

class CubicSpline:

    def __init__(self, x, y):
        self.bounds = x

        cnt_spls = len(self.bounds)-1

        self.pol = [[0] * 4 for i in range(cnt_spls)]

        a = [0]
        bt = 2*(x[2]-x[0])
        a.append(-(x[2]-x[1])/bt)
        b = [0]
        b.append(3*((y[2]-y[1])/(x[2]-x[1])-(y[1]-y[0])/(x[1]-x[0]))/bt)
        for i in range(2, cnt_spls-1):
            at = x[i]-x[i-1]
            bt = 2*(x[i+1]-x[i-1])
            a.append(-(x[i+1]-x[i])/(at*a[i-1]+bt))
            b.append((3*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1])) - at*b[i-1])/(at*a[i-1]+bt))
        
        self.pol[cnt_spls-1][2] = (3*((y[cnt_spls]-y[cnt_spls-1])/(x[cnt_spls]-x[cnt_spls-1])-(y[cnt_spls-1]-y[cnt_spls-2])/(x[cnt_spls-1]-x[cnt_spls-2])) - (x[cnt_spls-1]-x[cnt_spls-2])*b[cnt_spls-1])/((x[cnt_spls-1]-x[cnt_spls-2])*a[cnt_spls-1]+bt)
        self.pol[cnt_spls-1][0] = y[cnt_spls-1]
        self.pol[cnt_spls-1][1] = (y[cnt_spls]-y[cnt_spls-1])/(x[cnt_spls]-x[cnt_spls-1]) - 2*self.pol[cnt_spls-1][2]*(x[cnt_spls]-x[cnt_spls-1])/3
        self.pol[cnt_spls-1][3] = -self.pol[cnt_spls-1][2]/(3*(x[cnt_spls]-x[cnt_spls-1]))

        for i in range(cnt_spls-2, -1, -1):
            self.pol[i][2] = a[i]*self.pol[i+1][2]+b[i]
            self.pol[i][0] = y[i]
            self.pol[i][1] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (2*self.pol[i][2]+self.pol[i+1][2])*(x[i+1]-x[i])/3
            self.pol[i][3] = (self.pol[i+1][2]-self.pol[i][2])/(3*(x[i+1]-x[i]))


x = [int(x) for x in list(input('Введите значения по x: ')) if x != ' ']
y = [int(x) for x in list(input('Введите значения по y: ')) if x != ' ']

cs = CubicSpline(x, y)
print(cs.pol)

#cs = CubicSpline(x, y)
#xs = np.arange(-5, 5, 0.1)

#fig, ax = plt.subplots(figsize=(10,10)) 
#ax.plot(xs, cs(xs))
#plt.show()


