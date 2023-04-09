import numpy as np
import decimal
import matplotlib.pyplot as plt

#диффур
def func_l4v3(x,y,z):
    return (2*y+4*(x**2)*np.exp(x**2))

#аналитическое решение
def solution(x):
    return (np.exp(x**2)+np.exp(x*np.sqrt(2))+np.exp(-1*x*np.sqrt(2)))


def Eiler(f, y0, z0, a, b, h):
    y = [y0]
    z = [z0]
    x = []
    i = 0
    for xi in np.arange(a,b+h,h):
        x.append(xi)
        z.append(z[i]+h*f(x[i],y[i],z[i]))
        y.append(y[i]+h*z[i])
        i += 1
    z.pop()
    y.pop()
    return x,y,z


def runge_kutt(f, y0, z0, a, b, h):
    y = [y0]
    z = [z0]
    x = []
    i = 0
    for xi in np.arange(a,b+h,h):
        k = np.zeros(4)
        l = np.zeros(4)
        x.append(xi)
        for j in range(4):
            if j == 0:
                    l[j] = h * f(x[i], y[i], z[i])
                    k[j] = h * z[i]
            elif j == 3:
                    l[j] = h * f(x[i] + h, y[i] + k[j-1], z[i] + l[j-1])
                    k[j] = h * (z[i] + l[j-1])
            else:
                    l[j] = h * f(x[i] + 0.5*h, y[i] + 0.5 * k[j-1], z[i] + 0.5*l[j-1])
                    k[j] = h * (z[i] + 0.5*l[j-1])
        z.append(z[i] + (l[0] + 2 * (l[1] + l[2]) + l[3]) / 6)
        y.append(y[i] + (k[0] + 2 * (k[1] + k[2]) + k[3]) / 6)
        i+=1
    y.pop()
    z.pop()
    return x,y,z


def adams(f, y0, z0, a, b, h):
    x, y, z = runge_kutt(f, y0, z0, a, b, h)
    y = y[:4]
    z = z[:4]
    for i in range(3, len(x) - 1):
        f1 = f(x[i], y[i], z[i])
        f2 = f(x[i-1], y[i-1], z[i-1])
        f3 = f(x[i-2], y[i-2], z[i-2])
        f4 = f(x[i-3], y[i-3], z[i-3])
        z.append(z[i] + h / 24 * (55*f1 - 59*f2 + 37*f3 - 9*f4))
        y.append(y[i] + h / 24 * (55*z[i] - 59*z[i-1] + 37*z[i-2] - 9*z[i-3]))
    return x, y, z


def get_error(x,y,z,method, sol):
    print(method)
    a = 2
    for i in range(len(x)):
        print("x = {0:.15f};".format(x[i],), "y = {0:.15f},".format(y[i]), "y_ист = {0:.15f}".format(sol(x[i])), "eps: {0:.15f}".format(np.abs(y[i]-sol(x[i]))))


def pres(x,sol):
    pres = []
    for i in range(len(x)):
        pres.append(sol(x[i]))
    return pres


if __name__ == "__main__":
    
    y0 = 3
    z0 = 0
    a = 0
    b = 1
    h = 0.1
    x, y, z = Eiler(func_l4v3, y0, z0, a, b, h)
    get_error(x,y,z,"Метод Эйлера",solution)
    y_pres = pres(x,solution)
    fig, ax = plt.subplots()
    ax.grid()
    ax.plot(x,y_pres, color="red", label="Точное решение")
    ax.plot(x,y, label="Метод Эйлера", color="blue")
    ax.legend(loc='best')
    plt.show()
    x, y, z = runge_kutt(func_l4v3, y0, z0, a, b, h)
    get_error(x,y,z,"Метод Рунге-Кутты",solution)
    fig, ax = plt.subplots()
    ax.grid()
    ax.plot(x,y_pres, label="Точное решение", color="red")
    ax.plot(x,y, label="Метод Рунге-Кутты", color="blue")
    ax.legend(loc='best')
    plt.show()
    x, y, z = adams(func_l4v3, y0, z0, a, b, h)
    get_error(x,y,z,"Метод Адамса",solution)
    fig, ax = plt.subplots()
    ax.grid()
    ax.plot(x,y_pres, label="Точное решение", color="red")
    ax.plot(x,y, label="Метод Адамса", color="blue")
    ax.legend(loc='best')
    plt.show()

