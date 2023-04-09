import numpy as np
from matplotlib import pyplot as plt
import tkinter as tk
import tkinter.ttk as ttk

import warnings
warnings.filterwarnings("ignore")

a = -2
b = -2
c = -4
lx = np.pi/2
ly = np.pi/2
omega = 1.5
K = 100
N = 100
epsil = 0.0001

u0y = lambda y: np.exp(-y)*np.cos(y)
uly = lambda y: 0
ux0 = lambda x: np.exp(-x)*np.cos(x)
uxl = lambda x: 0
fResult = lambda x, y: np.exp(-x-y)*np.cos(x)*np.cos(y)

hx = lx / K
hy = ly / N
U = []


def gridFun(xCond, tCond):
    xlCond,xrCond = xCond
    tlCond, trCond = tCond
    return np.zeros((xrCond, trCond))


def plotSlice(f, X, t):
    plt.subplot(2, 1, 1)
    plt.plot(X, f(X, t))
    plt.grid


def showPostProcess(y, err):
    X = np.linspace(0, lx, K)
    plotSlice(fResult, X, hy * y)
    plt.subplot(2, 1, 1)
    plt.plot(X, U[:, y])
    plt.subplot(2, 1, 2)
    plt.plot(err)


def ellipticEquation(xConds, yConds, method):
    def interpolation():
        for i in range(1, K - 1):
            for j in range(1, N - 1):
                alpha = (j * hy) / ly
                U[i][j] = ux0(i * hx)*(1 - alpha) + uxl(i * hx) * alpha
        return U


    def Libman(U, epsil):
        n = 0
        errors = []
        while True:
            n += 1
            Uold = U.copy()
            for i in range(1, K - 1):
                for j in range(1, N - 1):
                    U[i][j] = delta * ((hhx + ahx) * Uold[i - 1][j] +
                                       (hhx - ahx) * Uold[i + 1][j] +
                                       (hhy + bhy) * Uold[i][j - 1] +
                                       (hhy - bhy) * Uold[i][j + 1])
            err = np.max(np.abs(Uold - U))
            errors.append(err)
            if (err < epsil):
                break
        print("Libman - ", n)
        return U, errors

    def Seidel(U, epsil):
        n = 0
        errors = []
        while True:
            n += 1
            Uold = U.copy()
            for i in range(1, K - 1):
                for j in range(1, N - 1):
                    U[i][j] = delta * ((hhx + ahx) * U[i - 1][j] +
                                       (hhx - ahx) * U[i + 1][j] +
                                       (hhy + bhy) * U[i][j - 1] +
                                       (hhy - bhy) * U[i][j + 1])
            err = np.max(np.abs(Uold - U))
            errors.append(err)
            if (err < epsil):
                break
        print("Seidel - ", n)
        return U, errors

    def upperRelaxation(U, epsil, omega):
        n = 0
        good = False
        errors = []
        while (n < 1000):
            n += 1
            Uold = U.copy()
            for i in range(1, K - 1):
                for j in range(1, N - 1):
                    U[i][j] = U[i][j] + omega * (delta * ((hhx + ahx) * U[i - 1][j] +
                                       (hhx - ahx) * Uold[i + 1][j] +
                                       (hhy + bhy) * U[i][j - 1] +
                                       (hhy - bhy) * Uold[i][j + 1]) - U[i][j])
            err = np.max(np.abs(Uold - U))
            errors.append(err)
            if (err < epsil):
                good = True
                break
        if (not good):
            print("Расходится!!!")
        print("upperRelaxation - ", n)
        return U, errors

    ux0, uxl = xConds
    u0y, uly = yConds
    U = gridFun((0, K), (0, N))

    for i in range(0, K):
        U[i][0] = ux0(hx * i)
        U[i][N-1] = uxl(hx * i)

    for j in range(0, N):
        U[0][j] = u0y(hy * j)
        U[K-1][j] = uly(hy * j)

    delta = 1/(2/hx**2 + 2/hy**2 + c)
    hhx = 1/hx**2
    ahx = a/2/hx
    hhy = 1/hy**2
    bhy = b/2/hy
    err = []
    U = interpolation()
    if method == 1:
        U, err = Libman(U, epsil)
    elif method == 2:
        U, err = Seidel(U, epsil)
    elif method == 3:
        U, err = upperRelaxation(U, epsil, omega)
    else:
        pass
    return U, err


def solver():
    global a, b, c, N, K, h, hx, hy, U, omega,epsil
    a = float(entrya.get())
    b = float(entryb.get())
    c = float(entryc.get())
    omega = float(scaleOm.get())
    K = int(scaleX.get())
    N = int(scaleY.get())
    tt = int(t0.get())
    hy = ly / N
    hx = lx / K
    U = []
    err = []
    epsil = float(entryep.get())

    method = combobox1.get()
    if method == 'Либмана':
        U, err = ellipticEquation((ux0, uxl), (u0y, uly), 1)
    elif method == 'Зейделя':
        U, err = ellipticEquation((ux0, uxl), (u0y, uly), 2)
    elif method == 'с Верхней релаксацией':
        U, err = ellipticEquation((ux0, uxl), (u0y, uly), 3)
    else:
        pass
    showPostProcess(tt, err)
    plt.show()


root = tk.Tk()
root.title("Лабораторная работа №3")
frame = ttk.Frame(root)
frame.grid()
combobox1 = ttk.Combobox(frame, values=["Либмана", "Зейделя", "с Верхней релаксацией"], height=3, width=50)
button = ttk.Button(root, text="Решить", command=solver)
image = tk.PhotoImage(file="nm3.PNG")
lab = ttk.Label(frame, image=image)
lab0 = ttk.Label(frame, text="Выберите метод")
labgrid = ttk.Label(frame, text="Выберите параметры сетки")
labtask = ttk.Label(frame, text="Выберите параметры задачи:\n\ta,b,c")
metPar = ttk.Label(frame, text="Выберите параметр для 3-го метода")
sliceTask = ttk.Label(frame, text="Выберите срез по Y")
scaleX = tk.Scale(frame, orient=tk.HORIZONTAL, length=200, from_=0, tickinterval=20, resolution=1, to=100)
scaleY = tk.Scale(frame, orient=tk.HORIZONTAL, length=200, from_=0, tickinterval=20, resolution=1, to=100)
scaleOm = tk.Scale(frame, orient=tk.HORIZONTAL, length=200, from_=0, tickinterval=0.4, resolution=0.1, to=2)
scaleX.set(60)
scaleY.set(60)
scaleOm.set(1.5)

entrya = tk.Entry(frame, width=10, bd=10)
entryb = tk.Entry(frame, width=10, bd=10)
entryc = tk.Entry(frame, width=10, bd=10)
entryep = tk.Entry(frame, width=10, bd=10)

t0 = tk.Entry(frame, width=10, bd=10)
t0.insert(0, 3)

entrya.insert(0, a)
entryb.insert(0, b)
entryc.insert(0, c)
entryep.insert(0, epsil)

combobox1.set("Либмана")

timeSlice = ttk.Label(frame, text="t0\t=", font="arial 20")
labtaska = ttk.Label(frame, text="a\t=", font="arial 20")
labtaskb = ttk.Label(frame, text="b\t=", font="arial 20")
labtaskc = ttk.Label(frame, text="c\t=", font="arial 20")
labtaskep = ttk.Label(frame, text="epsil\t=", font="arial 20")
labgridX = ttk.Label(frame, text="Nx\t=", font="arial 20")
labgridY = ttk.Label(frame, text="Ny\t=", font="arial 20")
metPar1 = ttk.Label(frame, text="omega\t=", font="arial 20")
labelgrid0 = ttk.Label(frame, background='#cc0')

lab.grid(row=0, column=0, columnspan=3)
labtask.grid(row=1, column=0)
labtaska.grid(row=1, column=1)
labtaskb.grid(row=2, column=1)
labtaskc.grid(row=3, column=1)
labtaskep.grid(row=4, column=1)

entrya.grid(row=1, column=2)
entryb.grid(row=2, column=2)
entryc.grid(row=3, column=2)
entryep.grid(row=4, column=2)

labgrid.grid(row=5,column=0)
labgridX.grid(row=5, column=1)
labgridY.grid(row=6, column=1)
metPar.grid(row=7, column=0)
metPar1.grid(row=7, column=1)

scaleX.grid(row=5, column=2)
scaleY.grid(row=6, column=2)
scaleOm.grid(row=7, column=2)

sliceTask.grid(row=8, column=0)
t0.grid(row=8, column=1)
lab0.grid(row=9, column=0)
combobox1.grid(row=9, column=1)
button.grid(row=10, column=0)

style = ttk.Style()
style.configure("TLabel", padding=3, background='#bb2', font="arial 12", foreground="black")
style.configure("TFrame", background='#CC0')
style.configure("TButton", width=20, height=5, font="arial 20", foreground='red')
root.mainloop()
