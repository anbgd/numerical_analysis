import numpy as np
from matplotlib import pyplot as plt
import tkinter as tk
import tkinter.ttk as ttk

a = 1
b = 1

lx = 1
ly = 1
lt = 5
omega = 1.5
nx = 100
ny = 100
nt = 100
eps = 0.001

u0yt = lambda y, t: 0
ulyt = lambda y, t: 0
ux0t = lambda x, t: 0
uxlt = lambda x, t: 0
uxy0 = lambda x, y: x*y
f = lambda x, y, t: -x*y*np.sin(t)

fResult = lambda x, y, t: x*y*np.cos(t)

alpha1 = 1
betta1 = -1
alpha2 = 1
betta2 = -1

hx = lx / nx
hy = ly / ny
ht = lt / nt

tt = 3
yy = 3
U = []


def plotSlice(f, X, y, t):
    plt.subplot(2, 1, 1)
    plt.plot(X, f(X, y, t))
    plt.grid


def showPostProcess(t, y):
    i = 2
    X = np.linspace(0, lx, nx)
    plotSlice(fResult, X, y * hy, ht * t * 2)
    plt.subplot(2, 1, 1)
    plt.plot(X, U[i*t, :, y])


def Error(U, fR):
    i = 2
    X = np.linspace(0, lx, nx)
    Y = np.linspace(0, ly, ny)
    T = list(range(0, int(U.shape[0]/2)))
    plt.subplot(2, 1, 2)
    XX, YY = np.meshgrid(X, Y)
    plt.plot(T, list(map(lambda t: np.max(np.abs(U[i*t] - np.transpose(fR(XX, YY, ht * t * i)))), T)))


def gridFun(iCond, bCond):
    xCond, yCond = bCond
    xlCond, xrCond = xCond
    ylCond, yrCond = yCond
    tlCond, trCond = iCond
    return np.zeros((trCond, xrCond, yrCond))


def parabolicEquation2D(iCond, bCond, method):

    def init():
        for m in range(nx):
            for n in range(ny):
                U[0][m][n] = uxy0(m * hx, n * hy)

    def solve():
        for k in range(0, nt - i, i):
            for n in range(1, ny - 1):
                A1 = getA1(k, n)
                B1 = getB1(k, n)
                U[int(k + i / 2), :, n] = np.linalg.solve(A1, B1)[:, 0]
            for m in range(1, nx - 1):
                A2 = getA2(k, m)
                B2 = getB2(k, m)
                U[k + i, m, :] = np.linalg.solve(A2, B2)[:, 0]
            for n in range(0, ny):
                U[k + i][0][n] = (u0yt(n * hy, (k + i) * ht) - c1 * U[k + i][1][n]) / b1
                U[k + i][nx - 1][n] = (ulyt(n * hy, (k + i) * ht) - aN1 * U[k + i][nx - 2][n]) / bN1

    def FractionSteps():

        aj1 = -a/hx**2
        bj1 = 1/ht/2 + 2*a/hx**2
        cj1 = -a/hx**2

        aj2 = -b/hy**2
        bj2 = 1/ht/2 + 2*b/hy**2
        cj2 = -b/hy**2

        def gA1(k, n):
            aa = np.zeros((nx, nx))
            aa[0][0] = b1
            aa[0][1] = c1
            for j in range(1, nx - 1):
                aa[j][j - 1] = aj1
                aa[j][j] = bj1
                aa[j][j + 1] = cj1
            aa[nx - 1][nx - 2] = aN1
            aa[nx - 1][nx - 1] = bN1
            return aa

        def gB1(k, n):
            bb = np.zeros((nx, 1))
            bb[0][0] = u0yt(n * hy, ht * (k + i))
            bb[nx - 1][0] = ulyt(n * hy, ht * (k + i))
            for m in range(1, nx - 1):
                bb[m][0] = f(m * hx, n * hy, k * ht)/2 + U[k][m][n] / ht / 2
            return bb

        def gA2(k, m):
            aa = np.zeros((ny, ny))
            aa[0][0] = b2
            aa[0][1] = c2
            for j in range(1, ny - 1):
                aa[j][j - 1] = aj2
                aa[j][j] = bj2
                aa[j][j + 1] = cj2
            aa[ny - 1][ny - 2] = aN2
            aa[ny - 1][ny - 1] = bN2
            return aa

        def gB2(k, m):
            bb = np.zeros((ny, 1))
            bb[0][0] = ux0t(m * hx, ht * (k + i))
            bb[ny - 1][0] = uxlt(m * hx, ht * (k + i))
            for n in range(1, ny - 1):
                bb[n][0] = f(m * hx, n * hy, (k + i) * ht)/2 + U[int(k + i/2)][m][n] / ht / 2
            return bb
        return (gA1, gB1), (gA2, gB2)

    def AlternatingDirection():

        aj1 = -a / hx ** 2
        bj1 = 1 / ht + 2 * a / hx ** 2
        cj1 = -a / hx ** 2

        aj2 = -b / hy ** 2
        bj2 = 1 / ht + 2 * b / hy ** 2
        cj2 = -b / hy ** 2

        def gA1(k, n):
            aa = np.zeros((nx, nx))
            aa[0][0] = b1
            aa[0][1] = c1
            for j in range(1, nx - 1):
                aa[j][j - 1] = aj1
                aa[j][j] = bj1
                aa[j][j + 1] = cj1
            aa[nx - 1][nx - 2] = aN1
            aa[nx - 1][nx - 1] = bN1
            return aa

        def gB1(k, n):
            bb = np.zeros((nx, 1))
            bb[0][0] = u0yt(n * hy, ht * (k + i))
            bb[nx - 1][0] = ulyt(n * hy, ht * (k + i))
            for m in range(1, nx - 1):
                bb[m][0] = f(m * hx, n * hy, int(k + i/2) * ht) + (1 / ht - 2*b/ hy**2) * U[k][m][n] + (b/hy**2) * (U[k][m][n + 1] + U[k][m][n - 1])
            return bb

        def gA2(k, m):
            aa = np.zeros((ny, ny))
            aa[0][0] = b2
            aa[0][1] = c2
            for j in range(1, ny - 1):
                aa[j][j - 1] = aj2
                aa[j][j] = bj2
                aa[j][j + 1] = cj2
            aa[ny - 1][ny - 2] = aN2
            aa[ny - 1][ny - 1] = bN2
            return aa

        def gB2(k, m):
            bb = np.zeros((ny, 1))
            bb[0][0] = ux0t(m * hx, ht * (k + i))
            bb[ny - 1][0] = uxlt(m * hx, ht * (k + i))
            for n in range(1, ny - 1):
                bb[n][0] = f(m * hx, n * hy, int(k + i/2) * ht) + (1 / ht - 2*a/ hx**2) * U[int(k + i/2)][m][n] + (a/hx**2) * (U[int(k + i/2)][m + 1][n] + U[int(k + i/2)][m - 1][n])
            return bb

        return (gA1, gB1), (gA2, gB2)
    if method == 1:
        (getA1, getB1), (getA2, getB2) = AlternatingDirection()
    elif method == 2:
        (getA1, getB1), (getA2, getB2) = FractionSteps()
    else:
        pass
    b1 = 1
    c1 = 0
    b2 = 1
    c2 = 0

    aN1 = -betta1 / hx
    bN1 = alpha1 + betta1 / hx
    aN2 = -betta2 / hy
    bN2 = alpha2 + betta2 / hy

    uxy0 = iCond
    xCond, yCond = bCond
    u0yt, ulyt = xCond
    ux0t, uxlt = yCond
    U = gridFun((0, nt), ((0, nx), (0, ny)))
    i = 2
    init()
    solve()
    return U


def solver():
    global a, b, nx, ny, nt, hx, hy, ht, U, tt, yy
    i = 2
    a = float(entryA.get())
    b = float(entryB.get())
    nx = int(scaleX.get())
    ny = int(scaleY.get())
    nt = 2 * (int(scaleT.get()))
    tt = int(t0.get())
    yy = int(y0.get())
    hx = lx / nx
    hy = ly / ny
    ht = lt / nt
    nx += 1
    ny += 1
    nt += i
    U = []

    method = combobox1.get()
    if method == 'Переменных направлений':
        U = parabolicEquation2D(uxy0, ((ux0t, uxlt), (u0yt, ulyt)), 1)
    elif method == 'Дробных шагов':
        U = parabolicEquation2D(uxy0, ((ux0t, uxlt), (u0yt, ulyt)), 2)
    else:
        pass
    showPostProcess(tt, yy)
    Error(U, fResult)
    plt.show()


root = tk.Tk()
root.title("Лабораторная работа №4")
frame = ttk.Frame(root)
frame.grid()
combobox1 = ttk.Combobox(frame, values=["Переменных направлений", 'Дробных шагов'], height=3, width=50)
button = ttk.Button(root, text="Решить", command=solver)
image = tk.PhotoImage(file="lab4.PNG")
lab = ttk.Label(frame, image=image)
lab0 = ttk.Label(frame, text="Выберите метод")
labgrid = ttk.Label(frame, text="Выберите параметры сетки")
labtask = ttk.Label(frame, text="Выберите параметры задачи:\n\ta,b")
sliceTaskT = ttk.Label(frame, text="Выберите срез по T")
sliceTaskY = ttk.Label(frame, text="Выберите срез по Y")
scaleX = tk.Scale(frame, orient=tk.HORIZONTAL, length=200, from_=0, tickinterval=20, resolution=1, to=100)
scaleY = tk.Scale(frame, orient=tk.HORIZONTAL, length=200, from_=0, tickinterval=20, resolution=1, to=100)
scaleT = tk.Scale(frame, orient=tk.HORIZONTAL, length=200, from_=0, tickinterval=20, resolution=1, to=100)
scaleX.set(10)
scaleY.set(10)
scaleT.set(10)

entryA = tk.Entry(frame, width=10, bd=10)
entryB = tk.Entry(frame, width=10, bd=10)

t0 = tk.Entry(frame, width=10, bd=10)
y0 = tk.Entry(frame, width=10, bd=10)
t0.insert(0, 0)
y0.insert(0, 1)

entryA.insert(0, a)
entryB.insert(0, b)

combobox1.set("Переменных направлений")

timeSlice = ttk.Label(frame, text="t0\t=", font="arial 20")
xSlice = ttk.Label(frame, text="x0\t=", font="arial 20")
labtaskA = ttk.Label(frame, text="a\t=", font="arial 20")
labtaskB = ttk.Label(frame, text="b\t=", font="arial 20")
labgridX = ttk.Label(frame, text="nx\t=", font="arial 20")
labgridY = ttk.Label(frame, text="ny\t=", font="arial 20")
labgridT = ttk.Label(frame, text="nt\t=", font="arial 20")
labelgrid0 = ttk.Label(frame, background='#cc0')

lab.grid(row=0, column=0, columnspan=3)
labtask.grid(row=1, column=0)
labtaskA.grid(row=1, column=1)
labtaskB.grid(row=2, column=1)

entryA.grid(row=1, column=2)
entryB.grid(row=2, column=2)

labgrid.grid(row=3,column=0)
labgridX.grid(row=4, column=1)
labgridY.grid(row=5, column=1)
labgridT.grid(row=6, column=1)

scaleX.grid(row=4, column=2)
scaleY.grid(row=5, column=2)
scaleT.grid(row=6, column=2)

sliceTaskT.grid(row=7, column=0)
sliceTaskY.grid(row=8, column=0)
t0.grid(row=7, column=1)
y0.grid(row=8, column=1)
lab0.grid(row=9, column=0)
combobox1.grid(row=9, column=1)
button.grid(row=10, column=0)

style = ttk.Style()
style.configure("TLabel", padding=3, background='#bb2', font="arial 12", foreground="black")
style.configure("TFrame", background='#CC0')
style.configure("TButton", width=20, height=5, font="arial 20", foreground='red')
root.mainloop()
