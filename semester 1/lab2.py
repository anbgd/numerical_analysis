import numpy as np
import matplotlib.pyplot as plt

#Выраженная функция
def fi_func(x):
    return np.log(np.sqrt(1-x*x)+0.1)

#Исходня функция
def f_func(x):
    return np.sqrt(1-x*x)-np.exp(x)+0.1

#Производня исходной функции
def f_derivative(x):
    return (-((x/(np.sqrt(1-x*x)))+np.exp(x)))

def fi_system(x):
    return [np.sqrt(16-(x[1]-2)*(x[1]-2))+2, 64/(x[0]*x[0]+16)]#[np.sqrt((64/x[1])-16), np.sqrt(16-(x[0]-2)*(x[0]-2))+2]

#Исходная система уравнений
def f_system(x):
    return [(x[0]*x[0]+16)*x[1]-64, (x[0]-2)*(x[0]-2)+(x[1]-2)*(x[1]-2)-16]

#Производная исходной системы по х1
def f_system_dx1(x):
    return [2*x[0]*x[1], 2*(x[0]-2)]

#Производная исходной системы по х2
def f_system_dx2(x):
    return [x[0]*x[0]+16, 2*(x[1]-2)]


#Метод простых итераций
def simple_iterations(fi_func, x0, epsilon, q):
    x = x0
    k=0
    cont = True
    while (cont):
        print("k={} x={} fi_func(x)={}".format(k, x, fi_func(x)))
        prev_x = x
        x = fi_func(x)
        k=k+1
        if ((q/(1-q)) * np.absolute(x - prev_x) < epsilon or k>1000):
            cont = False
    return x;

#Метод Ньютона
def newton_mathod(f_func, f_derivative, x0, epsilon):
    x = x0
    k=0
    cont = True
    while(cont):
        print("k={} x={} f_func(x)={} f_derivative(x)={} f_func(x)/f_derivative(x)={}".format(k, x, f_func(x), f_derivative(x), f_func(x)/f_derivative(x)))
        prev_x = x
        x -= f_func(x)/f_derivative(x)
        k=k+1
        if ( np.absolute(x-prev_x) <= epsilon or k>500):
            cont = False
    return x

#Метод итераций для системы 
def simple_iterations_fs(fi, x0, epsilon, q):
    k = 0
    cont = True
    while(cont):
        print("x1={} x2={} f1(x1,x2)={} f2(x1,x2)={}".format(x0[0],x0[1], fi(x0)[0], fi(x0)[1]))
        prevu_x = x0.copy()
        x0 = fi(x0)
        k=k+1
        if (np.sqrt(pow(x0[0]-prevu_x[0],2)+pow(x0[1]-prevu_x[1],2))<epsilon or k>150):
            cont = False
    print("{} итераций сделано".format(k));
    return x0

#Метод Ньютона для системы
def newton_method_fs(f, fdx1, fdx2, x0, epsilon):
    k = 0
    cont = True
    while(cont):
        prevu_x = x0.copy()
        fm = f(prevu_x)
        fdx1m = fdx1(prevu_x)
        fdx2m = fdx2(prevu_x)
        print("x1={} x2={} f1(x1,x2)={} f2(x1,x2)={} f1dx1(x1,x2)={} f1dx2(x1,x2)={} f2dx1(x1,x2)={} f2dx2(x1,x2)={}".format(x0[0], x0[1], fm[0], fm[1], fdx1m[0], fdx2m[0], fdx1m[1], fdx2m[1]))
        x0[0] -= (fm[0]*fdx2m[1]-fm[1]*fdx2m[0])/(fdx1m[0]*fdx2m[1]-fdx1m[1]*fdx2m[0])
        x0[1] -= (fm[1]*fdx1m[0]-fdx1m[1]*fm[0])/(fdx1m[0]*fdx2m[1]-fdx1m[1]*fdx2m[0])
        k=k+1
        if(np.sqrt(pow(x0[0]-prevu_x[0],2)+pow(x0[1]-prevu_x[1],2))<epsilon or k>150):
            cont = False
    print("{} итераций сделано".format(k));
    return x0


def draw_fs():
    X = np.arange(3,6.0000001,0.001)#-2,6
    Y1 = np.sqrt(16-(X-2)*(X-2))+2
    Y2 = -np.sqrt(16-(X-2)*(X-2))+2
    Y3 = 64/(X*X+16)
    plt.grid()
    plt.plot(X,Y2)
    plt.plot(X,Y3)
    plt.plot(X,Y1)
    plt.show()



def draw():
    X = np.arange(0,1,0.001)#-1,1
    Y1 = np.sqrt(1-X*X)+0.1
    Y2 = np.exp(X)
    plt.grid()
    plt.plot(X,Y1)
    plt.plot(X,Y2)
    plt.show()


if __name__ == '__main__':


    eps = float(0.0001)#input("\nenter precision: ")
    print("Заданный eps={}\n".format(eps))
    print("\nУравнения\n")

    print("Решение методом простых итераций")
    print("Итого\nx={}\n".format(simple_iterations(fi_func, 0, eps, 0.9)))
    print("Решение методом Ньютона")
    print("Итого\nx={}".format(newton_mathod(f_func, f_derivative, 0, eps)))
    draw()

    print("\nСистема уравнений\n")

    print("Решение методом простых итераций")
    res = simple_iterations_fs(fi_system, [6, 1.3], eps, 0.9)
    print("Итого\nx1={} x2={}\n".format(res[0], res[1]))
    print("Решение методом Ньютона")
    res = newton_method_fs(f_system, f_system_dx1, f_system_dx2, [6, 1.3], eps)
    print("Итого\nx1={} x2={}".format(res[0], res[1]))
    draw_fs()

