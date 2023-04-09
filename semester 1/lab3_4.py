import numpy as np


def first_derivative(i, x, y):
    return (y[i+1] - y[i]) / (x[i+1] - x[i])


def second_derivative(i, x, y):
    a = (y[i+2] - y[i+1]) / (x[i+2] - x[i+1])#правосторонняя
    b = (y[i+1] - y[i]) / (x[i+1] - x[i])#левосторонняя
    c = x[i+2] - x[i]
    return 2 * (a - b) /c


def main():
    print()
    x = [1.0, 1.5, 2.0, 2.5, 3.0]#my
    print('x =',x)
    #x = [0.0, 0.1, 0.2, 0.3, 0.4]
    y = [0, 0.40547, 0.69315, 0.91629, 1.0986]#my
    print('y =',y)
    #y = [1.0, 1.1052, 1.2214, 1.3499, 1.4918]
    x0 = 2.0#my
    #x0 = 0.2
    for i in range(len(x)):
        if ((x[i]<=x0) and (x0 <=x[i+1])):
            print("y' = {}".format(first_derivative(i, x, y)))
            print("y'' = {}".format(second_derivative(i, x, y)))
            break
    print()


if __name__ == "__main__":
    main()