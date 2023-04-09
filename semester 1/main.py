import numpy as np
import matplotlib.pyplot as plt
from DSolve import DSolve


def drawResult(X,Y,r,labels):
    labelY = labels[0]
    x = np.arange(0, 3.24, 0.1)
    labelYTrue = labels[1]
    #plt.scatter(X, Y)
    plt.plot(X, Y)
   # plt.plot(x,Y)
    #plt.plot(x, Y)
    plt.legend((labelY, labelYTrue))
    plt.show()

#//////////////////////////////////
res = lambda x: (1+x)*np.exp(x**2)


f = lambda x,y,z:y / (x**2 + y**2)**0.5#4*x*z-(4*x**2-2)*y
g = lambda x,y,z:-x / (x**2 + y**2)**0.5#z
#//////////////////////////////////


#//////////////////////////////////
res1 = lambda x: x + np.exp(-2*x)

fr = lambda x,y,z:z+2*y
f1 = lambda x,y,z:(4*y - 4*x*z)/(2*x+1)
g1 = lambda x,y,z:z
#//////////////////////////////////

dS = DSolve()

#///////////////////////////////////////////////////////////

dS.setStep(0.01)
X,Y = dS.Shoting([f1,g1],[(0,-1),(fr,3)])
dS.setStep(0.02)
X1,Y1 = dS.Shoting([f,g],[(0,-1),(fr,3)])

print("________________________________________")
print("Метод стрельбы\n")
print()
print("Метод - ",Y)
print()
print("Точное решение - ",list(map(res1, X)))
print()
print('Ошибка по Рунге-Ромбергу\n',dS.Runge_Romberg(Y,Y1,4))
print()
print('Ошибка с точным\n',abs(np.reshape(Y,(1,len(Y)))-np.reshape(list(map(res1, X)),(1,len(Y)))))
print("________________________________________")
drawResult(X,Y,res1,('with Shoting method', 'analytical solution'))


dS.setStep(0.01)
X,Y = dS.Ultimate_Difference([f1,g1],[(0,-1),(fr,3)])

print("________________________________________")
print("Конечно-разностный метод\n")
print()
print("Метод - ",Y)
print()
print("Точное решение - ",list(map(res1, X)))
print()
print('Ошибка по Рунге-Ромбергу\n',dS.Runge_Romberg(Y,Y1,4))
print()
print('Ошибка с точным\n',abs(np.reshape(Y,(1,len(Y)))-np.reshape(list(map(res1, X)),(1,len(Y)))))
print("________________________________________")

drawResult(X,Y,res1,('with ultimate difference method', 'analytical solution'))
#///////////////////////////////////////////////////////////
