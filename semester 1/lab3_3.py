import numpy as np
import matplotlib.pyplot as plt


def lu_decompose(mat):
    lu_mat = np.matrix(np.zeros([mat.shape[0],mat.shape[1]]))
    n = mat.shape[0]
    for k in range(n):
        for j in range(k, mat.shape[1]):
            lu_mat[k, j] = mat[k, j] - lu_mat[k, :k] * lu_mat[:k, j]
        for i in range(k+1, n):
            lu_mat[i,k] = (mat[i,k] - lu_mat[i, :k] * lu_mat[:k, k]) / lu_mat[k, k]
    return np.matrix(lu_mat)

def solve_LU(lu_mat, b):
    y = np.matrix(np.zeros([lu_mat.shape[0], 1]))
    for i in range(y.shape[0]):
        y[i, 0] = b[i] - lu_mat[i, :i] * y[:i]
    x = np.matrix(np.zeros([lu_mat.shape[0], 1]))
    for i in range(1, x.shape[0] + 1):
        x[-i, 0] = (y[-i] - lu_mat[-i, -i:] * x[-i:, 0] )/ lu_mat[-i, -i]
    return x


def find_a_vector(n, x, y):
    N = len(x)
    a = np.zeros([n,n]);
    b = np.zeros(n)
    s=""
    for i in range(n):
        for j in range(n):
            tmp_x = 0.0
            tmp_y = 0.0
            for k in range(N):
                tmp_x += pow(x[k], i+j)
                tmp_y += y[k]*pow(x[k], i)
            a[i][j] = tmp_x
            s+="a{}*{} + ".format(j, tmp_x)
        s+="\b\b= {}\n".format(tmp_y)
        b[i] = tmp_y
    A = np.transpose((solve_LU(lu_decompose(a), b)))
    #print(A0)
    #A = np.linalg.solve(a, b)
    #print(A)
    s+="\nF(x) = {}".format(A[0,0])
    for i in range(1,n):
        if(A[0,i]>0):
            s+=" + {}x^{}".format(A[0,i], i)
        else:
            s+=" - {}x^{}".format(-A[0,i], i)
    print(s)
    return A


def get_value(A, x, N):
    tmp=0
    for i in range(N):
        tmp+=A[0,i]*pow(x, i)
    return tmp



def plot_view(N):
    X = np.arange(-0.9,4.5,0.9)
    t = [-0.36892, 0, 0.36892, 0.85408, 1.7856, 6.3138]
    x = [-0.9, 0.0, 0.9, 1.8, 2.7, 3.6]
    #def Y1(x): return(t[])
    if N==2: Y = - 0.19013447619047597 + 1.2462082539682537*X
    if N==3: Y = - 0.4645011428571439- 0.12562507936508224*X + 0.5080864197530874*X*X
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.plot(X,Y, label = 'Приближенный многочлен')
    plt.plot(x,t,'r--', label = ' Табличные данные')
    plt.show()


def main(): 
    print()   
    x = [-0.9, 0.0, 0.9, 1.8, 2.7, 3.6]
    y = [-0.36892, 0, 0.36892, 0.85408, 1.7856, 6.3138]
    N = 3 #степень многочлена
    print("Степень многочлена: {}".format(N-1))
    tmp = 0
    A = find_a_vector(N, x, y)
    for i in range(len(x)):
        tmp += pow(y[i]-get_value(A, x[i], N),2);
    print("\n(y-F)^2 = {}".format(tmp));
    print()
    plot_view(N)



if __name__ == "__main__":
    main()