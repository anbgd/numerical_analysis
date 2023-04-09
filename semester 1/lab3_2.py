import numpy as np

def get_PQ(mat, answ):
    P1 = -mat[0,1]/mat[0,0]
    Q1 = answ[0]/mat[0,0]
    size_mat = mat[:,0].size
    i = 1
    P = [P1]
    Q = [Q1]
    while i < size_mat-1:
        P.append(-mat[i,i+1]/(mat[i,i]+mat[i,i-1]*P[i-1]))
        Q.append((answ[i]-mat[i,i-1]*Q[i-1])/(mat[i,i]+mat[i,i-1]*P[i-1]))
        i+=1
    P.append(0)
    Q.append((answ[i]-mat[i,i-1]*Q[i-1])/(mat[i,i]+mat[i,i-1]*P[i-1]))
    return P, Q

def solve_run(mat, b):
    P, Q = get_PQ(mat, b)
    size_mat = mat[:,0].size
    X = np.zeros(size_mat)
    X[size_mat-1] = Q[size_mat-1]
    i = size_mat-2
    while i > -1:
        X[i] = X[i+1]*P[i]+Q[i]
        i-=1
    return X

def h(i, x):
    return x[i]-x[i-1]

def get_c_value(x, y):
    n = len(x)
    X = np.zeros(n-1)
    A = np.zeros([n-2,n-2])
    B = np.zeros(n-2)
    A[0][0] = 2*(h(1, x)+h(2, x))
    A[0][1] = h(2, x)
    B[0] = 3*(((y[2] - y[1]) / h(2, x)) - ((y[1] - y[0]) / h(1, x)))
    j=0
    for i in range(3,n-1):
        A[i-2][j] = h(i-1,x)
        A[i-2][j+1] = 2*(h(i-1, x)+h(i, x))
        A[i-2][j+2] = h(i,x)
        B[i-2] = 3*(((y[i] - y[i-1]) / h(i, x)) - ((y[i-1] - y[i-2]) / h(i-1, x)))
        j+=1
    A[n-3][j] = h(n-2, x)
    A[n-3][j+1] = 2*(h(n-2, x)+h(n-1, x))
    B[n-3] = 3*(((y[n-1] - y[n-2]) / h(n-1, x)) - ((y[n-2] - y[n-3]) / h(n-2, x)))
    X[1:] = solve_run(A,B)
    #X[1:] = np.linalg.solve(A,B)
    return X


def get_a_b_d(x, y, c):
    n = len(x)-1
    a = np.zeros(n+1)
    b = np.zeros(n+1)
    d = np.zeros(n+1)
    for i in range(1, n):
        a[i-1] =  y[i-1]
        b[i-1] = ((y[i] - y[i-1])/h(i,x)) - ((h(i,x)*(c[i] + 2*c[i-1]))/3)
        d[i-1] = (c[i] - c[i-1]) / (3*h(i,x))
    a[n-1] = y[n-1]
    b[n-1] = ( ((y[n] - y[n-1]) / h(n-1,x)) - 2*h(n-1,x)*c[n-1]/3 )
    d[n-1] = (-c[n-1]/(3*h(n-1, x)))
    return a,b,d


def get_spline(X, a, b, c, d, x):
    n = len(x)
    j = 0;
    for k in range(int(x[0]),int(x[n-1]),1):
        if ((k<=X) and (k>=X-1)):
            break;
        j+=1
    for i in range(n-1):
        s = "Сплайн: {} + {}(x - {}) + {}(x - {})^2 + {}(x - {})^3".format(a[i], b[i], x[i], c[i], x[i], d[i], x[i])
        if (i==j):
            s+=" !!!!!"
        print(s)
    t = X-x[j];
    return a[j] + b[j]*t + c[j]*t*t + d[j]*t*t*t


def main():
    print()
    x = [0, 0.9, 1.8, 2.7, 3.6]
    print("x")
    print(x) 
    y = [0, 0.36892, 0.85408, 1.7856, 6.3138]
    print("y")
    print(y)
    X = 1.5
    print("X =", X)
    c = get_c_value(x,y)
    a, b, d = get_a_b_d(x, y, c)
    print("f({}) = {}".format(X,get_spline(X, a,b,c,d,x)))
    print()
    
    
if __name__ == "__main__":
    main()






