import numpy as np
from numpy.linalg import qr, solve

from scipy.linalg import norm
import scipy.stats as stats
import scipy.sparse as sparse
from scipy.sparse.linalg import eigsh
import lab1_5
from pytictoc import TicToc


#https://stackoverflow.com/questions/26895011/python-how-to-use-python-to-generate-a-random-sparse-symmetric-matrix
#создаем симметричную разреженную эрмитову матрицу
def SymmMat(n, density):
    #np.random.seed((0,0))     #откомментить, чтобы всегда генерировать одинаковые числа
    rvs = stats.norm().rvs
    X = sparse.random(n, n, density=density, data_rvs=rvs)
    #X.data[:] = 1              #замена ненулевых эл-тов единицами
    upper_X = sparse.triu(X) 
    result = upper_X + upper_X.T - sparse.diags(X.diagonal())
    return result


def special_qr(A, eps=1E-14, maxiter=5000):
    m = A.shape[0]
    Q = np.identity(m)
    residual = 10
    lprev = np.ones(m)
    ctr = 0
    while norm(residual) > eps:
        Q,R = qr(A@Q)
        lam = np.diagonal(Q.T @ A @ Q) #Матрица коэффициентов Рэлея
        #print(lam)
        residual = norm(lprev - np.sort(lam))
        lprev = np.sort(lam)
        ctr += 1
        if ctr == maxiter: break
    #print(ctr)
    print("Матрица собственных векторов")
    for i in range(1,len(Q)):
    	Q[:,i] *= -1
    print(Q[0,:])
        
    return(lam)


def LanczosTri(A):
    #триагонализация методом Ланцоша
    
    #проверка симметричности
    #if((A.transpose() != A).any()):
    #    print("WARNING: Введенная матрица не симметрична")
    n = A.shape[0]
    x = np.ones(n)                      #инициализационный вектор
    V = np.zeros((n,1))                 #триагонализирующая матрица

    #начало
    q = x/norm(x)
    V[:,0] = q
    r = A @ q
    a1 = q.T @ r
    r = r - a1*q
    b1 = norm(r)
    ctr = 0
    #print("a1 = %.12f, b1 = %.12f"%(a1,b1))
    for j in range(2,n+1):
        v = q
        q = r/b1
        r = A @ q - b1*v
        a1 = q.T @ r
        r -= a1*q
        b1 = norm(r)
        
        #вставляем новый столбец в конец V
        V = np.hstack((V,np.reshape(q,(n,1))))

        #Реортоганализируем все предыдущие v
        V = qr(V)[0]

        ctr+=1
        
        if b1 == 0: 
            print("WARNING: Метод Ланцоша приостановлен из-за b1 = 0")
            return V #нужно реотртонормализировать
        
        #print(np.trace(V.T@V)/j)
    #проверка ортонормальности V
    #print("|V.T@V - I| = ")
    #print(np.abs((V.T@V)-np.eye(n)))
    #if((V.T@V != np.eye(n)).any()):
    #    print("WARNING: V.T @ V != I: Ортонормальность преобразоваия потеряна")
        
    
    T = V.T @ A @ V
    
    return T


def main():
    #создаем трехдиагональную матрицу
    n = 75                                  #размер матрицы (nxn)
    A = SymmMat(n,density=1)              #Эрмитова, разреженная
    print("Матрица А")
    print(A)
    #print(type(A))
    
    #проверка на симметричность матрицы
    assert (A - A.T).nnz == 0
    
    #просто формат вывода
    np.set_printoptions(formatter={'float_kind':'{:f}'.format})
    
    #трансформация в трехдиагональную
    print("Преобразованная матрица - T")
    T = LanczosTri(A)

    print(T)
    
    print()
    print("Моя реализация QR(T)")
    t1 = TicToc()
    t1.tic()
    Q,R = lab1_5.QR_decomp(T)
    lab1_5.find_eigenvalues(Q,R)
    t1.toc()
    print()
	
    print("Моя реализация QR(A)")
    t5 = TicToc()
    t5.tic()
    Q,R = lab1_5.QR_decomp(A.toarray())
    lab1_5.find_eigenvalues(Q,R)
    t5.toc()
    print()
    
    print("Далее выводится только минимальное собственное значение и соответствующий собственный вектор")
    #СЗ с помощью QR
    print("С помощью встроенного QR")
    t2 = TicToc()
    t2.tic()
    lam = special_qr(T)
    t2.toc()
    print("Собственные значения(T) (special_qr): ",lam, "\n")#np.sort(lam)[:-1][0]
    
    #С использованием встроенных функций
    print("С использованием встроенных функций от матрицы T")
    t3 = TicToc()
    t3.tic()
    e_gs_T, gs_T = np.linalg.eigh(T)
    t3.toc()
    #e_gs_A = special_qr(A,maxiter=1000)
    print("Собственные значения(T) (np.eigsh): ", e_gs_T[0])
    print("Собственные векторы(T) (np.eigsh): ", gs_T[0, :], "\n")#[0, :]
    #print("Eigs(A): ",np.sort(e_gs_A[:-1]))
    

    print("С использованием встроенных функций от матрицы A")
    t4 = TicToc()
    t4.tic()
    e_gs_A, gs_A = eigsh(A,k=n-1,which='SA',maxiter=1000)
    t4.toc()
    #e_gs_A = special_qr(A,maxiter=1000)
    print("Собственные значения(A) (np.eigsh): ", e_gs_A[0])
    print("Собственные векторы(T) (np.eigsh): ", gs_A[:,0], "\n")
    #print("Eigs(A): ",np.sort(e_gs_A[:-1]))



    
    
if __name__ == '__main__':
    main()