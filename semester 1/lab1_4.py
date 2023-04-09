import numpy as np
import math


def get_angle(mat, i,j):
    if mat[i,i] == mat[j,j]:
        return math.PI/4
    else:
        return 0.5 * math.atan(2*mat[i,j]/(mat[i,i] -mat[j,j]))
    
def get_max(mat):
    index = [0,1]
    mat_size = mat[:,0].size
    mat_tmp = np.abs(mat.copy())
    max_elem = mat_tmp[index[0],index[1]]
    for i in np.arange(0,mat_size):
        for j in np.arange(0,mat_size):
            if i==j:
                continue
            else:
                if max_elem < mat_tmp[i,j]:
                    max_elem = mat_tmp[i,j]
                    index = [i,j]
                else:
                    pass
    
    return index

def norm_matrix_non_diag(matrix):
    norm = 0
    mat_size = matrix[:,0].size
    for i in np.arange(0,mat_size):
        for j in np.arange(0,mat_size):
            if i==j:
                continue
            else:
                norm += matrix[i,j]*matrix[i,j]
    return np.sqrt(norm)


def transpose(matr):
    res=[]
    n=len(matr)
    m=len(matr[0])
    for j in range(m):
        tmp=[]
        for i in range(n):
            tmp=tmp+[matr[i][j]]
        res=res+[tmp]
    return np.array(res)



if __name__ == '__main__':

	epsilon = 0.00001
	print()
	print("epsilon = ", epsilon)
	mat = [[5, 5, 3], [5, -4, 1], [3, 1, 2]]
	mat = np.asarray(mat)
	print('Исходная матрица: ')
	print(mat)
	print('\n')
	mat_test = mat.copy()
	mat_size = mat[:,0].size
	X = np.eye(mat_size)
	k = 0
	while True:
	    max_index = get_max(mat)
	    angle = get_angle(mat, max_index[0], max_index[1])
	    U = np.eye(mat_size).copy()
	    U[max_index[0],max_index[0]] =  math.cos(angle)
	    U[max_index[0],max_index[1]] = -math.sin(angle)
	    U[max_index[1],max_index[0]] =  math.sin(angle)
	    U[max_index[1],max_index[1]] =  math.cos(angle)
	    X = np.dot(X.copy(),U.copy())
	    mat = np.dot(transpose(U.copy()),mat.copy())
	    mat = np.dot(mat.copy(), U.copy())
	    p = norm_matrix_non_diag(mat)
	    k += 1
	    if epsilon > p or angle == 0 or k > 100:
	        break
	    
	#Ax=lambda*x x-СВ, lambda-СЗ
	X[:,0] = -1*X[:,0]
	print("Число итераций: ", k)
	print("Собственные значения:")
	print(np.diagonal(mat))
	print('\n')
	print("Собственные векторы:")
	print(X)
	print('\n')
	print("Проверка A*x-lambda*x=0:")
	#a, b = np.linalg.eig(mat_test)
	#print(a,'\n',b)
	mat_test = np.matrix(mat_test)
	X = np.matrix(X)
	mat = np.matrix(mat)
	for i in range(len(mat)):
		print("A*x",i+1)		
		print(mat_test * X[:, i])
		print("lambda",i+1,"*x",i+1)
		print(mat[i,i]* X[:, i])





