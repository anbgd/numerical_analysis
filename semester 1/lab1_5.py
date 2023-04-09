import numpy as np
import math

def Householder_transformation(A,i):

	A = np.matrix(A)
	A_n = A.copy()
	E = np.eye(A.shape[0])
	v = np.matrix(find_v(A_n, i))
	vvt = np.matrix(v*(v.transpose()))
	vtv = np.matrix(v.transpose()*v)
	H = E - 2*vvt/vtv
	A_n = H * A

	return(A_n, H)


def QR_decomp(A):

	A_n = A.copy()
	Q = np.matrix(np.eye(len(A)))
	for i in range(A.shape[0]-1):
		A_n, H = Householder_transformation(A_n,i)
		Q *= np.matrix(H)

	return(Q, A_n)


def find_v(A, col):

	v = np.zeros([A.shape[0], 1])
	#print(v)
	sum_v = 0

	for i in range(col,A.shape[0]):
		sum_v += A[i,col]*A[i,col]

	for i in range(A.shape[0]):
		if i != col:
			if i < col:
				v[i,0] = 0
			else:
				v[i,0] = A[i,col]
		else:
			v[i,0] = A[col,col] + np.sign(A[col,col])*math.sqrt(sum_v)

	return(v)


def get_lamda(A,j):

	aii = A[j,j]
	ajj = A[j+1,j+1]
	aij = A[j,j+1]
	aji = A[j+1,j]
	x = (aii+ajj)/2
	D = (-(aii + ajj)*(aii + ajj) + 4*(aii*ajj - aij*aji))
	D = -D
	
	return([x,D])


def sum_column(A,j):

	sum = 0
	for i in range(len(A)-1):
		sum += A[i,j]*A[i,j]

	return(math.sqrt(sum))


def mat_of_lambdas(A):

	s = ((len(A)),2)
	lds = np.matrix(np.zeros(s))
	i = 0
	l = len(A)
	while i < l:
		if (i < l-1): 
			x, D = get_lamda(A, i)
		else:
			D = 1

		if (D < 0): 
			D = -D
			y = math.sqrt(D)/2
			lds[i,0] = x
			lds[i,1] = y
			lds[i+1,0] = x
			lds[i+1,1] = y
			i+=1
		else:
			lds[i,0] = A[i,i]
		i+=1

	return lds


def is_True(old_lmbd,new_lmbd,eps):

	res = abs(old_lmbd[:,0] - new_lmbd[:,0])
	t = [x < eps for x in res]
	if all(t): 
		return False
	else: 
		return True


def find_eigenvalues(Q, R):

	eps = 0.0001
	A = Q*R
	old_lmbd = mat_of_lambdas(A)
	A_new = R*Q
	new_lmbd = mat_of_lambdas(A_new)
	i = 0
	while is_True(old_lmbd,new_lmbd,eps):
		Q_new, R_new = QR_decomp(A_new)
		A_new = R_new*Q_new
		old_lmbd = new_lmbd.copy()
		new_lmbd = mat_of_lambdas(A_new)
		i +=1

	print(i, "итераций")
	#print(new_lmbd)
	k = 0
	while k < new_lmbd.shape[0]:
		if new_lmbd[k,1]!=0:
			print("lambda", k, ":", new_lmbd[k,0], "+", str(new_lmbd[k,1]*1j))
			print("lambda", k+1, ":", new_lmbd[k+1,0], "-", str(new_lmbd[k,1]*1j))
			k+=1
		else:
			print("lambda", k, ":", new_lmbd[k,0])
		k+=1


if __name__ == '__main__':

	#A = np.matrix([[5, -5, -6],[-1, -8, -5], [2, 7, -3]])
	#A = np.matrix([[1, 3, 1],[1, 1, 4], [4, 3, 1]])
	#A = np.matrix([[1, 2, 5],[-8, 0, -6],[7, -9, -7]])
	#A = np.matrix([[17, 6],[6, 8]])
	A = np.matrix([[2.2, 1, 0.5, 2], [1, 1.3, 2, 1], [0.5, 2, 0.5, 1.6], [2, 1, 1.6, 2]])
	print("Матрица A:")
	print(A)
	Q, R = QR_decomp(A)
	print("Матрица Q:")
	print(Q)
	print("Матрица R:")
	print(R)
	print("Проверка A-Q*R:")
	print(A-Q*R)
	print("Поиск СЗ:")
	find_eigenvalues(Q,R)


