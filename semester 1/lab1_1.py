import numpy as np

def lu_decompose(mat):

	lu_mat = np.matrix(np.zeros([mat.shape[0],mat.shape[1]]))
	n = mat.shape[0]

	for k in range(n):
		for j in range(k, mat.shape[1]):
			lu_mat[k, j] = mat[k, j] - lu_mat[k, :k] * lu_mat[:k, j]
		for i in range(k+1, n):
			lu_mat[i,k] = (mat[i,k] - lu_mat[i, :k] * lu_mat[:k, k]) / lu_mat[k, k]

	return np.matrix(lu_mat)

def get_L(lu_mat):

	L = lu_mat.copy()

	for i in range(L.shape[0]):

		L[i, i] = 1
		L[i, i+1:] = 0

	#if np.allclose(L[:, L.shape[0]], np.zeros([L.shape[0],1])):
		#L = np.delete(L, L.shape[1]-1, 1)

	return np.matrix(L)

def get_U(lu_mat):

	U = lu_mat.copy()

	for i in range(1, U.shape[0]):

		U[i, :i] = 0

	return np.matrix(U)


def solve_LU(lu_mat,b):

	y = np.matrix(np.zeros([lu_mat.shape[0], 1]))
	for i in range(y.shape[0]):
		y[i, 0] = b[i, 0] - lu_mat[i, :i] * y[:i]

	x = np.matrix(np.zeros([lu_mat.shape[0], 1]))
	for i in range(1, x.shape[0] + 1):
		x[-i, 0] = (y[-i] - lu_mat[-i, -i:] * x[-i:, 0] )/ lu_mat[-i, -i]

	return x


def invertible_matrix(LU):

	E = np.matrix(np.eye(len(LU)))
	new_LU = LU.copy()
	for i in range(len(E)):
		new_LU[:,i] = solve_LU(LU,E[:, i])

	return(new_LU)


def is_zero(mat,b):

	if np.allclose(b, np.zeros([b.shape[0],1])):
		b = mat[:, mat.shape[1]-1]
		mat = np.delete(mat, mat.shape[1]-1, 1)

	return(mat,b)


def det(U):

	determ = 1
	for i in range(U.shape[0]):
		determ *= U[i,i]

	return determ



if __name__ == '__main__':

    mat_original = np.matrix([[9, -5, -6, 3], [1, -7, 1, 0], [3, -4, 9, 0], [6, -1, 9, 8]]) #мой вариант
    #b = np.matrix([-8, 38, 47, -8]).transpose()
    #-8, 38, 47, -8
    b_original = np.matrix([-8, 38, 47, -8]).transpose() #мой вариант
    #b_original = np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]).transpose()
    mat, b = is_zero(mat_original, b_original)
    print("Исходная матрица A:\n", mat)
    print("Исходная матрица b:\n", b)
    LU = lu_decompose(mat)
    print("Матрица LU:\n", LU)
    L = get_L(LU)
    U = get_U(LU)
    print("Матрица L:\n", np.matrix(L))
    print("Матрица U:\n", np.matrix(U))
    print("Проверка L*U = A:")
    print("Верно" if np.allclose(L*U, mat) else "Не верно")
    x = solve_LU(LU,b)
    print("Результат решения СЛАУ:\n", x)
    print("Проверка правильности решения уравнения:")
    print(mat*x - b)
    print("Определитель матрицыx: ", det(U))
    inv_mat = invertible_matrix(LU)
    print("Обратная матрица:")
    print(inv_mat)
    print("Проверка правильности найденной обратной м-цы:")
    print(mat_original*inv_mat)
    
    