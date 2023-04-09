import numpy as np

def norm_matrix(matrix):
    norm = 0
    matrix_tmp = matrix.reshape(1,-1)[0]
    for i in np.arange(0,matrix_tmp.size):
        norm += matrix_tmp[i]*matrix_tmp[i]
    return np.sqrt(norm)

def get_BC(alpha):
    alpha_size = alpha[:,0].size
    B = np.zeros(alpha_size*alpha_size).reshape(alpha_size,alpha_size)
    C = np.zeros(alpha_size*alpha_size).reshape(alpha_size,alpha_size)
    for i in np.arange(0, alpha_size):
        for j in np.arange(0, alpha_size):
            if i <= j:
                C[i,j] = alpha[i,j]
                B[i,j] = 0
            else:
                C[i,j] = 0
                B[i,j] = alpha[i,j]
    return B,C

def Zedel(alpha, betta, epsilon):
    k = 0
    if norm_matrix(alpha) < 1:
        x = betta.copy()
        alpha_size = alpha[:,0].size
        while True:
            pre_x = x.copy()
            for i in np.arange(0, alpha_size):
                tmp = 0
                for j in np.arange(0, alpha_size):
                    tmp += alpha[i,j]*x[j,0]
                x[i,0] = betta[i,0] + tmp
            epsilon_i = norm_matrix(pre_x - x)
            k += 1
            if epsilon_i < epsilon or k > 100:
                break

    else:
        print("Метод Зейделя не сходится")
    print("Число итераций: ", k)
    return x

def iteration(alpha, betta, epsilon):
    k = 0
    if norm_matrix(alpha) < 1:
        x = betta
        while True:
            pre_x = x
            x = betta + np.dot(alpha, pre_x)
            epsilon_i = norm_matrix(pre_x - x)
            k += 1
            if epsilon_i < epsilon or k > 100:
                break
    else:
        print("Метод итераций не сходится")
    print("Число итераций: ", k)
    return x

def get_alpha_betta(mat, answ):
    i = 0
    mat_size = mat[:,0].size
    alpha = np.zeros(mat_size*mat_size).reshape(mat_size,mat_size)
    betta = np.zeros(mat_size).reshape(-1,1)
    for row in mat:
        betta[i] = answ[i]/mat[i,i]
        alpha[i] = -row/mat[i,i]
        alpha[i,i] = 0
        i += 1
    return alpha, betta


if __name__ == '__main__':

    epsilon = 0.0001 #0.001 
    mat = [[-23, -7, 5, 2], [-7, -21, 4, 9], [9, 5, -31, -8], [0, 1, -2, 10]]
    answ = [[-26], [-55], [-58], [-24]]
    #print(answ)
    answ = np.asarray(answ)#.reshape(-1,1)
    mat = np.asarray(mat)
    print("mat = ")
    print(mat)
    print("answ = ")
    print(answ)
    print("eps = ", epsilon)
    alpha, betta = get_alpha_betta(mat,answ)
    print("Результат по методу итераций:")
    x1 = iteration(alpha, betta, epsilon)
    print(x1)
    print("Проверка:")
    print(mat@x1-answ)
    print("Результат по методу Зейделя:")
    x2 = Zedel(alpha, betta, epsilon)
    print(x2)
    print("Проверка:")
    print(mat@x2-answ)

	



