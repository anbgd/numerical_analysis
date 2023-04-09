import numpy as np

def get_PQ(mat):
    P1 = -mat[0,1]/mat[0,0]
    Q1 = answ[0,0]/mat[0,0]
    size_mat = mat[:,0].size
    i = 1
    P = [P1]
    Q = [Q1]
    while i < size_mat-1:
        P.append(-mat[i,i+1]/(mat[i,i]+mat[i,i-1]*P[i-1]))
        Q.append((answ[i,0]-mat[i,i-1]*Q[i-1])/(mat[i,i]+mat[i,i-1]*P[i-1]))
        i+=1
    P.append(0)
    Q.append((answ[i,0]-mat[i,i-1]*Q[i-1])/(mat[i,i]+mat[i,i-1]*P[i-1]))
    return P, Q

def method(size_mat, P, Q):
    X = np.zeros(size_mat)
    X[size_mat-1] = Q[size_mat-1]
    i = size_mat-2
    while i > -1:
        X[i] = X[i+1]*P[i]+Q[i]
        i-=1
    return X


if __name__ == '__main__':

	mat = [[13, -5, 0, 0, 0], [-4, 9, -5, 0, 0], [0, -1, -12, -6, 0], [0, 0, 6, 20, -5], [0, 0, 0, 4, 5]]
	answ = [[66], [47], [43], [74], [-14]]
	answ = np.asarray(answ)#.reshape(-1,1)
	mat = np.asarray(mat)
	print("mat = ")
	print(mat)
	print("answ = ")
	print(answ)
	P, Q = get_PQ(mat)
	size_mat = mat[:,0].size
	X = method(size_mat, P, Q)
	X = X.reshape(-1,1)
	print("X = ")
	print(X)
	print("Проверка правильности ответа: ")
	print(np.matrix(mat)*np.matrix(X) - answ)

	
	