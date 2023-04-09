import numpy as np

def my_3(x): return(np.tan(x))#np.tan(x)#np.sin(3.1415926*x/6)


def fi(x):

	res = np.zeros(len(x))
	for i in range(0,len(x)):
		res[i] = my_3(x[i])

	return(res)


def w_d(x):

	res = np.zeros(len(x))
	w_k = 1
	for i in range(0,len(x)):
		for k in range(0,len(x)):
			if k!=i: w_k *= x[i]-x[k]
		res[i] = w_k
		w_k = 1

	return res


def f_w(f,w):

	res = np.zeros(len(f))
	for i in range(0,len(f)):
		res[i] = f[i]/w[i]

	return res


def X_x(X,x_i):

	res = np.zeros(len(x_i))
	for i in range(0,len(x_i)):
		res[i] = X - x_i[i]

	return res


def p_s(x): 
	
	res = '' if x<0 else ' '
	return(res)


def print_table(X_i, fi_, w, fw, Xminusx):

	print('\n')
	print("i", "x_i\t", "   f_i\t", "   w(x_i)", "   f_i/w(x_i)", " X-x_i")
	print("================================================")
	for i in range(0,len(X_i)):
		print(i, '%.2f ' % X_i[i], p_s(fi_[i]), '%.5f ' % fi_[i], p_s(w[i]), '%.3f  ' % w[i], p_s(fw[i]), '%.5f  ' % fw[i], p_s(Xminusx[i]), '%.1f' % Xminusx[i])


def Lag_pol_print(fw,x_i):

	#L = 1
	print()
	print("Многочлен Лагранжа:")
	print("L = ", end='')
	for i in range(0,len(x_i)):
		print('%.5f' % fw[i], end='') if (fw[i]<0)  or (i==0) else print('+%.5f' % fw[i],end='')
		for k in range(0,len(x_i)):
			if k!=i: print('*(x-{0})'.format(x_i[k]), end='')
	print('\n')


def Lag_pol(fw,X_i,X):

	L = 1
	sum = 0
	for i in range(0,len(fw)):
		L *= fw[i]
		for k in range(0,len(fw)):
			if k!=i: 
				L *= (X-X_i[k])
		sum += L
		L = 1

	return(sum)


def Newton_pol(X_i,fi_):

	len_x = len(X_i)
	res = np.zeros(len(fi_)-1)
	len_r = len(res)
	i = 0
	while i<len_r:
		res[i] = (fi_[i]-fi_[i+1])/(X_i[i]-X_i[i+len_x-len_r])
		i+=1
	return res


def Newton_pol_res(X_i,fi_):

	res = np.zeros(len(X_i)-1)
	for i in range(0, len(res)):
		ff = Newton_pol(X_i, fi_)
		res[i] = ff[0]
		fi_ = ff.copy()

	return res


def print_Newton(ff,X_i,fi_):

	#print(ff)
	print("P(x) =", fi_[0],end='')
	for i in range(0, len(ff)):
		print("{0}x".format(ff[i]) if (ff[i]<0) else  "+{0}x".format(ff[i]),end='')
		if i!=0:
			k=1
			while k<=i:
				print("(x-{0})".format(X_i[k]),end='')
				k+=1
	print('\n')

def New(ff,X_i,x,fi_):

	res = fi_[0]+x*ff[0]
	for i in range(1, len(ff)):
		proizv = 1
		k = 1
		while k<=i:
			proizv *= (x-X_i[k])
			k+=1
		res += ff[i]*x*proizv

	return(res) 
		

if __name__ == '__main__':

	#X_i = [0.1, 0.5, 0.9, 1.3]#методичка
	pi = 3.1415926
	#cleX_i = [0, pi/8, 2*pi/8, 3*pi/8]#мой вариант а
	X_i = [0, pi/8, pi/3, 3*pi/8]# мой вариант б
	#X = 0.8#методичка	
	X = 3*pi/16
	fi_ = fi(X_i)
	w = w_d(X_i)
	fw = f_w(fi_,w)
	Xminusx = X_x(X, X_i)
	print_table(X_i,fi_,w,fw,Xminusx)
	Lag_pol_print(fw,X_i)
	L = Lag_pol(fw,X_i,X)
	print("L({0}) = {1}".format(X,L))
	Y = my_3(X)
	print("y({0}) = {1}".format(X,Y))
	print("delta(L) =", abs(L-Y))

	#X_i = [0, 1.0, 2.0, 3.0]
	#X = 1.5
	#fi_ = fi(X_i)
	ff = Newton_pol_res(X_i,fi_)
	print()
	print("Многочлен Ньютона:")
	print_Newton(ff,X_i,fi_)
	res = New(ff,X_i,X,fi_)
	print("P =",res)
	print("y({0}) = {1}".format(X,Y))
	print("delta(L) =", abs(res-Y))
	print()









