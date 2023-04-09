import scipy as np
import matplotlib.pyplot as plt
from numpy.linalg import linalg as lin
from numpy import abs
class DSolve:

    def __init__(self):
        self.boundaryCondition = []
        self.initialCondition = []
        self.h = 0.1
        self.DX = (0,3.24)

    def setStep(self,h):
        self.h = h


    def Runge_Kutta(self,sys,iCond):
        def getAlpha(func1,func2,dots):
            _x,_y,_z = dots
            _K = [0,0,0,0]
            _L = [0,0,0,0]
            _K[0] = h*func1(_x, _y, _z)
            _L[0] = h*func2(_x, _y, _z)

            _K[1] = h*func1(_x + 0.5*h,_y + 0.5*_K[0],_z + 0.5*_L[0])
            _L[1] = h*func2(_x + 0.5*h,_y + 0.5*_K[0],_z + 0.5*_L[0])

            _K[2] = h*func1(_x + 0.5*h,_y + 0.5*_K[1],_z + 0.5*_L[1])
            _L[2] = h*func2(_x + 0.5*h,_y + 0.5*_K[1],_z + 0.5*_L[1])

            _K[3] = h*func1(_x + 0.5*h,_y + 0.5*_K[2],_z + 0.5*_L[2])
            _L[3] = h*func2(_x + 0.5*h,_y + 0.5*_K[2],_z + 0.5*_L[2])
            return _K,_L
        self.initialCondition = iCond
        y0 = iCond[0][1]
        z0 = iCond[1][1]
        f = sys[0]
        g = sys[1]
        y = []
        z = []
        y.append(y0)
        z.append(z0)
        h = self.h
        x = np.arange(self.DX[0], self.DX[1], h)
        count = len(x)-1
        for i in range(count):
            K,L = getAlpha(f,g,(x[i],y[i],z[i]))
            dz = 1/6*(K[0]+2*K[1]+2*K[2]+K[3])
            dy = 1/6*(L[0]+2*L[1]+2*L[2]+L[3])
            currentZ = z[i] + dz
            currentY = y[i] + dy
            y.append(currentY)
            z.append(currentZ)
        return x, y ,z

    def Shoting(self,sys,bCond):
        def nextE(j):
            return E[j]-((E[j]-E[j-1])/(r[j]-r[j-1]))*(r[j]-rr)

        self.boundaryCondition = bCond
        eps = 0.0001
        l = bCond[0][1]
        rr = bCond[1][1]
        fr = bCond[1][0]
        f = sys[0]
        g = sys[1]
        h = self.h
        E = []
        r = []
        E1 = 3.0
        E2 = 2.8
        E.extend([E1,E2])
        x1,y1,z1 = self.Runge_Kutta(sys,[(bCond[0][0],E1),bCond[0]])
        x2,y2,z2 = self.Runge_Kutta(sys,[(bCond[0][0],E2),bCond[0]])
        r.extend([fr(x1[-1],y1[-1],z1[-1]),fr(x2[-1],y2[-1],z2[-1])])
        curE = E2
        curR = fr(x2[-1],y2[-1],z2[-1])
        curx = x2
        cury = y2
        curz = z2
        i = 1
        while abs(r[-1] - rr) > 0.001:
            curE = nextE(i)
            curx,cury,curz = self.Runge_Kutta(sys,[(bCond[0][0],curE),bCond[0]])
            x = curx[-1]
            y = cury[-1]
            z = curz[-1]
            r.append(fr(x,y,z))
            E.append(curE)
            i += 1
        return curx,cury

    def Ultimate_Difference(self,sys,bCond):
        def p(x):
            return 4*x/(2*x+1)
        def q(x):
            return -4/(2*x+1)
        def f(x):
            return 0
        def M():
            res = np.zeros((count,count))
            for i in range(1,count-1):
                res[i][i - 1] = (1 - (p(x[i + 1]) * h / 2))
                res[i][i] = (-2 + (h ** 2) * q(x[i + 1]))
                res[i][i + 1] = (1 + (p(x[i + 1]) * h / 2))
            res[0][0] = -1/h;
            res[0][1] = 1/h
            res[-1][count-2] = -1/h
            res[-1][count-1]= 1/h+2
            return res
        def b():
            fx = f(0)
            res = np.zeros((count,1))
            for i in range(1,count-1):
                res[i][0] = (h ** 2) * fx
            res[0][0] = -1
            res[-1][0] = 3
            return res

        self.boundaryCondition = bCond
        h = self.h
        fl = bCond[0][1]
        fr = bCond[1][1]
        x = np.arange(bCond[0][0],1, h)
        count = len(x)
        res = lin.solve(M(),b())
        return x,res

    def Runge_Romberg(self,yh,y2h,p):
        ans = []
        i = 0
        while i < len(yh):
            delta = abs(yh[i]-(y2h[i//2] if i%2 else 0))/(2**p-1)
            ans.append(delta)
            i += 1
        return ans
