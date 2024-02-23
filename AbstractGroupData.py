allAGindex=[[1,1], [2,1], [3,1], [4,1], [4,2], [6,1], [6,2], [8,1], [8,2], [8,3],
            [8,4], [8,5], [12,1], [12,2], [12,3], [12,4], [12,5], [12,6], [16,1], [16,2],
            [16,3], [16,4], [16,5], [16,6], [16,7], [16,8], [16,9], [16,10], [16,11], [16,12],
            [16,13], [16,14], [24,1], [24,2], [24,3], [24,4], [24,5], [24,6], [24,7], [24,8],
            [24,9], [24,10], [24,11], [24,12], [32,1], [32,2], [32,3], [32,4], [32,5], [32,6],
            [32,7], [32,8], [32,9], [32,10], [32,11], [32,12], [32,13], [32,14], [32,15], [32,16],
            [32,17], [48,1], [48,2], [48,3], [48,4], [48,5], [48,6], [48,7], [48,8], [48,9],
            [48,10], [48,11], [48,12], [48,13], [48,14], [48,15], [64,1], [64,2], [64,3], [64,4],
            [64,5], [96,1], [96,2], [96,3], [96,4], [96,5], [96,6], [96,7], [96,8], [96,9],
            [96,10], [192,1], [192,2]]

from sympy import Matrix,symbols, exp, sqrt, I, pi, conjugate, kronecker_product

# 定义复数
omega = exp(I * 2 * pi / 3)
theta = (1 + I) / sqrt(2)
sigma = exp(I * pi / 8)

otimesG21 = Matrix([[1, 1], [1, -1]])
otimesG41 = Matrix([[1, 1, 1, 1], [1, I, -1, -I], [1, -1, 1, -1], [1, -I, -1, I]])
otimesG42 =Matrix([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]])

AGCharTab = {
    (1, 1): Matrix([[1]]),
    (2, 1): Matrix([[1, 1], [1, -1]]),
    (3, 1): Matrix([[1, 1, 1], [1, omega, conjugate(omega)], [1, conjugate(omega), omega]]),
    (4, 1): Matrix([[1, 1, 1, 1], [1, I, -1, -I], [1, -1, 1, -1], [1, -I, -1, I]]),
    (4, 2): Matrix([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]]),
    (6, 1): Matrix([[1, 1, 1, 1, 1, 1], [1, -omega.conjugate(), omega, -1, omega.conjugate(), -omega],
                    [1, omega, conjugate(omega), 1, omega, conjugate(omega)],
                    [1, -1, 1, -1, 1, -1], [1, omega.conjugate(), omega, 1, omega.conjugate(), omega],
                    [1, -omega, conjugate(omega), -1, omega, -omega]]),
    (6, 2): Matrix([[1, 1, 1], [1, 1, -1], [2, -1, 0]]),
    (8, 1): Matrix([[1, 1, 1, 1, 1, 1, 1, 1], [1, (1 + I) / sqrt(2), I, -(1 + I) / sqrt(2), -1, -(1 + I) / sqrt(2), -I, (1 + I) / sqrt(2)],
                    [1, I, -1, -I, 1, I, -1, -I], [1, -(1 + I) / sqrt(2), -I, (1 + I) / sqrt(2), -1, (1 + I) / sqrt(2), I, -(1 + I) / sqrt(2)],
                    [1, -1, 1, -1, 1, -1, 1, -1], [1, (1 - I) / sqrt(2), I, (1 - I) / sqrt(2), -1, (1 - I) / sqrt(2), -I, (1 - I) / sqrt(2)],
                    [1, -I, -1, I, 1, -I, -1, I], [1, -(1 - I) / sqrt(2), -I, -(1 - I) / sqrt(2), -1, -(1 - I) / sqrt(2), I, -(1 - I) / sqrt(2)]]),
    (8, 2): Matrix([[1, 1, 1, 1, 1, 1, 1, 1], [1, I, -1, -I, 1, I, -1, -I],
                    [1, -1, 1, -1, 1, -1, 1, -1], [1, -I, -1, I, 1, -I, -1, I],
                    [1, 1, 1, 1, -1, -1, -1, -1], [1, I, -1, -I, -1, -I, 1, I],
                    [1, -1, 1, -1, -1, 1, -1, 1], [1, -I, -1, I, -1, I, 1, -I]]),
    (8, 3): Matrix([[1, 1, 1, 1, 1, 1, 1, 1], [1, -1, 1, -1, 1, -1, 1, -1],
                    [1, 1, -1, -1, 1, 1, -1, -1], [1, -1, -1, 1, 1, -1, -1, 1],
                    [1, 1, 1, 1, -1, -1, -1, -1], [1, -1, 1, -1, -1, 1, -1, 1],
                    [1, 1, -1, -1, -1, -1, 1, 1], [1, -1, -1, 1, -1, 1, 1, -1]]),
    (8, 4): Matrix([[1, 1, 1, 1, 1], [1, 1, 1, -1, -1], [1, 1, -1, 1, -1], [1, 1, -1, -1, 1], [2, -2, 0, 0, 0]]),
    (8, 5): Matrix([[1, 1, 1, 1, 1], [1, 1, 1, -1, -1], [1, 1, -1, 1, -1], [1, 1, -1, -1, 1], [2, -2, 0, 0, 0]]),
    (12, 1): Matrix([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                     [1, -I * omega.conjugate(), -conjugate(omega), I, omega, -I * omega.conjugate(), -1, I * omega, conjugate(omega), -I, -omega, I * omega.conjugate()],
                     [1, -conjugate(omega), omega, -1, conjugate(omega), -omega, 1, -conjugate(omega), omega, -1, conjugate(omega), -omega],
                     [1, I, -1, -I, 1, I, -1, -I, 1, I, -1, -I],
                     [1, omega, conjugate(omega), 1, omega, conjugate(omega), 1, omega, conjugate(omega), 1, omega, conjugate(omega)],
                     [1, -I * omega.conjugate(), -conjugate(omega), I, omega, -I * omega, -1, I * omega.conjugate(), conjugate(omega), -I, -conjugate(omega), I * omega],
                     [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
                     [1, I * omega, -conjugate(omega), -I, -omega, I * omega.conjugate(), -1, -I * omega, conjugate(omega), I, -conjugate(omega), -I * omega],
                     [1, omega, conjugate(omega), 1, omega, conjugate(omega), -1, -omega, -conjugate(omega), -1, -omega, -conjugate(omega)],
                     [1, -I, -1, I, 1, -I, -1, I, 1, -I, -1, I],
                     [1, -omega, conjugate(omega), -1, conjugate(omega), -omega, 1, -omega, conjugate(omega), -1, omega, -conjugate(omega)],
                     [1, I * omega.conjugate(), -conjugate(omega), -I, omega, I * omega, -1, -I * omega.conjugate(), conjugate(omega), I, -conjugate(omega), -I * omega]]),
    (12, 2): Matrix([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                     [1, omega, conjugate(omega), 1, omega, conjugate(omega), 1, omega, conjugate(omega), 1, omega, conjugate(omega)],
                     [1, conjugate(omega), omega, 1, conjugate(omega), omega, 1, conjugate(omega), omega, 1, conjugate(omega), omega],
                     [1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1],
                     [1, omega, conjugate(omega), -1, -omega, -conjugate(omega), 1, omega, conjugate(omega), -1, -omega, -conjugate(omega)],
                     [1, conjugate(omega), omega, -1, -conjugate(omega), -omega, 1, conjugate(omega), omega, -1, -conjugate(omega), -omega],
                     [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1],
                     [1, omega, conjugate(omega), 1, omega, conjugate(omega), -1, -omega, -conjugate(omega), -1, -omega, -conjugate(omega)],
                     [1, conjugate(omega), omega, 1, conjugate(omega), omega, -1, -conjugate(omega), -omega, -1, -conjugate(omega), -omega],
                     [1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1],
                     [1, omega, conjugate(omega), -1, -omega, -conjugate(omega), -1, -omega, -conjugate(omega), 1, omega, conjugate(omega)],
                     [1, conjugate(omega), omega, -1, -conjugate(omega), -omega, -1, -conjugate(omega), -omega, 1, omega, conjugate(omega)]]),
    (12, 3): Matrix([[1, 1, 1, 1, 1, 1], [1, 1, 1, 1, -1, -1], [1, -1, -1, 1, 1, -1],
                     [1, -1, -1, 1, -1, 1], [2, 2, -1, -1, 0, 0], [2, -2, 1, -1, 0, 0]]),
    (12, 4): Matrix([[1, 1, 1, 1, 1, 1], [1, 1, 1, 1, -1, -1], [1, -1, -1, 1, I, -I],
                     [1, -1, -1, 1, -I, I], [2, 2, -1, -1, 0, 0], [2, -2, 1, -1, 0, 0]]),
    (12, 5): Matrix([[1, 1, 1, 1], [1, 1, omega, conjugate(omega)], [1, 1, conjugate(omega), omega], [3, -1, 0, 0]]),
    (12, 6): Matrix([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                     [1, -conjugate(omega), omega, -1, conjugate(omega), -omega, 1, -conjugate(omega), omega, -1, conjugate(omega), -omega],
                     [1, omega, conjugate(omega), 1, omega, conjugate(omega), 1, omega, conjugate(omega), 1, omega, conjugate(omega)],
                     [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
                     [1, conjugate(omega), omega, 1, conjugate(omega), omega, 1, conjugate(omega), omega, 1, conjugate(omega), omega],
                     [1, -omega, conjugate(omega), -1, conjugate(omega), -omega, 1, -omega, conjugate(omega), -1, omega, -conjugate(omega)],
                     [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1],
                     [1, -conjugate(omega), omega, -1, conjugate(omega), -omega, -1, conjugate(omega), -omega, 1, -conjugate(omega), omega],
                     [1, omega, conjugate(omega), 1, omega, conjugate(omega), -1, -omega, -conjugate(omega), -1, -omega, -conjugate(omega)],
                     [1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1],
                     [1, conjugate(omega), omega, 1, conjugate(omega), omega, -1, -conjugate(omega), -omega, -1, -conjugate(omega), -omega],
                     [1, -omega, conjugate(omega), -1, conjugate(omega), -omega, -1, conjugate(omega), omega, 1, -omega, conjugate(omega)]]),
    (16, 1): Matrix([[sigma**((s-1)*(t-1)) for s in range(1, 17)] for t in range(1, 17)]),
}
AGCharTab[(16, 2)]=kronecker_product(otimesG21,AGCharTab[(8,1)])
AGCharTab[(16, 3)]=kronecker_product(otimesG41,AGCharTab[(4,1)])
AGCharTab[(16, 4)]=kronecker_product(otimesG21,AGCharTab[(8,2)])
AGCharTab[(16, 5)]=kronecker_product(otimesG21,AGCharTab[(8,3)])

AGCharTab[(16,6)]=Matrix([[1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,-1,-1,-1,-1],
                 [1,1,1,1,-1,-1,1,1,-1,-1],[1,1,1,1,-1,-1,-1,-1,1,1],
                 [1,-1,1,-1,I,-I,1,-1,I,-I],[1,-1,1,-1,I,-I,-1,1,-I,I],
                 [1,-1,1,-1,-I,I,1,-1,-I,I],[1,-1,1,-1,-I,I,-1,1,I,-I],
                 [2,2*I,-2,-2*I,0,0,0,0,0,0],[2,-2*I,-2,2*I,0,0,0,0,0,0]])
AGCharTab[(16,7)]=Matrix([[1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,-1,-1,-1,-1],
                 [1,1,1,1,-1,-1,1,1,-1,-1],[1,1,1,1,-1,-1,-1,-1,1,1],
                 [1,-1,1,-1,1,-1,1,-1,1,-1],[1,-1,1,-1,1,-1,-1,1,-1,1],
                 [1,-1,1,-1,-1,1,1,-1,-1,1],[1,-1,1,-1,-1,1,-1,1,1,-1],
                 [2,2*I,-2,-2*I,0,0,0,0,0,0],[2,-2*I,-2,2*I,0,0,0,0,0,0]])



AGCharTab[(16,8)]=Matrix([[1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,-1,-1,-1,-1],
                  [1,1,1,1,-1,-1,1,-1,1,-1],[1,1,1,1,-1,-1,-1,1,-1,1],
                  [1,1,-1,-1,-1,1,I,-I,-I,I],[1,1,-1,-1,-1,1,-I,I,I,-I],
                  [1,1,-1,-1,1,-1,-I,-I,I,I],[1,1,-1,-1,1,-1,I,I,-I,-I],
                  [2,-2,2,-2,0,0,0,0,0,0],[2,-2,-2,2,0,0,0,0,0,0]])
AGCharTab[(16,9)]=kronecker_product(otimesG21,AGCharTab[(8,4)])
AGCharTab[(16,10)]=Matrix([[1,1,1,1,1,1,1,1,1,1],[1,1,1,1,-1,-1,-1,-1,1,1],
                  [1,1,1,1,1,1,-1,-1,-1,-1],[1,1,1,1,-1,-1,1,1,-1,-1],
                  [1,1,-1,-1,1,-1,I,-I,I,-I],[1,1,-1,-1,1,-1,-I,I,-I,I],
                  [1,1,-1,-1,-1,1,-I,I,I,-I],[1,1,-1,-1,-1,1,I,-I,-I,I],
                  [2,-2,2,-2,0,0,0,0,0,0],[2,-2,-2,2,0,0,0,0,0,0]])
AGCharTab[(16,11)]=kronecker_product(otimesG21,AGCharTab[(8,5)])
AGCharTab[(16,12)]=Matrix([[1,1,1,1,1,1,1],[1,1,1,1,1,-1,-1],[1,1,1,-1,-1,1,-1],[1,1,1,-1,-1,-1,1],
                  [2,2,-2,0,0,0,0],[2,-2,0,sqrt(2),-sqrt(2),0,0],[2,-2,0,-sqrt(2),sqrt(2),0,0]])
AGCharTab[(16,13)]=Matrix([[1,1,1,1,1,1,1],[1,1,1,1,1,-1,-1],[1,1,1,-1,-1,1,-1],[1,1,1,-1,-1,-1,1],
                  [2,2,-2,0,0,0,0],[2,-2,0,I*sqrt(2),-I*sqrt(2),0,0],[2,-2,0,-I*sqrt(2),I*sqrt(2),0,0]])

AGCharTab[(16,14)]=Matrix([[1,1,1,1,1,1,1],[1,1,1,1,1,-1,-1],[1,1,1,-1,-1,1,-1],[1,1,1,-1,-1,-1,1],
                  [2,2,-2,0,0,0,0],[2,-2,0,sqrt(2),-sqrt(2),0,0],[2,-2,0,-sqrt(2),sqrt(2),0,0]])
AGCharTab[(24,1)]=Matrix([[1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,-1,-1],[1,1,1,1,-1,-1,-1,1,-1],
                 [1,1,1,1,-1,-1,-1,-1,1],[2,2,-1,-1,-2,1,1,0,0],[2,2,-1,-1,2,-1,-1,0,0],
                 [2,-2,-1,1,0,I*sqrt(3),-I*sqrt(3),0,0],[2,-2,-1,1,0,-I*sqrt(3),I*sqrt(3),0,0],[2,-2,2,-2,0,0,0,0,0]])
AGCharTab[(24,2)]=Matrix([[1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,-1,-1],[1,1,1,1,-1,-1,-1,1,-1],
                 [1,1,1,1,-1,-1,-1,-1,1],[2,2,-1,-1,-2,1,1,0,0],[2,2,-1,-1,2,-1,-1,0,0],
                 [2,-2,-1,1,0,sqrt(3),-sqrt(3),0,0],[2,-2,-1,1,0,-sqrt(3),sqrt(3),0,0],[2,-2,2,-2,0,0,0,0,0]])
AGCharTab[(24,3)]=kronecker_product(otimesG21,AGCharTab[(12,4)])
AGCharTab[(24,4)]=kronecker_product(otimesG41,AGCharTab[(6,2)])
AGCharTab[(24,5)]=kronecker_product(otimesG42,AGCharTab[(6,2)])
AGCharTab[(24,6)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                 [1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1],
                 [1,-omega.conjugate(),omega,-1,omega.conjugate(),-omega,1,omega,omega.conjugate(),1,-omega.conjugate(),omega,-1,omega.conjugate(),-omega],
                 [1,-omega.conjugate(),omega,-1,omega.conjugate(),-omega,1,omega,omega.conjugate(),-1,omega.conjugate(),-omega,1,-omega.conjugate(),omega],
                 [1,omega,omega.conjugate(),1,omega,omega.conjugate(),1,omega.conjugate(),omega,1,omega,omega.conjugate(),1,omega,omega.conjugate()],
                 [1,omega,omega.conjugate(),1,omega,omega.conjugate(),1,omega.conjugate(),omega,-1,-omega,-omega.conjugate(),-1,-omega,-omega.conjugate()],
                 [1,-1,1,-1,1,-1,1,1,1,1,-1,1,-1,1,-1],
                 [1,-1,1,-1,1,-1,1,1,1,-1,1,-1,1,-1,1],
                 [1,omega.conjugate(),omega,1,omega.conjugate(),omega,1,omega,omega.conjugate(),1,omega.conjugate(),omega,1,omega.conjugate(),omega],
                 [1,omega.conjugate(),omega,1,omega.conjugate(),omega,1,omega,omega.conjugate(),-1,-omega.conjugate(),-omega,-1,-omega.conjugate(),-omega],
                 [1,-omega,omega.conjugate(),-1,omega,-omega.conjugate(),1,omega.conjugate(),omega,1,-omega,omega.conjugate(),-1,omega,-omega.conjugate()],
                 [1,-omega,omega.conjugate(),-1,omega,-omega.conjugate(),1,omega.conjugate(),omega,-1,omega,-omega.conjugate(),1,-omega,omega.conjugate()],
                 [2,0,-2*omega.conjugate(),0,2*omega,0,-2,2*omega.conjugate(),-2*omega,0,0,0,0,0,0],
                 [2,0,-2,0,2,0,-2,2,-2,0,0,0,0,0,0],
                 [2,0,-2*omega,0,2*omega.conjugate(),0,-2,2*omega,-2*omega.conjugate(),0,0,0,0,0,0]])
AGCharTab[(24,7)]=Matrix([[1,1,1,1,1],[1,1,-1,-1,1],[2,2,0,0,-1],[3,-1,1,-1,0],[3,-1,-1,1,0]])
AGCharTab[(24,8)]=Matrix([[1,1,1,1,1,1,1,1],[1,1,1,1,omega,omega,omega.conjugate(),omega.conjugate()],
                 [1,1,1,1,omega.conjugate(),omega.conjugate(),omega,omega],[1,-1,-1,1,1,-1,1,-1],
                 [1,-1,-1,1,omega,-omega,omega.conjugate(),-omega.conjugate()],[1,-1,-1,1,omega.conjugate(),-omega.conjugate(),omega,-omega],
                 [3,3,-1,-1,0,0,0,0],[3,-3,1,-1,0,0,0,0]])
AGCharTab[(24,9)]=Matrix([[1,1,1,1,1,1,1],[1,1,1,omega,omega,omega.conjugate(),omega.conjugate()],[1,1,1,omega.conjugate(),omega.conjugate(),omega,omega],[2,-2,0,1,-1,-1,1],
                 [2,-2,0,omega,-omega,-omega.conjugate(),omega.conjugate()],[2,-2,0,omega.conjugate(),-omega.conjugate(),-omega,omega],[3,3,-1,0,0,0,0]])
AGCharTab[(24,10)]=Matrix([[1,1,1,1,1,1,1,1],[1,1,omega,omega.conjugate(),1,1,omega,omega.conjugate()],
                  [1,1,omega.conjugate(),omega,1,1,omega.conjugate(),omega],[3,-1,0,0,3,-1,0,0],
                  [1,1,1,1,-1,-1,-1,-1],[1,1,omega,omega.conjugate(),-1,-1,-omega,-omega.conjugate()],
                  [1,1,omega.conjugate(),omega,-1,-1,-omega.conjugate(),-omega],[3,-1,0,0,-3,1,0,0]])
AGCharTab[(24,11)]=Matrix([[1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,-1,-1],
                  [1,1,-1,-1,1,1,-1,1,-1],[1,1,-1,-1,1,1,-1,-1,1],
                  [2,2,-1,-1,-1,-1,2,0,0],[2,2,1,1,-1,-1,-2,0,0],
                  [2,-2,0,0,-2,2,0,0,0],[2,-2,sqrt(3),-sqrt(3),1,-1,0,0,0],
                  [2,-2,-sqrt(3),sqrt(3),1,-1,0,0,0]])
AGCharTab[(24,12)]=kronecker_product(otimesG21,AGCharTab[(12,1)])
AGCharTab[(32,1)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1],
                 [1,1,1,1,-1,-1,-1,-1, 1, 1, 1,-1,1,-1],[1,1,1,1,-1,-1,-1,-1,1,1,-1,1,-1,1],
                 [1,1,1,1,1,1,-1,-1,-1,-1,1,1,-1,-1],[1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1],
                 [1,1,1,1,-1,-1,1,1,-1,-1,1,-1,-1,1],[1,1,1,1,-1,-1,1,1,-1,-1,-1,1,1,-1],
                 [2,2,-2,-2,0,0,0,0,2,-2,0,0,0,0],[2,-2,2,-2,2,-2,0,0,0,0,0,0,0,0],
                 [2,-2,2,-2,-2,2,0,0,0,0,0,0,0,0],[2,-2,-2,2,0,0,2*I,-2*I,0,0,0,0,0,0],
                 [2,-2,-2,2,0,0,-2*I,2*I,0,0,0,0,0,0],[2,2,-2,-2,0,0,0,0,-2,2,0,0,0,0]])
AGCharTab[(32,2)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1],
                 [1,1,1,1,-1,-1,1,1, -1,- 1, 1,-1,1,-1],[1,1,1,1,-1,-1,1,1,-1,-1,-1,1,-1,1],
                 [1,1,1,1,1,1,-1,-1,-1,-1,1,1,-1,-1],[1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1],
                 [1,1,1,1,-1,-1,-1,-1,1,1,1,-1,-1,1],[1,1,1,1,-1,-1,-1,-1,1,1,-1,1,1,-1],
                 [2,-2,2,-2,2,-2,0,0,0,0,0,0,0,0],[2,2,-2,-2,0,0,2,-2,0,0,0,0,0,0],
                 [2,-2,-2,2,0,0,0,0,2,-2,0,0,0,0],[2,-2,2,-2,-2,2,0,0,0,0,0,0,0,0],
                 [2,2,-2,-2,0,0,-2,2,0,0,0,0,0,0],[2,-2,-2,2,0,0,0,0,-2,2,0,0,0,0]])

AGCharTab[(32,3)]=kronecker_product(otimesG42,AGCharTab[(8,4)])
AGCharTab[(32,4)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                 [1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1],
                 [1,-1,1,-1,-1,1,I,-I,I,-I,-I,I,-1,1],
                 [1,-1,1,-1,-1,1,I,-I,I,-I,I,-I,1,-1],
                 [1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1],
                 [1,1,1,1,1,1,-1,-1,-1,-1,1,1,-1,-1],
                 [1,-1,1,-1,-1,1,-I,I,-I,I,I,-I,-1,1],
                 [1,-1,1,-1,-1,1,-I,I,-I,I,-I,I,1,-1],
                 [2,2,2,2,-2,-2,0,0,0,0,0,0,0,0],
                 [2,-2,2,-2,2,-2,0,0,0,0,0,0,0,0],
                 [2,2*I,-2,-2*I,0,0,1-I,1+I,-1+I,-1-I,0,0,0,0],
                 [2,2*I,-2,-2*I,0,0,-1+I,-1-I,1-I,1+I,0,0,0,0],
                 [2,-2*I,-2,2*I,0,0,1+I,1-I,-1-I,-1+I,0,0,0,0],
                 [2,-2*I,-2,2*I,0,0,-1-I,-1+I,1+I,1-I,0,0,0,0]])
AGCharTab[(32,5)]=kronecker_product(otimesG21,AGCharTab[(16,10)])
AGCharTab[(32,6)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                 [1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1],
                 [1,1,-1,-1,1,1,1,-1,-1,1,1,-1,1,-1],
                 [1,1,-1,-1,1,1,1,-1,-1,1,-1,1,-1,1],
                 [1,1,1,1,1,1,-1,-1,-1,-1,1,1,-1,-1],
                 [1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1],
                 [1,1,-1,-1,1,1,-1,1,1,-1,1,-1,-1,1],
                 [1,1,-1,-1,1,1,-1,1,1,-1,-1,1,1,-1],
                 [2,-2,0,0,2,-2,0,2*I,-2*I,0,0,0,0,0],
                 [2,-2,0,0,2,-2,0,-2*I,2*I,0,0,0,0,0],
                 [2,2,0,0,-2,-2,2,0,0,-2,0,0,0,0],
                 [2,2,0,0,-2,-2,-2,0,0,2,0,0,0,0],
                 [2,-2,2*I,-2*I,-2,2,0,0,0,0,0,0,0,0],
                 [2,-2,-2*I,2*I,-2,2,0,0,0,0,0,0,0,0]])
AGCharTab[(32,7)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                 [1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1],
                 [1,1,1,1,1,-1,-1,-1,1,1,1,-1,-1,-1],
                 [1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,1],
                 [1,1,1,1,-1,1,1,-1,1,1,-1,1,1,-1],
                 [1,1,1,1,-1,1,1,-1,-1,-1,1,-1,-1,1],
                 [1,1,1,1,-1,-1,-1,1,1,1,-1,-1,-1,1],
                 [1,1,1,1,-1,-1,-1,1,-1,-1,1,1,1,-1],
                 [2,-2,2,-2,0,0,0,0,0,0,0,2,-2,0],
                 [2,-2,2,-2,0,0,0,0,0,0,0,-2,2,0],
                 [2,-2,-2,2,0,0,0,0,2,-2,0,0,0,0],
                 [2,-2,-2,2,0,0,0,0,-2,2,0,0,0,0],
                 [2,2,-2,-2,0,2*I,-2*I,0,0,0,0,0,0,0],
                 [2,2,-2,-2,0,-2*I,2*I,0,0,0,0,0,0,0]])
AGCharTab[(32,8)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                 [1,1,1,1,-1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,-1,-1],
                 [1,1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,1,-1,-1,-1,-1],
                 [1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,1,-1,-1,-1,-1,1,1],
                 [1,-1,1,-1,1,-1,I,-I,I,-I,1,-1,1,-1,1,-1,I,-I,I,-I],
                 [1,-1,1,-1,-1,1,I,-I,-I,I,1,-1,1,-1,-1,1,I,-I,-I,I],
                 [1,-1,1,-1,1,-1,-I,I,-I,I,1,-1,1,-1,1,-1,-I,I,-I,I],
                 [1,-1,1,-1,-1,1,-I,I,I,-I,1,-1,1,-1,-1,1,-I,I,I,-I],
                 [2,-2*I,-2,2*I,0,0,0,0,0,0,2,-2*I,-2,2*I,0,0,0,0,0,0],
                 [2,2*I,-2,-2*I,0,0,0,0,0,0,2,2*I,-2,-2*I,0,0,0,0,0,0],
                 [1,I,-1,-I,I,-1,theta,-conjugate(theta),-conjugate(theta),-theta,-1,-I,1,I,-I,1,-theta,conjugate(theta),conjugate(theta),theta],
                 [1,I,-1,-I,I,-1,-theta,conjugate(theta),conjugate(theta),theta,-1,-I,1,I,-I,1,theta,-conjugate(theta),-conjugate(theta),-theta],
                 [1,-I,-1,I,I,1,conjugate(theta),-theta,theta,conjugate(theta),-1,I,1,-I,-I,-1,-conjugate(theta),theta,-theta,-conjugate(theta)],
                 [1,-I,-1,I,I,1,-conjugate(theta),theta,-theta,-conjugate(theta),-1,I,1,-I,-I,-1,conjugate(theta),-theta,theta,conjugate(theta)],
                 [1,I,-1,-I,-I,1,theta,-conjugate(theta),conjugate(theta),theta,-1,-I,1,I,I,-1,-theta,conjugate(theta),-conjugate(theta),-theta],
                 [1,I,-1,-I,-I,1,-theta,conjugate(theta),-conjugate(theta),-theta,-1,-I,1,I,I,-1,theta,-conjugate(theta),conjugate(theta),theta],
                 [1,-I,-1,I,-I,-1,conjugate(theta),-theta,-theta,-conjugate(theta),-1,I,1,-I,I,1,-conjugate(theta),theta,theta,conjugate(theta)],
                 [1,-I,-1,I,-I,-1,-conjugate(theta),theta,theta,conjugate(theta),-1,I,1,-I,I,1,conjugate(theta),-theta,-theta,-conjugate(theta)],
                 [2,-2,2,-2,0,0,0,0,0,0,-2,2,-2,2,0,0,0,0,0,0],
                 [2,2,2,2,0,0,0,0,0,0,-2,-2,-2,-2,0,0,0,0,0,0]])
AGCharTab[(32,9)]=kronecker_product(otimesG21,AGCharTab[(16,14)])
AGCharTab[(32,10)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,-1,-1,1,1,1,1,1,-1,-1],
                  [1,1,1,-1,-1,1,-1,1,1,1,-1,-1,1,-1],[1,1,1,-1,-1,-1,1,1,1,1,-1,-1,-1,1],
                  [2,2,-2,0,0,0,0,2,2,-2,0,0,0,0],[2,-2,0,sqrt(2),-sqrt(2),0,0,2,-2,0,sqrt(2),-sqrt(2),0,0],
                  [2,-2,0,-sqrt(2),sqrt(2),0,0,2,-2,0,-sqrt(2),sqrt(2),0,0],[1,1,1,1,1,I,I,-1,-1,-1,-1,-1,-I,-I],
                  [1,1,1,-1,-1,I,-I,-1,-1,-1,1,1,-I,I],[1,1,1,1,1,-I,-I,-1,-1,-1,-1,-1,I,I],
                  [1,1,1,-1,-1,-I,I,-1,-1,-1,1,1,I,-I],[2,-2,0,sqrt(2),-sqrt(2),0,0,-2,2,0,-sqrt(2),sqrt(2),0,0],
                  [2,-2,0,-sqrt(2),sqrt(2),0,0,-2,2,0,sqrt(2),-sqrt(2),0,0],[2,2,-2,0,0,0,0,-2,-2,2,0,0,0,0]])

AGCharTab[(32,11)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1],
                  [1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1],[1,1,1,1,1,1,-1,-1,-1,-1,1,1,-1,-1],
                  [2,2,-2,2,2,-2,0,0,0,0,0,0,0,0],[2,-2,0,2,-2,0,I*sqrt(2),I*sqrt(2),-I*sqrt(2),-I*sqrt(2),0,0,0,0],
                  [2,-2,0,2,-2,0,-I*sqrt(2),-I*sqrt(2),I*sqrt(2),I*sqrt(2),0,0,0,0],[1,1,-1,-1,-1,1,I,-I,I,-I,-1,1,I,-I],
                  [1,1,-1,-1,-1,1,I,-I,I,-I,1,-1,-I,I],[1,1,-1,-1,-1,1,-I,I,-I,I,1,-1,I,-I],
                  [1,1,-1,-1,-1,1,-I,I,-I,I,-1,1,-I,I],[2,2,2,-2,-2,-2,0,0,0,0,0,0,0,0],
                  [2,-2,0,-2,2,0,-sqrt(2),sqrt(2),sqrt(2),-sqrt(2),0,0,0,0],[2,-2,0,-2,2,0,sqrt(2),-sqrt(2),-sqrt(2),sqrt(2),0,0,0,0]])
AGCharTab[(32,12)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                  [1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,-1,-1,1,1],
                  [1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1],
                  [1,1,1,1,1,1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,-1,-1],
                  [1,1,-1,-1,1,1,-1,-1,1,-1,1,-1,I,-I,I,-I,I,-I,I,-I],
                  [1,1,-1,-1,1,1,-1,-1,1,-1,1,-1,-I,I,-I,I,-I,I,-I,I],
                  [1,1,-1,-1,1,1,-1,-1,-1,1,-1,1,-I,I,I,-I,-I,I,I,-I],
                  [1,1,-1,-1,1,1,-1,-1,-1,1,-1,1,I,-I,-I,I,I,-I,-I,I],
                  [2,-2,2,-2,2,-2,2,-2,0,0,0,0,0,0,0,0,0,0,0,0],
                  [2,-2,-2,2,2,-2,-2,2,0,0,0,0,0,0,0,0,0,0,0,0],
                  [1,1,I,I,-1,-1,-I,-I,1,I,-1,-I,theta,-conjugate(theta),theta,-conjugate(theta),-theta,conjugate(theta),-theta,conjugate(theta)],
                  [1,1,I,I,-1,-1,-I,-I,-1,-I,1,I,-theta,conjugate(theta),theta,-conjugate(theta),theta,-conjugate(theta),-theta,conjugate(theta)],
                  [1,1,I,I,-1,-1,-I,-I,1,I,-1,-I,-theta,conjugate(theta),-theta,conjugate(theta),theta,-conjugate(theta),theta,-conjugate(theta)],
                  [1,1,I,I,-1,-1,-I,-I,-1,-I,1,I,theta,-conjugate(theta),-theta,conjugate(theta),-theta,conjugate(theta),theta,-conjugate(theta)],
                  [1,1,-I,-I,-1,-1,I,I,1,-I,-1,I,-conjugate(theta),theta,-conjugate(theta),theta,conjugate(theta),-theta,conjugate(theta),-theta],
                  [1,1,-I,-I,-1,-1,I,I,1,-I,-1,I,conjugate(theta),-theta,conjugate(theta),-theta,-conjugate(theta),theta,-conjugate(theta),theta],
                  [1,1,-I,-I,-1,-1,I,I,-1,I,1,-I,conjugate(theta),-theta,-conjugate(theta),theta,-conjugate(theta),theta,conjugate(theta),-theta],
                  [1,1,-I,-I,-1,-1,I,I,-1,I,1,-I,-conjugate(theta),theta,conjugate(theta),-theta,conjugate(theta),-theta,-conjugate(theta),theta],
                  [2,-2,2*I,-2*I,-2,2,-2*I,2*I,0,0,0,0,0,0,0,0,0,0,0,0],
                  [2,-2,-2*I,2*I,-2,2,2*I,-2*I,0,0,0,0,0,0,0,0,0,0,0,0]])
AGCharTab[(32,13)]=Matrix([[1,1,1,1,1,1,1,1,1,1,    1,1,1,1,1,1,1,1,1,1],
                  [1,1,1,1,1,1,1,1,1,1,    1,1,-1,-1,-1,-1,-1,-1,-1,-1],
                  [1,-1,1,-1,1,-1,1,-1,1,1,     1,1,1,-1,1,-1,1,-1,1,-1],
                  [1,-1,1,-1,1,-1,1,-1,1,1,     1,1,-1,1,-1,1,-1,1,-1,1],
                  [1,I,-1,-I,I,-1,-I,1,-1,1,    -I,I,1,I,-1,-I,I,-1,-I,1],
                  [1,I,-1,-I,I,-1,-I,1,-1,1,    -I,I,-1,-I,1,I,-I,1,I,-1],
                  [1,-I,-1,I,I,1,-I,-1,-1,1,    -I,I,1,-I,-1,I,I,1,-I,-1],
                  [1,-I,-1,I,I,1,-I,-1,-1,1,    -I,I,-1,I,1,-I,-I,-1,I,1],
                  [1,1,1,1,-1,-1,-1,-1,1,1,     -1,-1,1,1,1,1,-1,-1,-1,-1],
                  [1,1,1,1,-1,-1,-1,-1,1,1,     -1,-1,-1,-1,-1,-1,1,1,1,1],
                  [1,-1,1,-1,-1,1,-1,1,1,1,     -1,-1,1,-1,1,-1,-1,1,-1,1],
                  [1,-1,1,-1,-1,1,-1,1,1,1,     -1,-1,-1,1,-1,1,1,-1,1,-1],
                  [1,I,-1,-I,-I,1,I,-1,-1,1,     I,-I,1,I,-1,-I,-I,1,I,-1],
                  [1,I,-1,-I,-I,1,I,-1,-1,1,     I,-I,-1,-I,1,I,I,-1,-I,1],
                  [1,-I,-1,I,-I,-1,I,1,-1,1,     I,-I,1,-I,-1,I,-I,-1,I,1],
                  [1,-I,-1,I,-I,-1,I,1,-1,1,     I,-I,-1,I,1,-I,I,1,-I,-1],
                  [2,0,-2,0,2,0,-2,0,2,-2,       2,-2,0,0,0,0,0,0,0,0],
                  [2,0,2,0,2*I,0,2*I,0,-2,-2,   -2*I,-2*I,0,0,0,0,0,0,0,0],
                  [2,0,-2,0,-2,0,2,0,2,-2,      -2,2,0,0,0,0,0,0,0,00],
                  [2,0,2,0,-2*I,0,-2*I,0,-2,-2,   2*I,2*I,0,0,0,0,0,0,0,0]])

AGCharTab[(32,14)]=kronecker_product(otimesG21,AGCharTab[(16,8)])
          
AGCharTab[(32,15)]=kronecker_product(otimesG21,AGCharTab[(16,11)])

AGCharTab[(32,16)]=kronecker_product(otimesG41,AGCharTab[(8,5)])
          
AGCharTab[(32,17)]=kronecker_product(otimesG41,AGCharTab[(8,1)])
AGCharTab[(48,1)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,-1,1,1,-1,1,-1,-1,-1,1,-1,1,1,-1],
                 [1,1,-1,1,1,-1,-1,1,1,-1,1,-1,1,-1,1],[1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,-1],
                 [1,1,-1,1,1,-1,-1,1,1,1,-1,1,-1,1,-1],[1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,1,1],
                 [1,1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,-1,1],[1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1],
                 [2,2,2,-1,-1,-1,0,0,0,2,2,-1,-1,0,0],[2,2,-2,-1,-1,1,0,0,0,-2,2,1,-1,0,0],
                 [2,2,-2,-1,-1,1,0,0,0,2,-2,-1,1,0,0],[2,2,2,-1,-1,-1,0,0,0,-2,-2,1,1,0,0],
                 [2,-2,0,2,-2,0,0,2,-2,0,0,0,0,0,0],[2,-2,0,2,-2,0,0,-2,2,0,0,0,0,0,0],
                 [4,-4,0,-2,2,0,0,0,0,0,0,0,0,0,0]])
AGCharTab[(48,2)]=kronecker_product(otimesG21,AGCharTab[(24,1)])
          
AGCharTab[(48,3)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                 [1,1,1,1,1,1,omega,omega,omega,omega,omega.conjugate(),omega.conjugate(),omega.conjugate(),omega.conjugate()],
                 [1,1,1,1,1,1,omega.conjugate(),omega.conjugate(),omega.conjugate(),omega.conjugate(),omega,omega,omega,omega],
                 [1,-1,1,-1,-1,1,1,-1,1,-1,1,-1,1,-1],
                 [1,-1,1,-1,-1,1,omega,-omega,omega,-omega,omega.conjugate(),-omega.conjugate(),omega.conjugate(),-omega.conjugate()],
                 [1,-1,1,-1,-1,1,omega.conjugate(),-omega.conjugate(),omega.conjugate(),-omega.conjugate(),omega,-omega,omega,-omega],
                 [2,2*I,-2,-2*I,0,0,-1,-I,1,I,-1,-I,1,I],
                 [2,2*I,-2,-2*I,0,0,-omega,-I*omega,omega,I*omega,-omega.conjugate(),-I*omega.conjugate(),omega.conjugate(),I*omega.conjugate()],
                 [2,2*I,-2,-2*I,0,0,-omega.conjugate(),-I*omega.conjugate(),omega.conjugate(),I*omega.conjugate(),-omega,-I*omega,omega,I*omega],
                 [2,-2*I,-2,2*I,0,0,-1,I,1,-I,-1,I,1,-I],
                 [2,-2*I,-2,2*I,0,0,-omega.conjugate(),I*omega.conjugate(),omega.conjugate(),-I*omega.conjugate(),-omega,I*omega,omega,-I*omega],
                 [2,-2*I,-2,2*I,0,0,-omega,I*omega,omega,-I*omega,-omega.conjugate(),I*omega.conjugate(),omega.conjugate(),-I*omega.conjugate()],
                 [3,3,3,3,-1,-1,0,0,0,0,0,0,0,0],
                 [3,-3,3,-3,1,-1,0,0,0,0,0,0,0,0]])
AGCharTab[(48,4)]=kronecker_product(otimesG21,AGCharTab[(24,9)])

AGCharTab[(48,5)]=kronecker_product(otimesG21,AGCharTab[(24,8)])

AGCharTab[(48,6)]=Matrix([[1,1,1,1,1,1,1,1],[1,1,1,1,1,-1,-1,-1],
                 [2,2,2,-1,-1,0,0,0],[2,-2,0,-1,1,-I*sqrt(2),I*sqrt(2),0],
                 [2,-2,0,-1,1,I*sqrt(2),-I*sqrt(2),0],[3,3,-1,0,0,1,1,-1],
                 [3,3,-1,0,0,-1,-1,1],[4,-4,0,1,-1,0,0,0]])
AGCharTab[(48,7)]=Matrix([[1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,-1,-1,-1,-1],
                 [1,-1,-1,1,1,-1,-1,1,-1,1],[1,-1,-1,1,1,-1,1,-1,1,-1],
                 [2,2,2,2,-1,-1,0,0,0,0],[2,-2,-2,2,-1,1,0,0,0,0],
                 [3,3,-1,-1,0,0,1,1,-1,-1],[3,3,-1,-1,0,0,-1,-1,1,1],
                 [3,-3,1,-1,0,0,1,-1,-1,1],[3,-3,1,-1,0,0,-1,1,1,-1]])
AGCharTab[(48,8)]=Matrix([[1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,-1,-1,-1,-1],
                 [1,-1,1,-1,1,-1,I,-I,I,-I],[1,-1,1,-1,1,-1,-I,I,-I,I],
                 [2,2,-1,-1,2,2,0,0,0,0],[2,-2,-1,1,2,-2,0,0,0,0],
                 [3,3,0,0,-1,-1,1,1,-1,-1],[3,3,0,0,-1,-1,-1,-1,1,1],
                 [3,-3,0,0,-1,1,I,-I,-I,I],[3,-3,0,0,-1,1,-I,I,I,-I]])
AGCharTab[(48,9)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                 [1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1],
                 [1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,-1,-1],
                 [1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1],
                 [1,1,1,1,-1,-1,-1,-1,1,-1,1,-1,1,-1,I,-I,I,-I],
                 [1,1,1,1,-1,-1,-1,-1,1,-1,1,-1,1,-1,-I,I,-I,I],
                 [1,1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,1,I,-I,-I,I],
                 [1,1,1,1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-I,I,I,-I],
                 [2,-2,2,-2,2,-2,2,-2,0,0,0,0,0,0,0,0,0,0],
                 [2,-2,2,-2,-2,2,-2,2,0,0,0,0,0,0,0,0,0,0],
                 [2,2,-1,-1,2,2,-1,-1,2,2,-1,-1,-1,-1,0,0,0,0],
                 [2,2,-1,-1,2,2,-1,-1,-2,-2,1,1,1,1,0,0,0,0],
                 [2,-2,-1,1,2,-2,-1,1,0,0,-I*sqrt(3),-I*sqrt(3),I*sqrt(3),I*sqrt(3),0,0,0,0],
                 [2,-2,-1,1,2,-2,-1,1,0,0,I*sqrt(3),I*sqrt(3),-I*sqrt(3),-I*sqrt(3),0,0,0,0],
                 [2,2,-1,-1,-2,-2,1,1,2,-2,-1,1,-1,1,0,0,0,0],
                 [2,2,-1,-1,-2,-2,1,1,-2,2,1,-1,1,-1,0,0,0,0],
                 [2,-2,-1,1,-2,2,1,-1,0,0,-I*sqrt(3),I*sqrt(3),I*sqrt(3),-I*sqrt(3),0,0,0,0],
                 [2,-2,-1,1,-2,2,1,-1,0,0,I*sqrt(3),-I*sqrt(3),-I*sqrt(3),I*sqrt(3),0,0,0,0]])
AGCharTab[(48,10)]=Matrix([[1,1,1,1,1,1,1,1],[1,1,1,1,1,-1,-1,-1],
                  [2,2,2,-1,-1,0,0,0],[2,-2,0,1,-1,sqrt(2),-sqrt(2),0],
                  [2,-2,0,1,-1,-sqrt(2),sqrt(2),0],[3,3,-1,0,0,1,1,-1],
                  [3,3,-1,0,0,-1,-1,1],[4,-4,0,-1,1,0,0,0]])
AGCharTab[(48,11)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                  [1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1],
                  [1,-omega.conjugate(),omega,-1,omega.conjugate(),-omega,1,omega,omega.conjugate(),1,-omega.conjugate(),omega,-1,omega.conjugate(),-omega,1,-omega.conjugate(),omega,-1,omega.conjugate(),-omega,1,omega,omega.conjugate(),1,-omega.conjugate(),omega,-1,omega.conjugate(),-omega],
                  [1,-omega.conjugate(),omega,-1,omega.conjugate(),-omega,1,omega,omega.conjugate(),-1,omega.conjugate(),-omega,1,-omega.conjugate(),omega,1,-omega.conjugate(),omega,-1,omega.conjugate(),-omega,1,omega,omega.conjugate(),-1,omega.conjugate(),-omega,1,-omega.conjugate(),omega],
                  [1,omega,omega.conjugate(),1,omega,omega.conjugate(),1,omega.conjugate(),omega,1,omega,omega.conjugate(),1,omega,omega.conjugate(),1,omega,omega.conjugate(),1,omega,omega.conjugate(),1,omega.conjugate(),omega,1,omega,omega.conjugate(),1,omega,omega.conjugate()],
                  [1,omega,omega.conjugate(),1,omega,omega.conjugate(),1,omega.conjugate(),omega,-1,-omega,-omega.conjugate(),-1,-omega,-omega.conjugate(),1,omega,omega.conjugate(),1,omega,omega.conjugate(),1,omega.conjugate(),omega,-1,-omega,-omega.conjugate(),-1,-omega,-omega.conjugate()],
                  [1,-1,1,-1,1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,-1,1,-1,1,-1],
                  [1,-1,1,-1,1,-1,1,1,1,-1,1,-1,1,-1,1,1,-1,1,-1,1,-1,1,1,1,-1,1,-1,1,-1,1],
                  [1,omega.conjugate(),omega,1,omega.conjugate(),omega,1,omega,omega.conjugate(),1,omega.conjugate(),omega,1,omega.conjugate(),omega,1,omega.conjugate(),omega,1,omega.conjugate(),omega,1,omega,omega.conjugate(),1,omega.conjugate(),omega,1,omega.conjugate(),omega],
                  [1,omega.conjugate(),omega,1,omega.conjugate(),omega,1,omega,omega.conjugate(),-1,-omega.conjugate(),-omega,-1,-omega.conjugate(),-omega,1,omega.conjugate(),omega,1,omega.conjugate(),omega,1,omega,omega.conjugate(),-1,-omega.conjugate(),-omega,-1,-omega.conjugate(),-omega],
                  [1,-omega,omega.conjugate(),-1,omega,-omega.conjugate(),1,omega.conjugate(),omega,1,-omega,omega.conjugate(),-1,omega,-omega.conjugate(),1,-omega,omega.conjugate(),-1,omega,-omega.conjugate(),1,omega.conjugate(),omega,1,-omega,omega.conjugate(),-1,omega,-omega.conjugate()],
                  [1,-omega,omega.conjugate(),-1,omega,-omega.conjugate(),1,omega.conjugate(),omega,-1,omega,-omega.conjugate(),1,-omega,omega.conjugate(),1,-omega,omega.conjugate(),-1,omega,-omega.conjugate(),1,omega.conjugate(),omega,-1,omega,-omega.conjugate(),1,-omega,omega.conjugate()],
                  [1,-I,-1,I,1,-I,-1,1,-1,1,-I,-1,I,1,-I,-1,I,1,-I,-1,I,1,-1,1,-1,I,1,-I,-1,I],
                  [1,-I,-1,I,1,-I,-1,1,-1,-1,I,1,-I,-1,I,-1,I,1,-I,-1,I,1,-1,1,1,-I,-1,I,1,-I],
                  [1,I*omega.conjugate(),-omega,-I,omega.conjugate(),I*omega,-1,omega,-omega.conjugate(),1,I*omega.conjugate(),-omega,-I,omega.conjugate(),I*omega,-1,-I*omega.conjugate(),omega,I,-omega.conjugate(),-I*omega,1,-omega,omega.conjugate(),-1,-I*omega.conjugate(),omega,I,-omega.conjugate(),-I*omega],
                  [1,I*omega.conjugate(),-omega,-I,omega.conjugate(),I*omega,-1,omega,-omega.conjugate(),-1,-I*omega.conjugate(),omega,I,-omega.conjugate(),-I*omega,-1,-I*omega.conjugate(),omega,I,-omega.conjugate(),-I*omega,1,-omega,omega.conjugate(),1,I*omega.conjugate(),-omega,-I,omega.conjugate(),I*omega],
                  [1,-I*omega,-omega.conjugate(),I,omega,-I*omega.conjugate(),-1,omega.conjugate(),-omega,1,-I*omega,-omega.conjugate(),I,omega,-I*omega.conjugate(),-1,I*omega,omega.conjugate(),-I,-omega,I*omega.conjugate(),1,-omega.conjugate(),omega,-1,I*omega,omega.conjugate(),-I,-omega,I*omega.conjugate()],
                  [1,-I*omega,-omega.conjugate(),I,omega,-I*omega.conjugate(),-1,omega.conjugate(),-omega,-1,I*omega,omega.conjugate(),-I,-omega,I*omega.conjugate(),-1,I*omega,omega.conjugate(),-I,-omega,I*omega.conjugate(),1,-omega.conjugate(),omega,1,-I*omega,-omega.conjugate(),I,omega,-I*omega.conjugate()],
                  [1,I,-1,-I,1,I,-1,1,-1,1,I,-1,-I,1,I,-1,-I,1,I,-1,-I,1,-1,1,-1,-I,1,I,-1,-I],
                  [1,I,-1,-I,1,I,-1,1,-1,-1,-I,1,I,-1,-I,-1,-I,1,I,-1,-I,1,-1,1,1,I,-1,-I,1,I],
                  [1,-I*omega.conjugate(),-omega,I,omega.conjugate(),-I*omega,-1,omega,-omega.conjugate(),1,-I*omega.conjugate(),-omega,I,omega.conjugate(),-I*omega,-1,I*omega.conjugate(),omega,-I,-omega.conjugate(),I*omega,1,-omega,omega.conjugate(),-1,I*omega.conjugate(),omega,-I,-omega.conjugate(),I*omega],
                  [1,-I*omega.conjugate(),-omega,I,omega.conjugate(),-I*omega,-1,omega,-omega.conjugate(),-1,I*omega.conjugate(),omega,-I,-omega.conjugate(),I*omega,-1,I*omega.conjugate(),omega,-I,-omega.conjugate(),I*omega,1,-omega,omega.conjugate(),1,-I*omega.conjugate(),-omega,I,omega.conjugate(),-I*omega],
                  [1,I*omega,-omega.conjugate(),-I,omega,I*omega.conjugate(),-1,omega.conjugate(),-omega,1,I*omega,-omega.conjugate(),-I,omega,I*omega.conjugate(),-1,-I*omega,omega.conjugate(),I,-omega,-I*omega.conjugate(),1,-omega.conjugate(),omega,-1,-I*omega,omega.conjugate(),I,-omega,-I*omega.conjugate()],
                  [1,I*omega,-omega.conjugate(),-I,omega,I*omega.conjugate(),-1,omega.conjugate(),-omega,-1,-I*omega,omega.conjugate(),I,-omega,-I*omega.conjugate(),-1,-I*omega,omega.conjugate(),I,-omega,-I*omega.conjugate(),1,-omega.conjugate(),omega,1,I*omega,-omega.conjugate(),-I,omega,I*omega.conjugate()],
                  [2,0,-2*omega.conjugate(),0,2*omega,0,-2,2*omega.conjugate(),-2*omega,0,0,0,0,0,0,2,0,-2*omega.conjugate(),0,2*omega,0,-2,2*omega.conjugate(),-2*omega,0,0,0,0,0,0],
                  [2,0,-2,0,2,0,-2,2,-2,0,0,0,0,0,0,2,0,-2,0,2,0,-2,2,-2,0,0,0,0,0,0],
                  [2,0,-2*omega,0,2*omega.conjugate(),0,-2,2*omega,-2*omega.conjugate(),0,0,0,0,0,0,2,0,-2*omega,0,2*omega.conjugate(),0,-2,2*omega,-2*omega.conjugate(),0,0,0,0,0,0],
                  [2,0,2*omega.conjugate(),0,2*omega,0,2,2*omega.conjugate(),2*omega,0,0,0,0,0,0,-2,0,-2*omega.conjugate(),0,-2*omega,0,-2,-2*omega.conjugate(),-2*omega,0,0,0,0,0,0],
                  [2,0,2,0,2,0,2,2,2,0,0,0,0,0,0,-2,0,-2,0,-2,0,-2,-2,-2,0,0,0,0,0,0],
                  [2,0,2*omega,0,2*omega.conjugate(),0,2,2*omega,2*omega.conjugate(),0,0,0,0,0,0,-2,0,-2*omega,0,-2*omega.conjugate(),0,-2,-2*omega,-2*omega.conjugate(),0,0,0,0,0,0]])


AGCharTab[(48,12)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                  [1,1,1,1,1,1,1,-1,-1,1,1,1,1,1,1,1,-1,-1],
                  [1,1,1,1,-1,-1,-1,1,-1,1,1,1,1,-1,-1,-1,1,-1],
                  [1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1],
                  [2,2,-1,-1,-2,1,1,0,0,2,2,-1,-1,-2,1,1,0,0],
                  [2,2,-1,-1,2,-1,-1,0,0,2,2,-1,-1,2,-1,-1,0,0],
                  [2,-2,-1,1,0,sqrt(3),-sqrt(3),0,0,2,-2,-1,1,0,sqrt(3),-sqrt(3),0,0],
                  [2,-2,-1,1,0,-sqrt(3),sqrt(3),0,0,2,-2,-1,1,0,-sqrt(3),sqrt(3),0,0],
                  [2,-2,2,-2,0,0,0,0,0,2,-2,2,-2,0,0,0,0,0],
                  [1,1,1,1,1,1,1,I,I,-1,-1,-1,-1,-1,-1,-1,-I,-I],
                  [1,1,1,1,1,1,1,-I,-I,-1,-1,-1,-1,-1,-1,-1,I,I],
                  [1,1,1,1,-1,-1,-1,I,-I,-1,-1,-1,-1,1,1,1,-I,I],
                  [1,1,1,1,-1,-1,-1,-I,I,-1,-1,-1,-1,1,1,1,I,-I],
                  [2,2,-1,-1,-2,1,1,0,0,-2,-2,1,1,2,-1,-1,0,0],
                  [2,2,-1,-1,2,-1,-1,0,0,-2,-2,1,1,-2,1,1,0,0],
                  [2,-2,-1,1,0,sqrt(3),-sqrt(3),0,0,-2,2,1,-1,0,-sqrt(3),sqrt(3),0,0],
                  [2,-2,-1,1,0,-sqrt(3),sqrt(3),0,0,-2,2,1,-1,0,sqrt(3),-sqrt(3),0,0],
                  [2,-2,2,-2,0,0,0,0,0,-2,2,-2,2,0,0,0,0,0]])

AGCharTab[(48,13)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                  [1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1],
                  [1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,-1,-1],
                  [1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1],
                  [1,1,-1,-1,1,1,-1,-1,I,-I,I,I,-I,-I,1,-1,I,-I],
                  [1,1,-1,-1,1,1,-1,-1,I,-I,I,I,-I,-I,-1,1,-I,I],
                  [1,1,-1,-1,1,1,-1,-1,-I,I,-I,-I,I,I,1,-1,-I,I],
                  [1,1,-1,-1,1,1,-1,-1,-I,I,-I,-I,I,I,-1,1,I,-I],
                  [2,2,2,2,-1,-1,-1,-1,2,2,-1,-1,-1,-1,0,0,0,0],
                  [2,2,2,2,-1,-1,-1,-1,-2,-2,1,1,1,1,0,0,0,0],
                  [2,2,-2,-2,-1,-1,1,1,2*I,-2*I,-I,-I,I,I,0,0,0,0],
                  [2,2,-2,-2,-1,-1,1,1,-2*I,2*I,I,I,-I,-I,0,0,0,0],
                  [2,-2,2,-2,-2,2,-2,2,0,0,0,0,0,0,0,0,0,0],
                  [2,-2,2,-2,1,-1,1,-1,0,0,I*sqrt(3),-I*sqrt(3),I*sqrt(3),-I*sqrt(3),0,0,0,0],
                  [2,-2,2,-2,1,-1,1,-1,0,0,-I*sqrt(3),I*sqrt(3),-I*sqrt(3),I*sqrt(3),0,0,0,0],
                  [2,-2,-2,2,-2,2,2,-2,0,0,0,0,0,0,0,0,0,0],
                  [2,-2,-2,2,1,-1,-1,1,0,0,sqrt(3),-sqrt(3),-sqrt(3),sqrt(3),0,0,0,0],
                  [2,-2,-2,2,1,-1,-1,1,0,0,-sqrt(3),sqrt(3),sqrt(3),-sqrt(3),0,0,0,0]])
AGCharTab[(48,14)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                  [1,1,1,1,1,1,1,-1,-1,1,1,1,1,1,1,1,-1,-1],
                  [1,1,1,1,-1,-1,-1,1,-1,1,1,1,1,-1,-1,-1,1,-1],
                  [1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1],
                  [1,1,-1,-1,I,I,-I,1,I,-1,-1,1,1,-I,-I,I,-1,-I],
                  [1,1,-1,-1,I,I,-I,-1,-I,-1,-1,1,1,-I,-I,I,1,I],
                  [1,1,-1,-1,-I,-I,I,1,-I,-1,-1,1,1,I,I,-I,-1,I],
                  [1,1,-1,-1,-I,-I,I,-1,I,-1,-1,1,1,I,I,-I,1,-I],
                  [2,2,-1,-1,2,-1,-1,0,0,2,2,-1,-1,2,-1,-1,0,0],
                  [2,2,-1,-1,-2,1,1,0,0,2,2,-1,-1,-2,1,1,0,0],
                  [2,2,1,1,2*I,-I,I,0,0,-2,-2,-1,-1,-2*I,I,-I,0,0],
                  [2,2,1,1,-2*I,I,-I,0,0,-2,-2,-1,-1,2*I,-I,I,0,0],
                  [2,-2,-1,1,0,I*sqrt(3),-I*sqrt(3),0,0,2,-2,-1,1,0,I*sqrt(3),-I*sqrt(3),0,0],
                  [2,-2,-1,1,0,-I*sqrt(3),I*sqrt(3),0,0,2,-2,-1,1,0,-I*sqrt(3),I*sqrt(3),0,0],
                  [2,-2,1,-1,0,sqrt(3),sqrt(3),0,0,-2,2,-1,1,0,-sqrt(3),-sqrt(3),0,0],
                  [2,-2,1,-1,0,-sqrt(3),-sqrt(3),0,0,-2,2,-1,1,0,sqrt(3),sqrt(3),0,0],
                  [2,-2,2,-2,0,0,0,0,0,2,-2,2,-2,0,0,0,0,0],
                  [2,-2,-2,2,0,0,0,0,0,-2,2,2,-2,0,0,0,0,0]])
AGCharTab[(48,15)]=kronecker_product(otimesG21,AGCharTab[(24,11)])
AGCharTab[(64,1)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                 [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1],
                 [1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,-1,1,-1],
                 [1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,1],
                 [1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,-1,-1],
                 [1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1],
                 [1,1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,-1,-1,1],
                 [1,1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,-1,1,1,-1],
                 [2,2,2,2,-2,-2,0,0,0,0,0,0,0,0,2,2,-2,-2,0,0,0,0],
                 [2,2,-2,-2,2,-2,2,2,-2,-2,0,0,0,0,0,0,0,0,0,0,0,0],
                 [2,2,-2,-2,2,-2,-2,-2,2,2,0,0,0,0,0,0,0,0,0,0,0,0],
                 [2,2,-2,-2,-2,2,0,0,0,0,2*I,2*I,-2*I,-2*I,0,0,0,0,0,0,0,0],
                 [2,2,-2,-2,-2,2,0,0,0,0,-2*I,-2*I,2*I,2*I,0,0,0,0,0,0,0,0],
                 [2,2,2,2,-2,-2,0,0,0,0,0,0,0,0,-2,-2,2,2,0,0,0,0],
                 [2,-2,-2,2,0,0,sqrt(2),-sqrt(2),-sqrt(2),sqrt(2),sqrt(2),-sqrt(2),-sqrt(2),sqrt(2),2,-2,0,0,0,0,0,0],
                 [2,-2,-2,2,0,0,sqrt(2),-sqrt(2),-sqrt(2),sqrt(2),-sqrt(2),sqrt(2),sqrt(2),-sqrt(2),-2,2,0,0,0,0,0,0],
                 [2,-2,-2,2,0,0,-sqrt(2),sqrt(2),sqrt(2),-sqrt(2),-sqrt(2),sqrt(2),sqrt(2),-sqrt(2),2,-2,0,0,0,0,0,0],
                 [2,-2,-2,2,0,0,-sqrt(2),sqrt(2),sqrt(2),-sqrt(2),sqrt(2),-sqrt(2),-sqrt(2),sqrt(2),-2,2,0,0,0,0,0,0],
                 [2,-2,2,-2,0,0,sqrt(2),-sqrt(2),sqrt(2),-sqrt(2),I*sqrt(2),-I*sqrt(2),I*sqrt(2),-I*sqrt(2),0,0,2*I,-2*I,0,0,0,0],
                 [2,-2,2,-2,0,0,sqrt(2),-sqrt(2),sqrt(2),-sqrt(2),-I*sqrt(2),I*sqrt(2),-I*sqrt(2),I*sqrt(2),0,0,-2*I,2*I,0,0,0,0],
                 [2,-2,2,-2,0,0,-sqrt(2),sqrt(2),-sqrt(2),sqrt(2),-I*sqrt(2),I*sqrt(2),-I*sqrt(2),I*sqrt(2),0,0,2*I,-2*I,0,0,0,0],
                 [2,-2,2,-2,0,0,-sqrt(2),sqrt(2),-sqrt(2),sqrt(2),I*sqrt(2),-I*sqrt(2),I*sqrt(2),-I*sqrt(2),0,0,-2*I,2*I,0,0,0,0]])

AGCharTab[(64,2)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                 [1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1],
                 [1,1,1,1,1,1,-1,-1,1,1,1,-1,-1,1,-1,-1,1,-1,-1],
                 [1,1,1,1,1,1,-1,-1,1,1,1,-1,-1,-1,1,1,-1,1,1],
                 [1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,-1,-1,-1],
                 [1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1],
                 [1,1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,-1,-1,-1,1,1],
                 [1,1,1,1,1,1,-1,-1,-1,-1,-1,1,1,-1,1,1,1,-1,-1],
                 [2,2,-2,-2,2,-2,2,-2,0,0,0,0,0,0,0,0,0,0,0],
                 [2,2,2,2,-2,-2,0,0,2,2,-2,0,0,0,0,0,0,0,0],
                 [2,2,-2,-2,-2,2,0,0,0,0,0,2,-2,0,0,0,0,0,0],
                 [2,2,-2,-2,2,-2,-2,2,0,0,0,0,0,0,0,0,0,0,0],
                 [2,2,2,2,-2,-2,0,0,-2,-2,2,0,0,0,0,0,0,0,0],
                 [2,2,-2,-2,-2,2,0,0,0,0,0,-2,2,0,0,0,0,0,0],
                 [2,-2,2,-2,0,0,0,0,-2,2,0,0,0,0,sqrt(2),-sqrt(2),0,-sqrt(2),sqrt(2)],
                 [2,-2,2,-2,0,0,0,0,2,-2,0,0,0,0,sqrt(2),-sqrt(2),0,sqrt(2),-sqrt(2)],
                 [2,-2,2,-2,0,0,0,0,-2,2,0,0,0,0,-sqrt(2),sqrt(2),0,sqrt(2),-sqrt(2)],
                 [2,-2,2,-2,0,0,0,0,2,-2,0,0,0,0,-sqrt(2),sqrt(2),0,-sqrt(2),sqrt(2)],
                 [4,-4,-4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]])
AGCharTab[(64,3)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                 [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1],
                 [1,1,-1,-1,1,1,-1,-1,-1,-1,1,1,-I,-I,I,I,-I,-I,I,I,I,I,-I,-I,-1,-1,1,1],
                 [1,1,-1,-1,1,1,-1,-1,-1,-1,1,1,-I,-I,I,I,-I,-I,I,I,-I,-I,I,I,1,1,-1,-1],
                 [1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1],
                 [1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1],
                 [1,1,-1,-1,1,1,-1,-1,-1,-1,1,1,I,I,-I,-I,I,I,-I,-I,-I,-I,I,I,-1,-1,1,1],
                 [1,1,-1,-1,1,1,-1,-1,-1,-1,1,1,I,I,-I,-I,I,I,-I,-I,I,I,-I,-I,1,1,-1,-1],
                 [1,-1,I,-I,-1,1,-I,I,-I,I,1,-1,I*theta,-I*theta,-theta,theta,-I*theta,I*theta,theta,-theta,I*theta,-I*theta,-theta,theta,1,-1,I,-I],
                 [1,-1,I,-I,-1,1,-I,I,-I,I,1,-1,I*theta,-I*theta,-theta,theta,-I*theta,I*theta,theta,-theta,-I*theta,I*theta,theta,-theta,-1,1,-I,I],
                 [1,-1,I,-I,-1,1,-I,I,-I,I,1,-1,-I*theta,I*theta,theta,-theta,I*theta,-I*theta,-theta,theta,I*theta,-I*theta,-theta,theta,-1,1,-I,I],
                 [1,-1,I,-I,-1,1,-I,I,-I,I,1,-1,-I*theta,I*theta,theta,-theta,I*theta,-I*theta,-theta,theta,-I*theta,I*theta,theta,-theta,1,-1,I,-I],
                 [1,-1,-I,I,-1,1,I,-I,I,-I,1,-1,theta,-theta,-I*theta,I*theta,-theta,theta,I*theta,-I*theta,theta,-theta,-I*theta,I*theta,1,-1,-I,I],
                 [1,-1,-I,I,-1,1,I,-I,I,-I,1,-1,theta,-theta,-I*theta,I*theta,-theta,theta,I*theta,-I*theta,-theta,theta,I*theta,-I*theta,-1,1,I,-I],
                 [1,-1,-I,I,-1,1,I,-I,I,-I,1,-1,-theta,theta,I*theta,-I*theta,theta,-theta,-I*theta,I*theta,theta,-theta,-I*theta,I*theta,-1,1,I,-I],
                 [1,-1,-I,I,-1,1,I,-I,I,-I,1,-1,-theta,theta,I*theta,-I*theta,theta,-theta,-I*theta,I*theta,-theta,theta,I*theta,-I*theta,1,-1,-I,I],
                 [2,2,2,2,2,2,2,2,-2,-2,-2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                 [2,2,-2,-2,2,2,-2,-2,2,2,-2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                 [2,-2,2*I,-2*I,-2,2,-2*I,2*I,2*I,-2*I,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                 [2,-2,-2*I,2*I,-2,2,2*I,-2*I,-2*I,2*I,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                 [2,2,-2*I,-2*I,-2,-2,2*I,2*I,0,0,0,0,1+I,1+I,1-I,1-I,-1-I,-1-I,-1+I,-1+I,0,0,0,0,0,0,0,0],
                 [2,2,-2*I,-2*I,-2,-2,2*I,2*I,0,0,0,0,-1-I,-1-I,-1+I,-1+I,1+I,1+I,1-I,1-I,0,0,0,0,0,0,0,0],
                 [2,2,2*I,2*I,-2,-2,-2*I,-2*I,0,0,0,0,1-I,1-I,1+I,1+I,-1+I,-1+I,-1-I,-1-I,0,0,0,0,0,0,0,0],
                 [2,2,2*I,2*I,-2,-2,-2*I,-2*I,0,0,0,0,-1+I,-1+I,-1-I,-1-I,1-I,1-I,1+I,1+I,0,0,0,0,0,0,0,0],
                 [2,-2,2,-2,2,-2,2,-2,0,0,0,0,sqrt(2),-sqrt(2),sqrt(2),-sqrt(2),sqrt(2),-sqrt(2),sqrt(2),-sqrt(2),0,0,0,0,0,0,0,0],
                 [2,-2,2,-2,2,-2,2,-2,0,0,0,0,-sqrt(2),sqrt(2),-sqrt(2),sqrt(2),-sqrt(2),sqrt(2),-sqrt(2),sqrt(2),0,0,0,0,0,0,0,0],
                 [2,-2,-2,2,2,-2,-2,2,0,0,0,0,I*sqrt(2),-I*sqrt(2),-I*sqrt(2),I*sqrt(2),I*sqrt(2),-I*sqrt(2),-I*sqrt(2),I*sqrt(2),0,0,0,0,0,0,0,0],
                 [2,-2,-2,2,2,-2,-2,2,0,0,0,0,-I*sqrt(2),I*sqrt(2),I*sqrt(2),-I*sqrt(2),-I*sqrt(2),I*sqrt(2),I*sqrt(2),-I*sqrt(2),0,0,0,0,0,0,0,0]])

AGCharTab[(64,4)]=kronecker_product(otimesG21,AGCharTab[(32,11)])

AGCharTab[(64,5)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                 [1,1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,1,-1,-1,-1],
                 [1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1],
                 [1,1,1,1,-1,-1,1,1,-1,1,1,1,1,-1,-1,1,-1,1,1],
                 [1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1],
                 [1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1],
                 [1,1,1,1,-1,-1,1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1],
                 [1,1,1,1,-1,-1,1,1,-1,1,1,-1,-1,1,1,-1,1,-1,-1],
                 [2,2,-2,-2,0,0,2,-2,0,0,0,0,0,2*I,-2*I,0,0,0,0],
                 [2,2,-2,-2,0,0,2,-2,0,0,0,0,0,-2*I,2*I,0,0,0,0],
                 [2,2,2,2,0,0,-2,-2,0,0,0,2,2,0,0,-2,0,0,0],
                 [2,2,2,2,0,0,-2,-2,0,0,0,-2,-2,0,0,2,0,0,0],
                 [2,2,-2,-2,2*I,-2*I,-2,2,0,0,0,0,0,0,0,0,0,0,0],
                 [2,2,-2,-2,-2*I,2*I,-2,2,0,0,0,0,0,0,0,0,0,0,0],
                 [2,-2,2,-2,0,0,0,0,0,sqrt(2),-sqrt(2),2,-2,0,0,0,0,sqrt(2),-sqrt(2)],
                 [2,-2,2,-2,0,0,0,0,0,sqrt(2),-sqrt(2),-2,2,0,0,0,0,-sqrt(2),sqrt(2)],
                 [2,-2,2,-2,0,0,0,0,0,-sqrt(2),sqrt(2),2,-2,0,0,0,0,-sqrt(2),sqrt(2)],
                 [2,-2,2,-2,0,0,0,0,0,-sqrt(2),sqrt(2),-2,2,0,0,0,0,sqrt(2),-sqrt(2)],
                 [4,-4,-4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]])
AGCharTab[(96,1)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                 [1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1],
                 [1,-1,1,-1,-1,1,1,-1,1,-1,I,-I,I,-I,I,-I],
                 [1,-1,1,-1,-1,1,1,-1,1,-1,-I,I,-I,I,-I,I],
                 [2,2,2,2,2,2,-1,-1,-1,-1,0,0,0,0,0,0],
                 [2,-2,2,-2,-2,2,-1,1,-1,1,0,0,0,0,0,0],
                 [2,2*I,-2,-2*I,0,0,-1,-I,1,I,1-I,1+I,-1+I,-1-I,0,0],
                 [2,2*I,-2,-2*I,0,0,-1,-I,1,I,-1+I,-1-I,1-I,1+I,0,0],
                 [2,-2*I,-2,2*I,0,0,-1,I,1,-I,1+I,1-I,-1-I,-1+I,0,0],
                 [2,-2*I,-2,2*I,0,0,-1,I,1,-I,-1-I,-1+I,1+I,1-I,0,0],
                 [3,3,3,3,-1,-1,0,0,0,0,1,1,1,1,-1,-1],
                 [3,3,3,3,-1,-1,0,0,0,0,-1,-1,-1,-1,1,1],
                 [3,-3,3,-3,1,-1,0,0,0,0,I,-I,I,-I,-I,I],
                 [3,-3,3,-3,1,-1,0,0,0,0,-I,I,-I,I,I,-I],
                 [4,4*I,-4,-4*I,0,0,1,I,-1,-I,0,0,0,0,0,0],
                 [4,-4*I,-4,4*I,0,0,1,-I,-1,I,0,0,0,0,0,0]])
AGCharTab[(96,2)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,-1,-1,1,1,1,1,-1,-1,-1,-1],
                 [1,1,1,1,1,1,1,1,-1,-1,1,1,-1,-1],[1,1,1,1,-1,-1,1,1,-1,-1,-1,-1,1,1],
                 [2,2,2,2,2,2,-1,-1,0,0,-1,-1,0,0],[2,2,2,2,-2,-2,-1,-1,0,0,1,1,0,0],
                 [2,-2,2,-2,0,0,2,-2,0,0,0,0,0,0],[2,-2,2,-2,0,0,-1,1,0,0,I*sqrt(3),-I*sqrt(3),0,0],
                 [2,-2,2,-2,0,0,-1,1,0,0,-I*sqrt(3),I*sqrt(3),0,0],[3,3,-1,-1,3,-1,0,0,1,-1,0,0,1,-1],
                 [3,3,-1,-1,-3,1,0,0,1,-1,0,0,-1,1],[3,3,-1,-1,3,-1,0,0,-1,1,0,0,-1,1],
                 [3,3,-1,-1,-3,1,0,0,-1,1,0,0,1,-1],[6,-6,-2,2,0,0,0,0,0,0,0,0,0,0]])
AGCharTab[(96,3)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,-1,-1,1,1,1,1,-1,-1,-1,-1],
                 [1,1,1,1,1,1,1,1,-1,-1,1,1,-1,-1],[1,1,1,1,-1,-1,1,1,-1,-1,-1,-1,1,1],
                 [2,2,2,2,2,2,-1,-1,0,0,-1,-1,0,0],[2,2,2,2,-2,-2,-1,-1,0,0,1,1,0,0],
                 [2,-2,2,-2,0,0,2,-2,0,0,0,0,0,0],[2,-2,2,-2,0,0,-1,1,0,0,I*sqrt(3),-I*sqrt(3),0,0],
                 [2,-2,2,-2,0,0,-1,1,0,0,-I*sqrt(3),I*sqrt(3),0,0],[3,3,-1,-1,3,-1,0,0,1,-1,0,0,1,-1],
                 [3,3,-1,-1,-3,1,0,0,1,-1,0,0,-1,1],[3,3,-1,-1,3,-1,0,0,-1,1,0,0,-1,1],
                 [3,3,-1,-1,-3,1,0,0,-1,1,0,0,1,-1],[6,-6,-2,2,0,0,0,0,0,0,0,0,0,0]])

AGCharTab[(96,4)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,-1,-1,1,1,1,1,-1,-1,-1,-1],
                 [1,1,1,1,1,1,1,1,-1,-1,1,1,-1,-1],[1,1,1,1,-1,-1,1,1,-1,-1,-1,-1,1,1],
                 [2,2,2,2,2,2,-1,-1,0,0,-1,-1,0,0],[2,2,2,2,-2,-2,-1,-1,0,0,1,1,0,0],
                 [2,-2,2,-2,0,0,2,-2,0,0,0,0,0,0],[2,-2,2,-2,0,0,-1,1,0,0,I*sqrt(3),-I*sqrt(3),0,0],
                 [2,-2,2,-2,0,0,-1,1,0,0,-I*sqrt(3),I*sqrt(3),0,0],[3,3,-1,-1,3,-1,0,0,1,-1,0,0,1,-1],
                 [3,3,-1,-1,-3,1,0,0,1,-1,0,0,-1,1],[3,3,-1,-1,3,-1,0,0,-1,1,0,0,-1,1],
                 [3,3,-1,-1,-3,1,0,0,-1,1,0,0,1,-1],[6,-6,-2,2,0,0,0,0,0,0,0,0,0,0]])

AGCharTab[(96,5)]=kronecker_product(otimesG41,AGCharTab[(24,9)])


AGCharTab[(96,6)]=kronecker_product(otimesG21,AGCharTab[(48,4)])

AGCharTab[(96,7)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                 [1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1],
                 [1,-1,-1,1,-1,1,-1,1,1,-1,I,-I,-I,I,-I,I],
                 [1,-1,-1,1,-1,1,-1,1,1,-1,-I,I,I,-I,I,-I],
                 [2,2,2,2,2,2,-1,-1,-1,-1,0,0,0,0,0,0],
                 [2,2,-2,-2,0,0,-1,-1,1,1,-I*sqrt(2),-I*sqrt(2),I*sqrt(2),I*sqrt(2),0,0],
                 [2,2,-2,-2,0,0,-1,-1,1,1,I*sqrt(2),I*sqrt(2),-I*sqrt(2),-I*sqrt(2),0,0],
                 [2,-2,-2,2,-2,2,1,-1,-1,1,0,0,0,0,0,0],
                 [2,-2,2,-2,0,0,1,-1,1,-1,sqrt(2),-sqrt(2),sqrt(2),-sqrt(2),0,0],
                 [2,-2,2,-2,0,0,1,-1,1,-1,-sqrt(2),sqrt(2),-sqrt(2),sqrt(2),0,0],
                 [3,3,3,3,-1,-1,0,0,0,0,1,1,1,1,-1,-1],
                 [3,3,3,3,-1,-1,0,0,0,0,-1,-1,-1,-1,1,1],
                 [3,-3,-3,3,1,-1,0,0,0,0,I,-I,-I,I,I,-I],
                 [3,-3,-3,3,1,-1,0,0,0,0,-I,I,I,-I,-I,I],
                 [4,4,-4,-4,0,0,1,1,-1,-1,0,0,0,0,0,0],
                 [4,-4,4,-4,0,0,-1,1,-1,1,0,0,0,0,0,0]])

AGCharTab[(96,8)]=kronecker_product(otimesG21,AGCharTab[(48,10)])

AGCharTab[(96,9)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                 [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
                 [1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1],
                 [1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1],
                 [1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,-1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,-1],
                 [1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,-1,1],
                 [1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,1],
                 [1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,-1],
                 [2,2,2,2,-1,-1,-1,-1,2,2,-1,-1,-1,-1,0,0,2,2,-1,-1,-1,-1,2,2,-1,-1,-1,-1,0,0],
                 [2,2,2,2,-1,-1,-1,-1,2,2,-1,-1,-1,-1,0,0,-2,-2,1,1,1,1,-2,-2,1,1,1,1,0,0],
                 [2,2,2,2,-1,-1,-1,-1,-2,-2,1,1,1,1,0,0,2,2,-1,-1,-1,-1,-2,-2,1,1,1,1,0,0],
                 [2,2,2,2,-1,-1,-1,-1,-2,-2,1,1,1,1,0,0,-2,-2,1,1,1,1,2,2,-1,-1,-1,-1,0,0],
                 [2,2,-2,-2,2,2,-2,-2,0,0,0,0,0,0,0,0,2,-2,2,2,-2,-2,0,0,0,0,0,0,0,0],
                 [2,2,-2,-2,2,2,-2,-2,0,0,0,0,0,0,0,0,-2,2,-2,-2,2,2,0,0,0,0,0,0,0,0],
                 [2,2,-2,-2,-1,-1,1,1,0,0,I*sqrt(3),I*sqrt(3),-I*sqrt(3),-I*sqrt(3),0,0,2,-2,-1,-1,1,1,0,0,I*sqrt(3),I*sqrt(3),-I*sqrt(3),-I*sqrt(3),0,0],
                 [2,2,-2,-2,-1,-1,1,1,0,0,I*sqrt(3),I*sqrt(3),-I*sqrt(3),-I*sqrt(3),0,0,-2,2,1,1,-1,-1,0,0,-I*sqrt(3),-I*sqrt(3),I*sqrt(3),I*sqrt(3),0,0],
                 [2,2,-2,-2,-1,-1,1,1,0,0,-I*sqrt(3),-I*sqrt(3),I*sqrt(3),I*sqrt(3),0,0,2,-2,-1,-1,1,1,0,0,-I*sqrt(3),-I*sqrt(3),I*sqrt(3),I*sqrt(3),0,0],
                 [2,2,-2,-2,-1,-1,1,1,0,0,-I*sqrt(3),-I*sqrt(3),I*sqrt(3),I*sqrt(3),0,0,-2,2,1,1,-1,-1,0,0,I*sqrt(3),I*sqrt(3),-I*sqrt(3),-I*sqrt(3),0,0],
                 [2,-2,2,-2,-2,2,-2,2,I*2,-I*2,I*2,-I*2,-I*2,I*2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                 [2,-2,2,-2,-2,2,-2,2,-I*2,I*2,-I*2,I*2,I*2,-I*2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                 [2,-2,2,-2,1,-1,1,-1,I*2,-I*2,-I,I,I,-I,0,0,0,0,sqrt(3),-sqrt(3),sqrt(3),-sqrt(3),0,0,I*sqrt(3),-I*sqrt(3),I*sqrt(3),-I*sqrt(3),0,0],
                 [2,-2,2,-2,1,-1,1,-1,I*2,-I*2,-I,I,I,-I,0,0,0,0,-sqrt(3),sqrt(3),-sqrt(3),sqrt(3),0,0,-I*sqrt(3),I*sqrt(3),-I*sqrt(3),I*sqrt(3),0,0],
                 [2,-2,2,-2,1,-1,1,-1,-I*2,I*2,I,-I,-I,I,0,0,0,0,sqrt(3),-sqrt(3),sqrt(3),-sqrt(3),0,0,-I*sqrt(3),I*sqrt(3),-I*sqrt(3),I*sqrt(3),0,0],
                 [2,-2,2,-2,1,-1,1,-1,-I*2,I*2,I,-I,-I,I,0,0,0,0,-sqrt(3),sqrt(3),-sqrt(3),sqrt(3),0,0,I*sqrt(3),-I*sqrt(3),I*sqrt(3),-I*sqrt(3),0,0],
                 [2,-2,-2,2,-2,2,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,-2,2,-2,-2,2,0,0],
                 [2,-2,-2,2,-2,2,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,2,-2,2,2,-2,0,0],
                 [2,-2,-2,2,1,-1,-1,1,0,0,sqrt(3),-sqrt(3),sqrt(3),-sqrt(3),0,0,0,0,-sqrt(3),sqrt(3),sqrt(3),-sqrt(3),2,-2,-1,1,1,-1,0,0],
                 [2,-2,-2,2,1,-1,-1,1,0,0,sqrt(3),-sqrt(3),sqrt(3),-sqrt(3),0,0,0,0,sqrt(3),-sqrt(3),-sqrt(3),sqrt(3),-2,2,1,-1,-1,1,0,0],
                 [2,-2,-2,2,1,-1,-1,1,0,0,-sqrt(3),sqrt(3),-sqrt(3),sqrt(3),0,0,0,0,sqrt(3),-sqrt(3),-sqrt(3),sqrt(3),2,-2,-1,1,1,-1,0,0],
                 [2,-2,-2,2,1,-1,-1,1,0,0,-sqrt(3),sqrt(3),-sqrt(3),sqrt(3),0,0,0,0,-sqrt(3),sqrt(3),sqrt(3),-sqrt(3),-2,2,1,-1,-1,1,0,0]])

AGCharTab[(96,10)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                  [1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
                  [1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1],
                  [1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1],
                  [1,1,1,1,-1,1,1,1,1,-1,-1,1,-1,-1,1,1,-1,1,1,-1,-1,-1,-1,1],
                  [1,1,1,1,-1,1,1,1,1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,1,1,1,1,-1],
                  [1,1,1,1,-1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,-1,1,1,-1],
                  [1,1,1,1,-1,1,1,1,1,-1,-1,-1,1,1,-1,-1,1,-1,-1,1,1,-1,-1,1],
                  [2,2,2,2,2,-1,-1,-1,-1,-1,-1,0,0,0,2,2,2,-1,-1,-1,-1,0,0,0],
                  [2,2,2,2,2,-1,-1,-1,-1,-1,-1,0,0,0,-2,-2,-2,1,1,1,1,0,0,0],
                  [2,2,2,2,-2,-1,-1,-1,-1,1,1,0,0,0,2,2,-2,-1,-1,1,1,0,0,0],
                  [2,2,2,2,-2,-1,-1,-1,-1,1,1,0,0,0,-2,-2,2,1,1,-1,-1,0,0,0],
                  [2,-2,2,-2,0,2,2,-2,-2,0,0,0,2,-2,0,0,0,0,0,0,0,0,0,0],
                  [2,-2,2,-2,0,2,2,-2,-2,0,0,0,-2,2,0,0,0,0,0,0,0,0,0,0],
                  [2,2,-2,-2,0,2,-2,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,2*I,-2*I,0],
                  [2,2,-2,-2,0,2,-2,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,-2*I,2*I,0],
                  [2,-2,-2,2,0,2,-2,-2,2,0,0,0,0,0,2,-2,0,-2,2,0,0,0,0,0],
                  [2,-2,-2,2,0,2,-2,-2,2,0,0,0,0,0,-2,2,0,2,-2,0,0,0,0,0],
                  [2,-2,-2,2,0,-1,1,1,-1,sqrt(3),-sqrt(3),0,0,0,2,-2,0,1,-1,sqrt(3),-sqrt(3),0,0,0],
                  [2,-2,-2,2,0,-1,1,1,-1,sqrt(3),-sqrt(3),0,0,0,-2,2,0,-1,1,-sqrt(3),sqrt(3),0,0,0],
                  [2,-2,-2,2,0,-1,1,1,-1,-sqrt(3),sqrt(3),0,0,0,2,-2,0,1,-1,-sqrt(3),sqrt(3),0,0,0],
                  [2,-2,-2,2,0,-1,1,1,-1,-sqrt(3),sqrt(3),0,0,0,-2,2,0,-1,1,sqrt(3),-sqrt(3),0,0,0],
                  [4,-4,4,-4,0,-2,-2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                  [4,4,-4,-4,0,-2,2,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]])

AGCharTab[(192,1)]=Matrix([
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1],
                    [1,-1,1,-1,-1,1,1,-1,1,-1,I,-I,I,-I,I,-I,1,-1,1,-1,-1,1,1,-1,1,-1,I,-I,I,-I,I,-I],
                    [1,-1,1,-1,-1,1,1,-1,1,-1,-I,I,-I,I,-I,I,1,-1,1,-1,-1,1,1,-1,1,-1,-I,I,-I,I,-I,I],
                    [2,2,2,2,2,2,-1,-1,-1,-1,0,0,0,0,0,0,2,2,2,2,2,2,-1,-1,-1,-1,0,0,0,0,0,0],
                    [2,-2,2,-2,-2,2,-1,1,-1,1,0,0,0,0,0,0,2,-2,2,-2,-2,2,-1,1,-1,1,0,0,0,0,0,0],
                    [2,2*I,-2,-2*I,0,0,-1,-I,1,I,sqrt(2)*conjugate(theta),sqrt(2)*theta,-sqrt(2)*conjugate(theta),-sqrt(2)*theta,0,0,2,2*I,-2,-2*I,0,0,-1,-I,1,I,sqrt(2)*conjugate(theta),sqrt(2)*theta,-sqrt(2)*conjugate(theta),-sqrt(2)*theta,0,0],
                    [2,2*I,-2,-2*I,0,0,-1,-I,1,I,-sqrt(2)*conjugate(theta),-sqrt(2)*theta,sqrt(2)*conjugate(theta),sqrt(2)*theta,0,0,2,2*I,-2,-2*I,0,0,-1,-I,1,I,-sqrt(2)*conjugate(theta),-sqrt(2)*theta,sqrt(2)*conjugate(theta),sqrt(2)*theta,0,0],
                    [2,-2*I,-2,2*I,0,0,-1,I,1,-I,sqrt(2)*theta,sqrt(2)*conjugate(theta),-sqrt(2)*theta,-sqrt(2)*conjugate(theta),0,0,2,-2*I,-2,2*I,0,0,-1,I,1,-I,sqrt(2)*theta,sqrt(2)*conjugate(theta),-sqrt(2)*theta,-sqrt(2)*conjugate(theta),0,0],
                    [2,-2*I,-2,2*I,0,0,-1,I,1,-I,-sqrt(2)*theta,-sqrt(2)*conjugate(theta),sqrt(2)*theta,sqrt(2)*conjugate(theta),0,0,2,-2*I,-2,2*I,0,0,-1,I,1,-I,-sqrt(2)*theta,-sqrt(2)*conjugate(theta),sqrt(2)*theta,sqrt(2)*conjugate(theta),0,0],
                    [3,3,3,3,-1,-1,0,0,0,0,1,1,1,1,-1,-1,3,3,3,3,-1,-1,0,0,0,0,1,1,1,1,-1,-1],
                    [3,3,3,3,-1,-1,0,0,0,0,-1,-1,-1,-1,1,1,3,3,3,3,-1,-1,0,0,0,0,-1,-1,-1,-1,1,1],
                    [3,-3,3,-3,1,-1,0,0,0,0,I,-I,I,-I,-I,I,3,-3,3,-3,1,-1,0,0,0,0,I,-I,I,-I,-I,I],
                    [3,-3,3,-3,1,-1,0,0,0,0,-I,I,-I,I,I,-I,3,-3,3,-3,1,-1,0,0,0,0,-I,I,-I,I,I,-I],
                    [4,4*I,-4,-4*I,0,0,1,I,-1,-I,0,0,0,0,0,0,4,4*I,-4,-4*I,0,0,1,I,-1,-I,0,0,0,0,0,0],
                    [4,-4*I,-4,4*I,0,0,1,-I,-1,I,0,0,0,0,0,0,4,-4*I,-4,4*I,0,0,1,-I,-1,I,0,0,0,0,0,0],
                    [1,-I,-1,I,-I,-1,-1,I,1,-I,-theta,-conjugate(theta),theta,conjugate(theta),theta,conjugate(theta),-1,I,1,-I,I,1,1,-I,-1,I,theta,conjugate(theta),-theta,-conjugate(theta),-theta,-conjugate(theta)],
                    [1,-I,-1,I,-I,-1,-1,I,1,-I,theta,conjugate(theta),-theta,-conjugate(theta),-theta,-conjugate(theta),-1,I,1,-I,I,1,1,-I,-1,I,-theta,-conjugate(theta),theta,conjugate(theta),theta,conjugate(theta)],
                    [1,I,-1,-I,I,-1,-1,-I,1,I,-conjugate(theta),-theta,conjugate(theta),theta,conjugate(theta),theta,-1,-I,1,I,-I,1,1,I,-1,-I,conjugate(theta),theta,-conjugate(theta),-theta,-conjugate(theta),-theta],
                    [1,I,-1,-I,I,-1,-1,-I,1,I,conjugate(theta),theta,-conjugate(theta),-theta,-conjugate(theta),-theta,-1,-I,1,I,-I,1,1,I,-1,-I,-conjugate(theta),-theta,conjugate(theta),theta,conjugate(theta),theta],
                    [2,-2*I,-2,2*I,-2*I,-2,1,-I,-1,I,0,0,0,0,0,0,-2,2*I,2,-2*I,2*I,2,-1,I,1,-I,0,0,0,0,0,0],
                    [2,2*I,-2,-2*I,2*I,-2,1,I,-1,-I,0,0,0,0,0,0,-2,-2*I,2,2*I,-2*I,2,-1,-I,1,I,0,0,0,0,0,0],
                    [2,-2,2,-2,0,0,1,-1,1,-1,-I*sqrt(2),I*sqrt(2),-I*sqrt(2),I*sqrt(2),0,0,-2,2,-2,2,0,0,-1,1,-1,1,I*sqrt(2),-I*sqrt(2),I*sqrt(2),-I*sqrt(2),0,0],
                    [2,-2,2,-2,0,0,1,-1,1,-1,I*sqrt(2),-I*sqrt(2),I*sqrt(2),-I*sqrt(2),0,0,-2,2,-2,2,0,0,-1,1,-1,1,-I*sqrt(2),I*sqrt(2),-I*sqrt(2),I*sqrt(2),0,0],
                    [2,2,2,2,0,0,1,1,1,1,-sqrt(2),-sqrt(2),-sqrt(2),-sqrt(2),0,0,-2,-2,-2,-2,0,0,-1,-1,-1,-1,sqrt(2),sqrt(2),sqrt(2),sqrt(2),0,0],
                    [2,2,2,2,0,0,1,1,1,1,sqrt(2),sqrt(2),sqrt(2),sqrt(2),0,0,-2,-2,-2,-2,0,0,-1,-1,-1,-1,-sqrt(2),-sqrt(2),-sqrt(2),-sqrt(2),0,0],
                    [3,-3*I,-3,3*I,I,1,0,0,0,0,-theta,-conjugate(theta),theta,conjugate(theta),-theta,-conjugate(theta),-3,3*I,3,-3*I,-I,-1,0,0,0,0,theta,conjugate(theta),-theta,-conjugate(theta),theta,conjugate(theta)],
                    [3,-3*I,-3,3*I,I,1,0,0,0,0,theta,conjugate(theta),-theta,-conjugate(theta),theta,conjugate(theta),-3,3*I,3,-3*I,-I,-1,0,0,0,0,-theta,-conjugate(theta),theta,conjugate(theta),-theta,-conjugate(theta)],
                    [3,3*I,-3,-3*I,-I,1,0,0,0,0,-conjugate(theta),-theta,conjugate(theta),theta,-conjugate(theta),-theta,-3,-3*I,3,3*I,I,-1,0,0,0,0,conjugate(theta),theta,-conjugate(theta),-theta,conjugate(theta),theta],
                    [3,3*I,-3,-3*I,-I,1,0,0,0,0,conjugate(theta),theta,-conjugate(theta),-theta,conjugate(theta),theta,-3,-3*I,3,3*I,I,-1,0,0,0,0,-conjugate(theta),-theta,conjugate(theta),theta,-conjugate(theta),-theta],
                    [4,-4,4,-4,0,0,-1,1,-1,1,0,0,0,0,0,0,-4,4,-4,4,0,0,1,-1,1,-1,0,0,0,0,0,0],
                    [4,4,4,4,0,0,-1,-1,-1,-1,0,0,0,0,0,0,-4,-4,-4,-4,0,0,1,1,1,1,0,0,0,0,0,0]])


AGCharTab[(192,2)]=Matrix([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                  [1,1,1,1,1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1],
                  [1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,-1,-1],
                  [1,1,1,1,1,1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,1,1,1],
                  [2,2,2,2,2,2,2,2,2,-1,-1,-1,-1,0,0,0,-1,-1,-1,-1,0,0,0],
                  [2,2,2,2,2,2,-2,-2,-2,-1,-1,-1,-1,0,0,0,1,1,1,1,0,0,0],
                  [2,2,-2,-2,2,-2,0,0,0,2,2,-2,-2,0,0,0,0,0,0,0,0,0,0],
                  [2,2,-2,-2,2,-2,0,0,0,-1,-1,1,1,0,0,0,I*sqrt(3),I*sqrt(3),-I*sqrt(3),-I*sqrt(3),0,0,0],
                  [2,2,-2,-2,2,-2,0,0,0,-1,-1,1,1,0,0,0,-I*sqrt(3),-I*sqrt(3),I*sqrt(3),I*sqrt(3),0,0,0],
                  [2,-2,2,-2,0,0,2,-2,0,1,-1,1,-1,sqrt(2),-sqrt(2),0,1,-1,1,-1,sqrt(2),-sqrt(2),0],
                  [2,-2,2,-2,0,0,-2,2,0,1,-1,1,-1,sqrt(2),-sqrt(2),0,-1,1,-1,1,-sqrt(2),sqrt(2),0],
                  [2,-2,2,-2,0,0,2,-2,0,1,-1,1,-1,-sqrt(2),sqrt(2),0,1,-1,1,-1,-sqrt(2),sqrt(2),0],
                  [2,-2,2,-2,0,0,-2,2,0,1,-1,1,-1,-sqrt(2),sqrt(2),0,-1,1,-1,1,sqrt(2),-sqrt(2),0],
                  [3,3,3,3,-1,-1,3,3,-1,0,0,0,0,1,1,-1,0,0,0,0,1,1,-1],
                  [3,3,3,3,-1,-1,-3,-3,1,0,0,0,0,1,1,-1,0,0,0,0,-1,-1,1],
                  [3,3,3,3,-1,-1,3,3,-1,0,0,0,0,-1,-1,1,0,0,0,0,-1,-1,1],
                  [3,3,3,3,-1,-1,-3,-3,1,0,0,0,0,-1,-1,1,0,0,0,0,1,1,-1],
                  [4,-4,-4,4,0,0,0,0,0,2,-2,-2,2,0,0,0,0,0,0,0,0,0,0],
                  [4,-4,4,-4,0,0,4,-4,0,-1,1,-1,1,0,0,0,-1,1,-1,1,0,0,0],
                  [4,-4,4,-4,0,0,-4,4,0,-1,1,-1,1,0,0,0,1,-1,1,-1,0,0,0],
                  [4,-4,-4,4,0,0,0,0,0,-1,1,1,-1,0,0,0,I*sqrt(3),-I*sqrt(3),-I*sqrt(3),I*sqrt(3),0,0,0],
                  [4,-4,-4,4,0,0,0,0,0,-1,1,1,-1,0,0,0,-I*sqrt(3),I*sqrt(3),I*sqrt(3),-I*sqrt(3),0,0,0],
                  [6,6,-6,-6,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]])


AGClasses={}
AGClasses[(1,1)]=[[[0]]];

AGClasses[(2,1)]=[[[0]],[[1]]];

AGClasses[(3,1)]=[[[0]],[[1]],[[2]]];

AGClasses[(4,1)]=[[[0]],[[1]],[[2]],[[3]]];

AGClasses[(4,2)]=[[[0,0]],[[1,0]],[[0,1]],[[1,1]]];

AGClasses[(6,1)]=[[[0]],[[1]],[[2]],[[3]],[[4]],[[5]]];

AGClasses[(6,2)]=[[[0,0]],[[1,0],[2,0]],[[0,1],[1,1],[2,1]]];

AGClasses[(8,1)]=[[[0]],[[1]],[[2]],[[3]],[[4]],[[5]],[[6]],[[7]]];

AGClasses[(8,2)]=[[[0,0]],[[1,0]],[[2,0]],[[3,0]],[[0,1]],[[1,1]],[[2,1]],[[3,1]]];

AGClasses[(8,3)]=[[[0,0,0]],[[1,0,0]],[[0,1,0]],[[1,1,0]],[[0,0,1]],[[1,0,1]],[[0,1,1]],[[1,1,1]]];

AGClasses[(8,4)]=[[[0,0]],[[2,0]],[[1,0],[3,0]],[[0,1],[2,1]],[[1,1],[3,1]]];

AGClasses[(8,5)]=[[[0,0]],[[2,0]],[[1,0],[3,0]],[[0,1],[2,1]],[[1,1],[3,1]]];

AGClasses[(12,1)]=[[[0]],[[1]],[[2]],[[3]],[[4]],[[5]],[[6]],[[7]],[[8]],[[9]],[[10]],[[11]]];

AGClasses[(12,2)]=[[[0,0,0]],[[1,0,0]],[[2,0,0]],[[0,1,0]],[[1,1,0]],[[2,1,0]],
                 [[0,0,1]],[[1,0,1]],[[2,0,1]],[[0,1,1]],[[1,1,1]],[[2,1,1]]];

AGClasses[(12,3)]=[[[0,0]],[[3,0]],[[1,0],[5,0]],[[2,0],[4,0]],[[0,1],[2,1],[4,1]],[[1,1],[3,1],[5,1]]];

AGClasses[(12,4)]=[[[0,0]],[[3,0]],[[1,0],[5,0]],[[2,0],[4,0]],[[0,1],[2,1],[4,1]],[[1,1],[3,1],[5,1]]];

AGClasses[(12,5)]=[[[0,0,0]],[[0,1,0],[0,0,1],[0,1,1]],[[1,0,0],[1,1,0],[1,0,1],[1,1,1]],
                 [[2,0,0],[2,1,0],[2,0,1],[2,1,1]]];

AGClasses[(12,6)]=[[[0,0]],[[1,0]],[[2,0]],[[3,0]],[[4,0]],[[5,0]],[[0,1]],
                 [[1,1]],[[2,1]],[[3,1]],[[4,1]],[[5,1]]];

AGClasses[(16,1)]=[[[s - 1] for s in range(1, 17)]];

AGClasses[(16,2)]=[[[0,0]],[[1,0]],[[2,0]],[[3,0]],[[4,0]],[[5,0]],[[6,0]],[[7,0]],
                 [[0,1]],[[1,1]],[[2,1]],[[3,1]],[[4,1]],[[5,1]],[[6,1]],[[7,1]]];

AGClasses[(16,3)]=[[[0,0]],[[1,0]],[[2,0]],[[3,0]],[[0,1]],[[1,1]],[[2,1]],[[3,1]],
                 [[0,2]],[[1,2]],[[2,2]],[[3,2]],[[0,3]],[[1,3]],[[2,3]],[[3,3]]];

AGClasses[(16,4)]=[[[0,0,0]],[[1,0,0]],[[2,0,0]],[[3,0,0]],
                 [[0,1,0]],[[1,1,0]],[[2,1,0]],[[3,1,0]],
                 [[0,0,1]],[[1,0,1]],[[2,0,1]],[[3,0,1]],
                 [[0,1,1]],[[1,1,1]],[[2,1,1]],[[3,1,1]]];

AGClasses[(16,5)]=[[[0,0,0,0]],[[1,0,0,0]],[[0,1,0,0]],[[1,1,0,0]],
                 [[0,0,1,0]],[[1,0,1,0]],[[0,1,1,0]],[[1,1,1,0]],
                 [[0,0,0,1]],[[1,0,0,1]],[[0,1,0,1]],[[1,1,0,1]],
                 [[0,0,1,1]],[[1,0,1,1]],[[0,1,1,1]],[[1,1,1,1]]];

AGClasses[(16,6)]=[[[0,0]],[[2,0]],[[4,0]],[[6,0]],[[1,0],[5,0]],[[3,0],[7,0]],
                  [[0,1],[4,1]],[[2,1],[6,1]],[[1,1],[5,1]],[[3,1],[7,1]]];


AGClasses[(16,7)]=[[[0,0,0]],[[1,0,0]],[[2,0,0]],[[3,0,0]],
                 [[0,1,0],[2,1,0]],[[1,1,0],[3,1,0]],
                 [[0,0,1],[2,0,1]],[[1,0,1],[3,0,1]],
                 [[0,1,1],[2,1,1]],[[1,1,1],[3,1,1]]];

AGClasses[(16,8)]=[[[0,0]],[[2,0]],[[0,2]],[[2,2]],[[1,0],[3,0]],[[1,2],[3,2]],
                  [[0,1],[2,1]],[[1,1],[3,1]],[[0,3],[2,3]],[[1,3],[3,3]]];

AGClasses[(16,9)]=[[[0,0,0]],[[2,0,0]],[[1,0,0],[3,0,0]],
                 [[0,1,0],[2,1,0]],[[1,1,0],[3,1,0]],
                 [[0,0,1]],[[2,0,1]],[[1,0,1],[3,0,1]],
                 [[0,1,1],[2,1,1]],[[1,1,1],[3,1,1]]];

AGClasses[(16,10)]=[[[0,0,0]],[[0,1,0]],[[2,0,0]],[[2,1,0]],[[0,0,1],[0,1,1]],[[2,0,1],[2,1,1]],
                  [[1,0,0],[1,1,0]],[[3,0,0],[3,1,0]],[[1,0,1],[1,1,1]],[[3,0,1],[3,1,1]]];


AGClasses[(16,11)]=[[[0,0,0]],[[2,0,0]],[[1,0,0],[3,0,0]],[[0,1,0],[2,1,0]],[[1,1,0],[3,1,0]],
                  [[0,0,1]],[[2,0,1]],[[1,0,1],[3,0,1]],[[0,1,1],[2,1,1]],[[1,1,1],[3,1,1]]];

AGClasses[(16,12)]=[[[0,0]],[[4,0]],[[2,0],[6,0]],[[1,0],[7,0]],[[3,0],[5,0]],
                  [[0,1],[2,1],[4,1],[6,1]],[[1,1],[3,1],[5,1],[7,1]]];

AGClasses[(16,13)]=[[[0,0]],[[4,0]],[[2,0],[6,0]],[[1,0],[3,0]],[[5,0],[7,0]],
                  [[0,1],[2,1],[4,1],[6,1]],[[1,1],[3,1],[5,1],[7,1]]];

AGClasses[(16,14)]=[[[0,0]],[[4,0]],[[2,0],[6,0]],[[1,0],[7,0]],[[3,0],[5,0]],
                  [[0,1],[2,1],[4,1],[6,1]],[[1,1],[3,1],[5,1],[7,1]]];

AGClasses[(24,1)]=[[[0,0,0]],[[2,0,0]],[[0,1,0],[0,2,0]],[[2,1,0],[2,2,0]],
                 [[0,0,1],[2,0,1]],[[0,2,1],[2,1,1]],[[0,1,1],[2,2,1]],
                 [[1,0,0],[1,1,0],[1,2,0],[3,0,0],[3,1,0],[3,2,0]],
                 [[1,0,1],[1,1,1],[1,2,1],[3,0,1],[3,1,1],[3,2,1]]];

AGClasses[(24,2)]=[[[0,0]],[[6,0]],[[4,0],[8,0]],[[2,0],[10,0]],[[3,0],[9,0]],
                  [[1,0],[11,0]],[[5,0],[7,0]],[[0,1],[2,1],[4,1],[6,1],[8,1],[10,1]],
                  [[1,1],[3,1],[5,1],[7,1],[9,1],[11,1]]];

AGClasses[(24,3)]=[[[0,0,0]],[[3,0,0]],[[1,0,0],[5,0,0]],[[2,0,0],[4,0,0]],
                 [[0,1,0],[2,1,0],[4,1,0]],[[1,1,0],[3,1,0],[5,1,0]],[[0,0,1]],
                 [[3,0,1]],[[1,0,1],[5,0,1]],[[2,0,1],[4,0,1]],[[0,1,1],[2,1,1],
                 [4,1,1]],[[1,1,1],[3,1,1],[5,1,1]]];


AGClasses[(24,4)]=[[[0,0,0]],[[1,0,0],[2,0,0]],[[0,1,0],[1,1,0],[2,1,0]],[[0,0,1]],
                 [[1,0,1],[2,0,1]],[[0,1,1],[1,1,1],[2,1,1]],[[0,0,2]],[[1,0,2],[2,0,2]],
                 [[0,1,2],[1,1,2],[2,1,2]],[[0,0,3]],[[1,0,3],[2,0,3]],[[0,1,3],[1,1,3],[2,1,3]]];

AGClasses[(24,5)]=[[[0,0,0,0]],[[1,0,0,0],[2,0,0,0]],[[0,1,0,0],[1,1,0,0],[2,1,0,0]],
                 [[0,0,1,0]],[[1,0,1,0],[2,0,1,0]],[[0,1,1,0],[1,1,1,0],[2,1,1,0]],
                 [[0,0,0,1]],[[1,0,0,1],[2,0,0,1]],[[0,1,0,1],[1,1,0,1],[2,1,0,1]],
                 [[0,0,1,1]],[[1,0,1,1],[2,0,1,1]],[[0,1,1,1],[1,1,1,1],[2,1,1,1]]];


AGClasses[(24,6)]=[[[0,0]],[[1,0],[7,0]],[[2,0]],[[3,0],[9,0]],[[4,0]],[[5,0],[11,0]],
                 [[6,0]],[[8,0]],[[10,0]],[[0,1],[6,1]],[[1,1],[7,1]],[[2,1],[8,1]],
                 [[3,1],[9,1]],[[4,1],[10,1]],[[5,1],[11,1]]];

AGClasses[(24,7)]=[[[0,0,0,0]],[[0,0,1,0],[0,1,0,0],[0,1,1,0]],
                 [[1,0,0,1],[0,0,1,1],[2,0,0,1],[1,1,1,1],[0,1,1,1],[2,0,1,1]],
                 [[0,0,0,1],[0,1,0,1],[1,1,0,1],[1,0,1,1],[2,1,1,1],[2,1,0,1]],
                 [[1,0,0,0],[2,0,0,0],[1,1,0,0],[1,0,1,0],[1,1,1,0],[2,1,1,0],[2,0,1,0],[2,1,0,0]]];


AGClasses[(24,8)]=[[[0,0,0]],[[3,0,0]],[[0,1,0],[0,0,1],[3,1,1]],[[3,1,0],[3,0,1],[0,1,1]],
                 [[1,1,0],[1,0,1],[4,0,0],[4,1,1]],[[4,1,0],[4,0,1],[1,0,0],[1,1,1]],
                 [[2,0,0],[2,1,1],[5,1,0],[5,0,1]],[[5,0,0],[5,1,1],[2,1,0],[2,0,1]]];

AGClasses[(24,9)]=[[[0,0,0]],[[3,0,0]],[[0,1,0],[0,0,1],[3,1,0],[3,0,1],[0,1,1],[3,1,1]],
                 [[1,0,0],[1,1,0],[1,1,1],[1,0,1]],[[4,0,0],[4,1,0],[4,1,1],[4,0,1]],
                 [[2,0,0],[5,0,1],[5,1,1],[5,1,0]],[[5,0,0],[2,0,1],[2,1,1],[2,1,0]]];


AGClasses[(24,10)]=[[[0,0,0]],[[0,1,0],[0,0,1],[0,1,1]],[[2,0,0],[2,1,0],[2,0,1],[2,1,1]],
                  [[4,0,0],[4,1,0],[4,0,1],[4,1,1]],[[3,0,0]],[[3,1,0],[3,0,1],[3,1,1]],
                  [[5,0,0],[5,1,0],[5,0,1],[5,1,1]],[[1,0,0],[1,1,0],[1,0,1],[1,1,1]]];


AGClasses[(24,11)]=[[[0,0]],[[6,0]],[[1,0],[11,0]],[[5,0],[7,0]],[[2,0],[10,0]],
                  [[4,0],[8,0]],[[3,0],[9,0]],[[0,1],[2,1],[4,1],[6,1],[8,1],[10,1]],
                  [[1,1],[3,1],[5,1],[7,1],[9,1],[11,1]]];


AGClasses[(24,12)]=[[[0,0]],[[1,0]],[[2,0]],[[3,0]],[[4,0]],[[5,0]],[[6,0]],[[7,0]],
                  [[8,0]],[[9,0]],[[10,0]],[[11,0]],[[0,1]],[[1,1]],[[2,1]],[[3,1]],
                  [[4,1]],[[5,1]],[[6,1]],[[7,1]],[[8,1]],[[9,1]],[[10,1]],[[11,1]]];


AGClasses[(32,1)]=[[[0,0,0]],[[0,2,0]],[[0,0,2]],[[0,2,2]],[[0,0,1],[0,0,3]],[[0,2,1],[0,2,3]],
                 [[1,0,1],[1,2,3]],[[1,2,1],[1,0,3]],[[1,0,0],[1,2,0]],[[1,0,2],[1,2,2]],
                 [[1,1,0],[1,3,0],[1,1,2],[1,3,2]],[[1,1,1],[1,1,3],[1,3,1],[1,3,3]],
                 [[0,1,0],[0,3,0],[0,1,2],[0,3,2]],[[0,3,1],[0,1,1],[0,3,3],[0,1,3]]];


AGClasses[(32,2)]=[[[0,0,0]],[[0,2,0]],[[2,0,0]],[[2,2,0]],[[3,0,1],[1,0,1]],[[3,2,1],[1,2,1]],
                 [[0,1,1],[0,3,1]],[[2,3,1],[2,1,1]],[[1,1,0],[3,3,0]],[[1,3,0],[3,1,0]],
                 [[0,0,1],[2,0,1],[0,2,1],[2,2,1]],[[1,0,0],[3,0,0],[1,2,0],[3,2,0]],
                 [[0,1,0],[2,1,0],[0,3,0],[2,3,0]],[[1,1,1],[3,1,1],[1,3,1],[3,3,1]]];


AGClasses[(32,3)]=[[[0,0,0,0]],[[2,0,0,0]],[[1,0,0,0],[3,0,0,0]],[[0,1,0,0],[2,1,0,0]],
                 [[1,1,0,0],[3,1,0,0]],[[0,0,1,0]],[[2,0,1,0]],[[1,0,1,0],[3,0,1,0]],
                 [[0,1,1,0],[2,1,1,0]],[[1,1,1,0],[3,1,1,0]],[[0,0,0,1]],[[2,0,0,1]],
                 [[1,0,0,1],[3,0,0,1]],[[0,1,0,1],[2,1,0,1]],[[1,1,0,1],[3,1,0,1]],
                 [[0,0,1,1]],[[2,0,1,1]],[[1,0,1,1],[3,0,1,1]],[[0,1,1,1],[2,1,1,1]],
                 [[1,1,1,1],[3,1,1,1]]];


AGClasses[(32,4)]=[[[0,0,0]],[[0,1,0]],[[0,2,0]],[[0,3,0]],[[2,0,0],[2,2,0]],[[2,1,0],[2,3,0]],
                 [[1,0,0],[3,3,0]],[[1,1,0],[3,0,0]],[[1,2,0],[3,1,0]],[[1,3,0],[3,2,0]],
                 [[1,1,1],[3,0,1],[1,3,1],[3,2,1]],[[1,0,1],[1,2,1],[3,1,1],[3,3,1]],
                 [[0,1,1],[0,3,1],[2,0,1],[2,2,1]],[[0,0,1],[0,2,1],[2,1,1],[2,3,1]]];


AGClasses[(32,5)]=[[[0,0,0,0]],[[0,1,0,0]],[[2,0,0,0]],[[2,1,0,0]],[[0,0,1,0],[0,1,1,0]],
                 [[2,0,1,0],[2,1,1,0]],[[1,0,0,0],[1,1,0,0]],[[3,0,0,0],[3,1,0,0]],
                 [[1,0,1,0],[1,1,1,0]],[[3,0,1,0],[3,1,1,0]],[[0,0,0,1]],[[0,1,0,1]],
                 [[2,0,0,1]],[[2,1,0,1]],[[0,0,1,1],[0,1,1,1]],[[2,0,1,1],[2,1,1,1]],
                 [[1,0,0,1],[1,1,0,1]],[[3,0,0,1],[3,1,0,1]],[[1,0,1,1],[1,1,1,1]],
                 [[3,0,1,1],[3,1,1,1]]];


AGClasses[(32,6)]=[[[0,0,0]],[[2,0,0]],[[1,0,0],[3,2,0]],[[3,0,0],[1,2,0]],[[0,2,0]],[[2,2,0]],
                  [[0,0,1],[2,0,1]],[[1,0,1],[1,2,1]],[[3,0,1],[3,2,1]],[[0,2,1],[2,2,1]],
                  [[1,3,0],[3,3,0],[1,1,0],[3,1,0]],[[0,1,0],[0,3,0],[2,1,0],[2,3,0]],
                  [[1,3,1],[3,3,1],[1,1,1],[3,1,1]],[[0,1,1],[0,3,1],[2,1,1],[2,3,1]]];


AGClasses[(32,7)]=[[[0,0,0]],[[2,0,0]],[[0,2,0]],[[2,2,0]],[[1,0,0],[3,0,0],[1,2,0],[3,2,0]],
                 [[0,1,0],[2,1,0]],[[0,3,0],[2,3,0]],[[1,1,0],[3,1,0],[1,3,0],[3,3,0]],
                 [[0,0,1],[2,2,1]],[[2,0,1],[0,2,1]],[[1,0,1],[3,0,1],[1,2,1],[3,2,1]],
                 [[0,1,1],[0,3,1]],[[2,1,1],[2,3,1]],[[1,1,1],[3,1,1],[1,3,1],[3,3,1]]];


AGClasses[(32,8)]=[[[0,0]],[[6,2]],[[4,0]],[[2,2]],[[0,1],[4,3]],[[6,3],[2,1]],
                 [[1,0],[5,2]],[[7,2],[3,0]],[[1,1],[5,3]],[[7,3],[3,1]],
                 [[0,2]],[[6,0]],[[4,2]],[[2,0]],[[0,3],[4,1]],[[6,1],[2,3]],
                 [[1,2],[5,0]],[[7,0],[3,2]],[[1,3],[5,1]],[[7,1],[3,3]]];

AGClasses[(32,9)]=[[[0,0,0]],[[4,0,0]],[[2,0,0],[6,0,0]],[[1,0,0],[7,0,0]],[[3,0,0],[5,0,0]],
                 [[0,1,0],[2,1,0],[4,1,0],[6,1,0]],[[1,1,0],[3,1,0],[5,1,0],[7,1,0]],
                 [[0,0,1]],[[4,0,1]],[[2,0,1],[6,0,1]],[[1,0,1],[7,0,1]],[[3,0,1],[5,0,1]],
                 [[0,1,1],[2,1,1],[4,1,1],[6,1,1]],[[1,1,1],[3,1,1],[5,1,1],[7,1,1]]];


AGClasses[(32,10)]=[[[0,0]],[[4,0]],[[2,0],[6,0]],[[1,0],[7,0]],[[3,0],[5,0]],
                  [[0,1],[2,1],[4,1],[6,1]],[[1,1],[3,1],[5,1],[7,1]],
                  [[0,2]],[[4,2]],[[2,2],[6,2]],[[1,2],[7,2]],[[3,2],[5,2]],
                  [[0,3],[2,3],[4,3],[6,3]],[[1,3],[3,3],[5,3],[7,3]]];
AGClasses[(32,11)]=[[[0,0]],[[4,0]],[[2,0],[6,0]],[[0,2]],[[4,2]],[[2,2],[6,2]],
                  [[1,0],[3,2]],[[1,2],[3,0]],[[5,0],[7,2]],[[5,2],[7,0]],
                  [[1,1],[3,3],[5,1],[7,3]],[[1,3],[3,1],[5,3],[7,1]],
                  [[0,1],[2,3],[4,1],[6,3]],[[0,3],[2,1],[4,3],[6,1]]];

AGClasses[(32,12)]=[[[0,0,0]],[[0,1,0]],[[2,0,0]],[[2,1,0]],[[4,0,0]],[[4,1,0]],
                  [[6,0,0]],[[6,1,0]],[[0,0,1],[0,1,1]],[[2,0,1],[2,1,1]],
                  [[4,0,1],[4,1,1]],[[6,0,1],[6,1,1]],[[1,0,0],[1,1,0]],
                  [[3,0,0],[3,1,0]],[[1,0,1],[1,1,1]],[[3,0,1],[3,1,1]],
                  [[5,0,0],[5,1,0]],[[7,0,0],[7,1,0]],[[5,0,1],[5,1,1]],[[7,0,1],[7,1,1]]];


AGClasses[(32,13)]=[[[0,0,0]],[[1,0,0],[3,2,0]],[[2,0,0]],[[3,0,0],[1,2,0]],
                  [[0,1,0]],[[1,1,0],[3,3,0]],[[2,1,0]],[[3,1,0],[1,3,0]],
                  [[0,2,0]],[[2,2,0]],[[0,3,0]],[[2,3,0]],[[0,0,1],[2,2,1]],
                  [[1,0,1],[3,2,1]],[[2,0,1],[0,2,1]],[[3,0,1],[1,2,1]],
                  [[0,1,1],[2,3,1]],[[1,1,1],[3,3,1]],[[2,1,1],[0,3,1]],[[3,1,1],[1,3,1]]];

AGClasses[(32,14)]=[[[0,0,0]],[[2,0,0]],[[0,2,0]],[[2,2,0]],[[1,0,0],[3,0,0]],[[1,2,0],[3,2,0]],
                  [[0,1,0],[2,1,0]],[[1,1,0],[3,1,0]],[[0,3,0],[2,3,0]],[[1,3,0],[3,3,0]],
                  [[0,0,1]],[[2,0,1]],[[0,2,1]],[[2,2,1]],[[1,0,1],[3,0,1]],[[1,2,1],[3,2,1]],
                  [[0,1,1],[2,1,1]],[[1,1,1],[3,1,1]],[[0,3,1],[2,3,1]],[[1,3,1],[3,3,1]]];

AGClasses[(32,15)]=[[[0,0,0,0]],[[2,0,0,0]],[[1,0,0,0],[3,0,0,0]],[[0,1,0,0],[2,1,0,0]],
                  [[1,1,0,0],[3,1,0,0]],[[0,0,1,0]],[[2,0,1,0]],[[1,0,1,0],[3,0,1,0]],
                  [[0,1,1,0],[2,1,1,0]],[[1,1,1,0],[3,1,1,0]],[[0,0,0,1]],[[2,0,0,1]],
                  [[1,0,0,1],[3,0,0,1]],[[0,1,0,1],[2,1,0,1]],[[1,1,0,1],[3,1,0,1]],
                  [[0,0,1,1]],[[2,0,1,1]],[[1,0,1,1],[3,0,1,1]],[[0,1,1,1],[2,1,1,1]],
                  [[1,1,1,1],[3,1,1,1]]];

AGClasses[(32,16)]=[[[0,0,0]],[[2,0,0]],[[1,0,0],[3,0,0]],[[0,1,0],[2,1,0]],[[1,1,0],[3,1,0]],
                  [[0,0,1]],[[2,0,1]],[[1,0,1],[3,0,1]],[[0,1,1],[2,1,1]],[[1,1,1],[3,1,1]],
                  [[0,0,2]],[[2,0,2]],[[1,0,2],[3,0,2]],[[0,1,2],[2,1,2]],[[1,1,2],[3,1,2]],
                  [[0,0,3]],[[2,0,3]],[[1,0,3],[3,0,3]],[[0,1,3],[2,1,3]],[[1,1,3],[3,1,3]]];

AGClasses[(32,17)]=[[[0,0]],[[1,0]],[[2,0]],[[3,0]],[[4,0]],[[5,0]],[[6,0]],[[7,0]],
                  [[0,1]],[[1,1]],[[2,1]],[[3,1]],[[4,1]],[[5,1]],[[6,1]],[[7,1]],
                  [[0,2]],[[1,2]],[[2,2]],[[3,2]],[[4,2]],[[5,2]],[[6,2]],[[7,2]],
                  [[0,3]],[[1,3]],[[2,3]],[[3,3]],[[4,3]],[[5,3]],[[6,3]],[[7,3]]];

AGClasses[(48,1)]=[[[0,0,0]],[[6,0,0]],[[3,0,1],[9,0,1]],[[4,0,0],[8,0,0]],[[2,0,0],[10,0,0]],
                 [[1,0,1],[5,0,1],[7,0,1],[11,0,1]],[[1,1,0],[3,1,0],[5,1,0],[7,1,0],[9,1,0],[11,1,0]],
                 [[0,1,1],[4,1,1],[8,1,1]],[[2,1,1],[6,1,1],[10,1,1]],[[0,0,1],[6,0,1]],
                 [[3,0,0],[9,0,0]],[[2,0,1],[4,0,1],[8,0,1],[10,0,1]],
                 [[1,0,0],[5,0,0],[7,0,0],[11,0,0]],
                 [[0,1,0],[2,1,0],[4,1,0],[6,1,0],[8,1,0],[10,1,0]],
                 [[1,1,1],[3,1,1],[5,1,1],[7,1,1],[9,1,1],[11,1,1]]];

AGClasses[(48,2)]=[[[0,0,0,0]],[[2,0,0,0]],[[0,1,0,0],[0,2,0,0]],[[2,1,0,0],[2,2,0,0]],
                 [[0,0,1,0],[2,0,1,0]],[[0,2,1,0],[2,1,1,0]],[[0,1,1,0],[2,2,1,0]],
                 [[1,0,0,0],[1,1,0,0],[1,2,0,0],[3,0,0,0],[3,1,0,0],[3,2,0,0]],
                 [[1,0,1,0],[1,1,1,0],[1,2,1,0],[3,0,1,0],[3,1,1,0],[3,2,1,0]],[[0,0,0,1]],
                 [[2,0,0,1]],[[0,1,0,1],[0,2,0,1]],[[2,1,0,1],[2,2,0,1]],[[0,0,1,1],[2,0,1,1]],
                 [[0,2,1,1],[2,1,1,1]],[[0,1,1,1],[2,2,1,1]],
                 [[1,0,0,1],[1,1,0,1],[1,2,0,1],[3,0,0,1],[3,1,0,1],[3,2,0,1]],
                 [[1,0,1,1],[1,1,1,1],[1,2,1,1],[3,0,1,1],[3,1,1,1],[3,2,1,1]]];


AGClasses[(48,3)]=[[[0,0,0]],[[3,0,0]],[[6,0,0]],[[9,0,0]],
                 [[0,1,0],[0,0,1],[6,1,0],[6,0,1],[3,1,1],[9,1,1]],
                 [[0,1,1],[3,1,0],[3,0,1],[6,1,1],[9,1,0],[9,0,1]],[[4,0,0],[1,1,0],[7,0,1],[4,1,1]],
                 [[7,0,0],[4,1,0],[10,0,1],[7,1,1]],[[10,0,0],[7,1,0],[1,0,1],[10,1,1]],
                 [[1,0,0],[10,1,0],[4,0,1],[1,1,1]],[[8,0,0],[11,1,0],[5,0,1],[2,1,1]],
                 [[11,0,0],[2,1,0],[8,0,1],[5,1,1]],[[2,0,0],[5,1,0],[11,0,1],[8,1,1]],
                 [[5,0,0],[8,1,0],[2,0,1],[11,1,1]]];

AGClasses[(48,4)]=[[[0,0,0,0]],[[3,0,0,0]],
                 [[0,1,0,0],[0,0,1,0],[3,1,0,0],[3,0,1,0],[0,1,1,0],[3,1,1,0]],
                 [[1,0,0,0],[1,1,0,0],[1,1,1,0],[1,0,1,0]],[[4,0,0,0],[4,1,0,0],[4,1,1,0],[4,0,1,0]],
                 [[2,0,0,0],[5,0,1,0],[5,1,1,0],[5,1,0,0]],[[5,0,0,0],[2,0,1,0],[2,1,1,0],[2,1,0,0]],
                 [[0,0,0,1]],[[3,0,0,1]],
                 [[0,1,0,1],[0,0,1,1],[3,1,0,1],[3,0,1,1],[0,1,1,1],[3,1,1,1]],
                 [[1,0,0,1],[1,1,0,1],[1,1,1,1],[1,0,1,1]],[[4,0,0,1],[4,1,0,1],[4,1,1,1],[4,0,1,1]],
                 [[2,0,0,1],[5,0,1,1],[5,1,1,1],[5,1,0,1]],[[5,0,0,1],[2,0,1,1],[2,1,1,1],[2,1,0,1]]];

AGClasses[(48,5)]=[[[0,0,0,0]],[[3,0,0,0]],[[0,1,0,0],[0,0,1,0],[3,1,1,0]],
                 [[3,1,0,0],[3,0,1,0],[0,1,1,0]],[[1,1,0,0],[1,0,1,0],[4,0,0,0],[4,1,1,0]],
                 [[4,1,0,0],[4,0,1,0],[1,0,0,0],[1,1,1,0]],[[2,0,0,0],[2,1,1,0],[5,1,0,0],[5,0,1,0]],
                 [[5,0,0,0],[5,1,1,0],[2,1,0,0],[2,0,1,0]],
                 [[0,0,0,1]],[[3,0,0,1]],[[0,1,0,1],[0,0,1,1],[3,1,1,1]],
                 [[3,1,0,1],[3,0,1,1],[0,1,1,1]],[[1,1,0,1],[1,0,1,1],[4,0,0,1],[4,1,1,1]],
                 [[4,1,0,1],[4,0,1,1],[1,0,0,1],[1,1,1,1]],[[2,0,0,1],[2,1,1,1],[5,1,0,1],[5,0,1,1]],
                 [[5,0,0,1],[5,1,1,1],[2,1,0,1],[2,0,1,1]]];

AGClasses[(48,6)]=[[[0,0,0,0]],[[2,0,0,0]],
                 [[1,0,0,0],[0,1,0,0],[1,1,0,0],[3,0,0,0],[2,1,0,0],[3,1,0,0]],
                 [[0,0,1,0],[0,0,2,0],[3,0,1,0],[2,1,1,0],[1,1,1,0],[0,1,2,0],[1,0,2,0],[3,1,2,0]],
                 [[2,0,1,0],[2,0,2,0],[1,0,1,0],[0,1,1,0],[3,1,1,0],[2,1,2,0],[3,0,2,0],[1,1,2,0]],
                 [[0,1,0,1],[3,0,0,1],[2,1,1,1],[3,1,1,1],[1,0,2,1],[1,1,2,1]],
                 [[2,1,0,1],[1,0,0,1],[0,1,1,1],[1,1,1,1],[3,0,2,1],[3,1,2,1]],
                 [[1,1,0,1],[3,1,0,1],[0,0,1,1],[2,0,1,1],[0,0,2,1],[2,0,2,1],
                  [0,0,0,1],[2,0,0,1],[3,0,1,1],[1,0,1,1],[0,1,2,1],[2,1,2,1]]];

AGClasses[(48,7)]=[[[0,0,0,0]],[[3,0,0,0]],[[0,1,0,0],[0,0,1,0],[3,1,1,0]],
                 [[3,1,0,0],[3,0,1,0],[0,1,1,0]],
                 [[1,0,1,0],[4,0,0,0],[1,1,0,0],[4,1,1,0],[5,1,0,0],[2,0,0,0],[2,1,1,0],[5,0,1,0]],
                 [[4,0,1,0],[1,0,0,0],[4,1,0,0],[1,1,1,0],[2,1,0,0],[5,0,0,0],[5,1,1,0],[2,0,1,0]],
                 [[4,0,0,1],[5,1,0,1],[3,1,0,1],[1,0,1,1],[2,0,0,1],[3,0,1,1]],
                 [[1,0,0,1],[2,1,0,1],[0,1,0,1],[4,0,1,1],[5,0,0,1],[0,0,1,1]],
                 [[0,0,0,1],[0,1,1,1],[2,1,1,1],[4,1,1,1],[5,0,1,1],[1,1,0,1]],
                 [[3,0,0,1],[3,1,1,1],[5,1,1,1],[1,1,1,1],[2,0,1,1],[4,1,0,1]]];

AGClasses[(48,8)]=[[[0,0,0]],[[0,2,0]],
                 [[3,1,0],[1,3,0],[0,0,1],[2,2,1],[1,3,1],[2,2,2],[0,0,2],[3,1,1]],
                 [[3,3,0],[1,1,0],[0,2,1],[2,0,1],[1,1,1],[2,0,2],[0,2,2],[3,3,1]],
                 [[2,2,0],[3,1,2],[1,3,2]],[[2,0,0],[3,3,2],[1,1,2]],
                 [[0,1,0],[1,0,1],[0,1,1],[2,3,2],[3,2,2],[0,1,2]],
                 [[0,3,0],[1,2,1],[0,3,1],[2,1,2],[3,0,2],[0,3,2]],
                 [[1,0,0],[1,0,2],[2,3,0],[3,2,0],[2,3,1],[3,2,1]],
                 [[1,2,0],[1,2,2],[2,1,0],[3,0,0],[2,1,1],[3,0,1]]];


AGClasses[(48,9)]=[[[0,0,0]],[[3,0,0]],[[2,0,0],[4,0,0]],[[1,0,0],[5,0,0]],[[0,2,0]],[[3,2,0]],
                 [[2,2,0],[4,2,0]],[[1,2,0],[5,2,0]],[[0,0,1],[3,0,1]],[[0,2,1],[3,2,1]],
                 [[2,0,1],[1,0,1]],[[2,2,1],[1,2,1]],[[5,0,1],[4,0,1]],[[5,2,1],[4,2,1]],
                 [[0,1,0],[1,1,0],[2,1,0],[3,1,0],[4,1,0],[5,1,0]],
                 [[0,3,0],[1,3,0],[2,3,0],[3,3,0],[4,3,0],[5,3,0]],
                 [[0,1,1],[1,1,1],[2,1,1],[3,1,1],[4,1,1],[5,1,1]],
                 [[0,3,1],[1,3,1],[2,3,1],[3,3,1],[4,3,1],[5,3,1]]];


AGClasses[(48,10)]=[[[0,0,0]],[[4,0,0]],[[2,0,0],[5,1,1],[7,1,1],[6,0,0],[1,1,1],[3,1,1]],
                  [[4,1,0],[6,1,0],[5,0,1],[7,0,1],[4,2,0],[1,2,1],[3,2,1],[2,2,0]],
                  [[0,1,0],[2,1,0],[1,0,1],[3,0,1],[0,2,0],[5,2,1],[7,2,1],[6,2,0]],
                  [[1,0,0],[3,2,0],[6,0,1],[7,0,0],[2,2,1],[5,1,0]],
                  [[5,0,0],[7,2,0],[2,0,1],[3,0,0],[6,2,1],[1,1,0]],
                  [[3,1,0],[0,0,1],[1,2,0],[6,1,1],[4,2,1],[4,1,1],
                   [7,1,0],[4,0,1],[5,2,0],[2,1,1],[0,2,1],[0,1,1]]];

AGClasses[(48,11)]=[[[0,0,0]],[[1,0,0],[7,0,1]],[[2,0,0]],[[3,0,0],[9,0,1]],
                  [[4,0,0]],[[5,0,0],[11,0,1]],[[6,0,0]],[[8,0,0]],[[10,0,0]],
                  [[0,1,0],[6,1,1]],[[1,1,0],[7,1,1]],[[2,1,0],[8,1,1]],
                  [[3,1,0],[9,1,1]],[[4,1,0],[10,1,1]],[[5,1,0],[11,1,1]],
                  [[0,0,1]],[[1,0,1],[7,0,0]],[[2,0,1]],[[3,0,1],[9,0,0]],
                  [[4,0,1]],[[5,0,1],[11,0,0]],[[6,0,1]],[[8,0,1]],[[10,0,1]],
                  [[0,1,1],[6,1,0]],[[1,1,1],[7,1,0]],[[2,1,1],[8,1,0]],
                  [[3,1,1],[9,1,0]],[[4,1,1],[10,1,0]],[[5,1,1],[11,1,0]]];

AGClasses[(48,12)]=[[[0,0]],[[6,0]],[[4,0],[8,0]],[[2,0],[10,0]],[[3,0],[9,0]],
                  [[1,0],[11,0]],[[5,0],[7,0]],[[0,1],[2,1],[4,1],[6,1],[8,1],[10,1]],
                  [[1,1],[3,1],[5,1],[7,1],[9,1],[11,1]],[[0,2]],[[6,2]],[[4,2],[8,2]],
                  [[2,2],[10,2]],[[3,2],[9,2]],[[1,2],[11,2]],[[5,2],[7,2]],
                  [[0,3],[2,3],[4,3],[6,3],[8,3],[10,3]],[[1,3],[3,3],[5,3],[7,3],[9,3],[11,3]]];


AGClasses[(48,13)]=[[[0,0,0]],[[3,0,0]],[[0,0,2]],[[3,0,2]],[[1,0,0],[5,0,0]],[[2,0,0],[4,0,0]],
                  [[1,0,2],[5,0,2]],[[2,0,2],[4,0,2]],[[0,0,1],[3,0,1]],[[0,0,3],[3,0,3]],
                  [[1,0,1],[2,0,1]],[[4,0,1],[5,0,1]],[[1,0,3],[2,0,3]],[[4,0,3],[5,0,3]],
                  [[0,1,0],[1,1,0],[2,1,0],[3,1,0],[4,1,0],[5,1,0]],
                  [[0,1,2],[1,1,2],[2,1,2],[3,1,2],[4,1,2],[5,1,2]],
                  [[0,1,1],[1,1,1],[2,1,1],[3,1,1],[4,1,1],[5,1,1]],
                  [[0,1,3],[1,1,3],[2,1,3],[3,1,3],[4,1,3],[5,1,3]]];


AGClasses[(48,14)]=[[[0,0,0]],[[2,0,0]],[[0,1,0],[0,5,0]],[[2,1,0],[2,5,0]],[[0,0,1],[2,0,1]],
                  [[0,2,1],[2,4,1]],[[0,1,1],[2,5,1]],
                  [[1,0,0],[1,4,0],[1,2,0],[3,0,0],[3,4,0],[3,2,0]],
                  [[1,0,1],[1,4,1],[1,2,1],[3,0,1],[3,4,1],[3,2,1]],
                  [[0,3,0]],[[2,3,0]],[[0,4,0],[0,2,0]],[[2,4,0],[2,2,0]],
                  [[0,3,1],[2,3,1]],[[0,5,1],[2,1,1]],[[0,4,1],[2,2,1]],
                  [[1,3,0],[1,1,0],[1,5,0],[3,3,0],[3,1,0],[3,5,0]],
                  [[1,3,1],[1,1,1],[1,5,1],[3,3,1],[3,1,1],[3,5,1]]];

AGClasses[(48,15)]=[[[0,0,0]],[[6,0,0]],[[1,0,0],[11,0,0]],[[5,0,0],[7,0,0]],[[2,0,0],[10,0,0]],
                  [[4,0,0],[8,0,0]],[[3,0,0],[9,0,0]],
                  [[0,1,0],[2,1,0],[4,1,0],[6,1,0],[8,1,0],[10,1,0]],
                  [[1,1,0],[3,1,0],[5,1,0],[7,1,0],[9,1,0],[11,1,0]],
                  [[0,0,1]],[[6,0,1]],[[1,0,1],[11,0,1]],[[5,0,1],[7,0,1]],[[2,0,1],[10,0,1]],
                  [[4,0,1],[8,0,1]],[[3,0,1],[9,0,1]],
                  [[0,1,1],[2,1,1],[4,1,1],[6,1,1],[8,1,1],[10,1,1]],
                  [[1,1,1],[3,1,1],[5,1,1],[7,1,1],[9,1,1],[11,1,1]]];

AGClasses[(64,1)]=[[[0,0,0]],[[0,0,4]],[[0,2,0]],[[0,2,4]],[[0,0,2],[0,0,6]],[[0,2,2],[0,2,6]],
                 [[0,0,1],[0,0,7]],[[0,0,5],[0,0,3]],[[0,2,1],[0,2,7]],[[0,2,5],[0,2,3]],
                 [[1,0,1],[1,2,3]],[[1,0,5],[1,2,7]],[[1,2,1],[1,0,3]],[[1,2,5],[1,0,7]],
                 [[1,0,0],[1,2,4]],[[1,0,4],[1,2,0]],[[1,0,2],[1,2,2]],[[1,0,6],[1,2,6]],
                 [[1,1,0],[1,1,2],[1,3,4],[1,3,6],[1,1,4],[1,1,6],[1,3,0],[1,3,2]],
                 [[1,1,1],[1,1,3],[1,3,5],[1,3,7],[1,1,5],[1,1,7],[1,3,1],[1,3,3]],
                 [[0,1,0],[0,1,2],[0,3,4],[0,3,6],[0,1,4],[0,1,6],[0,3,0],[0,3,2]],
                 [[0,1,3],[0,1,1],[0,3,7],[0,3,5],[0,1,7],[0,1,5],[0,3,3],[0,3,1]]];


AGClasses[(64,2)]=[[[0,0,0]],[[4,0,0]],[[4,2,0]],[[0,2,0]],[[2,0,0],[6,0,0]],[[2,2,0],[6,2,0]],
                  [[3,0,1],[1,0,1],[7,0,1],[5,0,1]],[[3,2,1],[1,2,1],[7,2,1],[5,2,1]],
                  [[0,1,1],[4,3,1]],[[4,1,1],[0,3,1]],[[2,3,1],[6,1,1],[6,3,1],[2,1,1]],
                  [[1,1,0],[7,3,0],[5,1,0],[3,3,0]],[[1,3,0],[7,1,0],[5,3,0],[3,1,0]],
                  [[0,0,1],[2,0,1],[4,2,1],[6,2,1],[4,0,1],[6,0,1],[0,2,1],[2,2,1]],
                  [[1,0,0],[7,0,0],[5,2,0],[3,2,0]],[[5,0,0],[3,0,0],[1,2,0],[7,2,0]],
                  [[0,1,0],[2,1,0],[4,3,0],[6,3,0],[4,1,0],[6,1,0],[0,3,0],[2,3,0]],
                  [[1,1,1],[7,1,1],[5,3,1],[3,3,1]],[[5,1,1],[3,1,1],[1,3,1],[7,3,1]]];
AGClasses[(64,3)]=[[[0,0]],[[4,0]],[[4,2]],[[0,2]],[[0,4]],[[4,4]],[[4,6]],[[0,6]],
                 [[2,0],[6,4]],[[6,0],[2,4]],[[6,2],[2,6]],[[2,2],[6,6]],[[1,0],[3,6]],
                 [[5,0],[7,6]],[[5,2],[7,0]],[[1,2],[3,0]],[[1,4],[3,2]],[[5,4],[7,2]],
                 [[5,6],[7,4]],[[1,6],[3,4]],[[0,1],[4,5],[2,7],[6,3]],[[4,1],[0,5],[6,7],[2,3]],
                 [[4,3],[0,7],[6,1],[2,5]],[[0,3],[4,7],[2,1],[6,5]],[[5,3],[1,7],[7,1],[3,5]],
                 [[1,3],[5,7],[3,1],[7,5]],[[1,5],[5,1],[3,3],[7,7]],[[5,5],[1,1],[7,3],[3,7]]];

AGClasses[(64,4)]=[[[0,0,0]],[[4,0,0]],[[2,0,0],[6,0,0]],[[0,2,0]],[[4,2,0]],[[2,2,0],[6,2,0]],
                 [[1,0,0],[3,2,0]],[[1,2,0],[3,0,0]],[[5,0,0],[7,2,0]],[[5,2,0],[7,0,0]],
                 [[1,1,0],[3,3,0],[5,1,0],[7,3,0]],[[1,3,0],[3,1,0],[5,3,0],[7,1,0]],
                 [[0,1,0],[2,3,0],[4,1,0],[6,3,0]],[[0,3,0],[2,1,0],[4,3,0],[6,1,0]],
                 [[0,0,1]],[[4,0,1]],[[2,0,1],[6,0,1]],[[0,2,1]],[[4,2,1]],[[2,2,1],[6,2,1]],
                 [[1,0,1],[3,2,1]],[[1,2,1],[3,0,1]],[[5,0,1],[7,2,1]],[[5,2,1],[7,0,1]],
                 [[1,1,1],[3,3,1],[5,1,1],[7,3,1]],[[1,3,1],[3,1,1],[5,3,1],[7,1,1]],
                 [[0,1,1],[2,3,1],[4,1,1],[6,3,1]],[[0,3,1],[2,1,1],[4,3,1],[6,1,1]]];


AGClasses[(64,5)]=[[[0,0,0]],[[4,0,0]],[[4,2,0]],[[0,2,0]],[[7,1,1],[5,3,1],[3,1,1],[1,3,1]],
                 [[3,3,1],[1,1,1],[7,3,1],[5,1,1]],[[2,0,0],[6,0,0]],[[6,2,0],[2,2,0]],
                 [[4,3,1],[6,3,1],[0,1,1],[2,1,1],[0,3,1],[2,3,1],[4,1,1],[6,1,1]],
                 [[1,0,0],[7,0,0],[5,2,0],[3,2,0]],[[5,0,0],[3,0,0],[1,2,0],[7,2,0]],
                 [[0,0,1],[4,2,1]],[[4,0,1],[0,2,1]],[[7,1,0],[1,1,0],[3,1,0],[5,1,0]],
                 [[3,3,0],[5,3,0],[7,3,0],[1,3,0]],[[2,0,1],[6,0,1],[6,2,1],[2,2,1]],
                 [[0,1,0],[2,1,0],[4,3,0],[6,3,0],[4,1,0],[6,1,0],[0,3,0],[2,3,0]],
                 [[5,2,1],[3,2,1],[1,0,1],[7,0,1]],[[1,2,1],[7,2,1],[5,0,1],[3,0,1]]];


AGClasses[(96,1)]=[[[0,0,0,0]],[[3,0,0,0]],[[6,0,0,0]],[[9,0,0,0]],
                 [[0,1,0,0],[0,0,1,0],[9,1,1,0],[6,1,0,0],[6,0,1,0],[3,1,1,0]],
                 [[3,1,0,0],[9,1,0,0],[3,0,1,0],[9,0,1,0],[0,1,1,0],[6,1,1,0]],
                 [[7,0,1,0],[4,0,0,0],[1,1,0,0],[4,1,1,0],[8,0,0,0],[11,1,0,0],[2,1,1,0],[5,0,1,0]],
                 [[10,0,1,0],[7,0,0,0],[4,1,0,0],[7,1,1,0],[11,0,0,0],[2,1,0,0],[5,1,1,0],[8,0,1,0]],
                 [[1,0,1,0],[10,0,0,0],[7,1,0,0],[10,1,1,0],[2,0,0,0],[5,1,0,0],[8,1,1,0],[11,0,1,0]],
                 [[4,0,1,0],[1,0,0,0],[10,1,0,0],[1,1,1,0],[5,0,0,0],[8,1,0,0],[11,1,1,0],[2,0,1,0]],
                 [[0,0,0,1],[8,1,1,1],[10,1,1,1],[3,0,1,1],[2,0,0,1],[1,0,1,1]],
                 [[3,0,0,1],[11,1,1,1],[1,1,1,1],[6,0,1,1],[5,0,0,1],[4,0,1,1]],
                 [[6,0,0,1],[2,1,1,1],[4,1,1,1],[9,0,1,1],[8,0,0,1],[7,0,1,1]],
                 [[9,0,0,1],[5,1,1,1],[7,1,1,1],[0,0,1,1],[11,0,0,1],[10,0,1,1]],
                 [[4,0,0,1],[10,0,0,1],[7,1,0,1],[1,1,0,1],[11,1,0,1],[5,1,0,1],
                  [9,1,0,1],[3,1,0,1],[5,0,1,1],[11,0,1,1],[6,1,1,1],[0,1,1,1]],
                 [[7,0,0,1],[1,0,0,1],[10,1,0,1],[4,1,0,1],[2,1,0,1],[8,1,0,1],
                  [0,1,0,1],[6,1,0,1],[8,0,1,1],[2,0,1,1],[9,1,1,1],[3,1,1,1]]];


AGClasses[(96,2)]=[[[0,0,0,0,0]],[[3,0,0,0,0]],[[0,0,1,0,0],[0,1,0,0,0],[0,1,1,0,0]],
                 [[3,0,1,0,0],[3,1,0,0,0],[3,1,1,0,0]],[[0,0,0,0,1],[3,0,0,0,1]],
                 [[0,0,1,0,1],[0,1,0,0,1],[0,1,1,0,1],[3,0,1,0,1],[3,1,0,0,1],[3,1,1,0,1]],
                 [[2,0,0,0,0],[4,0,0,0,0],[2,1,0,0,0],[2,0,1,0,0],
                  [2,1,1,0,0],[4,1,1,0,0],[4,0,1,0,0],[4,1,0,0,0]],
                 [[1,0,0,0,0],[5,0,0,0,0],[5,1,0,0,0],[5,0,1,0,0],
                  [5,1,1,0,0],[1,1,1,0,0],[1,0,1,0,0],[1,1,0,0,0]],
                 [[1,0,0,1,0],[0,0,1,1,0],[2,0,0,1,0],[1,1,1,1,0],[0,1,1,1,0],[2,0,1,1,0],
                  [4,0,0,1,0],[3,0,1,1,0],[5,0,0,1,0],[4,1,1,1,0],[3,1,1,1,0],[5,0,1,1,0]],
                 [[0,0,0,1,0],[0,1,0,1,0],[1,1,0,1,0],[1,0,1,1,0],[2,1,1,1,0],[2,1,0,1,0],
                  [3,0,0,1,0],[3,1,0,1,0],[4,1,0,1,0],[4,0,1,1,0],[5,1,1,1,0],[5,1,0,1,0]],
                 [[5,0,0,0,1],[5,1,0,0,1],[5,0,1,0,1],[5,1,1,0,1],
                  [4,0,0,0,1],[4,1,0,0,1],[4,0,1,0,1],[4,1,1,0,1]],
                 [[1,0,0,0,1],[1,1,0,0,1],[1,0,1,0,1],[1,1,1,0,1],
                  [2,0,0,0,1],[2,1,0,0,1],[2,0,1,0,1],[2,1,1,0,1]],
                 [[1,0,0,1,1],[0,0,1,1,1],[2,0,0,1,1],[1,1,1,1,1],[0,1,1,1,1],[2,0,1,1,1],
                  [4,0,0,1,1],[3,0,1,1,1],[5,0,0,1,1],[4,1,1,1,1],[3,1,1,1,1],[5,0,1,1,1]],
                 [[0,0,0,1,1],[0,1,0,1,1],[1,1,0,1,1],[1,0,1,1,1],[2,1,1,1,1],[2,1,0,1,1],
                  [3,0,0,1,1],[3,1,0,1,1],[4,1,0,1,1],[4,0,1,1,1],[5,1,1,1,1],[5,1,0,1,1]]];


AGClasses[(96,3)]=[[[0,0,0,0,0]],[[3,0,0,0,0]],[[0,0,1,0,0],[0,1,0,0,0],[0,1,1,0,0]],
                 [[3,0,1,0,0],[3,1,0,0,0],[3,1,1,0,0]],[[0,0,0,0,1],[3,0,0,0,1]],
                 [[0,0,1,0,1],[0,1,0,0,1],[0,1,1,0,1],[3,0,1,0,1],[3,1,0,0,1],[3,1,1,0,1]],
                 [[2,0,0,0,0],[4,0,0,0,0],[2,1,0,0,0],[2,0,1,0,0],
                  [2,1,1,0,0],[4,1,1,0,0],[4,0,1,0,0],[4,1,0,0,0]],
                 [[1,0,0,0,0],[5,0,0,0,0],[5,1,0,0,0],[5,0,1,0,0],
                  [5,1,1,0,0],[1,1,1,0,0],[1,0,1,0,0],[1,1,0,0,0]],
                 [[1,0,0,1,0],[0,0,1,1,0],[2,0,0,1,0],[1,1,1,1,0],[0,1,1,1,0],[2,0,1,1,0],
                  [4,0,0,1,0],[3,0,1,1,0],[5,0,0,1,0],[4,1,1,1,0],[3,1,1,1,0],[5,0,1,1,0]],
                 [[0,0,0,1,0],[0,1,0,1,0],[1,1,0,1,0],[1,0,1,1,0],[2,1,1,1,0],[2,1,0,1,0],
                  [3,0,0,1,0],[3,1,0,1,0],[4,1,0,1,0],[4,0,1,1,0],[5,1,1,1,0],[5,1,0,1,0]],
                 [[5,0,0,0,1],[5,1,0,0,1],[5,0,1,0,1],[5,1,1,0,1],
                  [4,0,0,0,1],[4,1,0,0,1],[4,0,1,0,1],[4,1,1,0,1]],
                 [[1,0,0,0,1],[1,1,0,0,1],[1,0,1,0,1],[1,1,1,0,1],
                  [2,0,0,0,1],[2,1,0,0,1],[2,0,1,0,1],[2,1,1,0,1]],
                 [[1,0,0,1,1],[0,0,1,1,1],[2,0,0,1,1],[1,1,1,1,1],[0,1,1,1,1],[2,0,1,1,1],
                  [4,0,0,1,1],[3,0,1,1,1],[5,0,0,1,1],[4,1,1,1,1],[3,1,1,1,1],[5,0,1,1,1]],
                 [[0,0,0,1,1],[0,1,0,1,1],[1,1,0,1,1],[1,0,1,1,1],[2,1,1,1,1],[2,1,0,1,1],
                  [3,0,0,1,1],[3,1,0,1,1],[4,1,0,1,1],[4,0,1,1,1],[5,1,1,1,1],[5,1,0,1,1]]];


AGClasses[(96,4)]=[[[0,0,0,0,0]],[[3,0,0,0,0]],[[3,1,0,0,0],[3,0,1,0,0],[0,1,1,0,0]],
                  [[0,1,0,0,0],[0,0,1,0,0],[3,1,1,0,0]],[[0,0,0,0,1],[3,0,0,0,1]],
                  [[0,1,0,0,1],[0,0,1,0,1],[3,1,1,0,1],[3,1,0,0,1],[3,0,1,0,1],[0,1,1,0,1]],
                  [[2,0,0,0,0],[4,0,0,0,0],[5,1,0,0,0],[5,0,1,0,0],
                   [2,1,1,0,0],[4,1,1,0,0],[1,1,0,0,0],[1,0,1,0,0]],
                  [[1,0,0,0,0],[5,0,0,0,0],[2,1,0,0,0],[2,0,1,0,0],
                   [5,1,1,0,0],[1,1,1,0,0],[4,1,0,0,0],[4,0,1,0,0]],
                  [[1,0,0,1,0],[0,0,1,1,0],[2,0,0,1,0],[4,0,0,1,0],[3,0,1,1,0],[5,0,0,1,0],
                   [0,1,0,1,0],[3,1,0,1,0],[1,0,1,1,0],[4,0,1,1,0],[2,1,0,1,0],[5,1,0,1,0]],
                  [[1,1,1,1,0],[0,1,1,1,0],[2,0,1,1,0],[4,1,1,1,0],[3,1,1,1,0],[5,0,1,1,0],
                   [0,0,0,1,0],[3,0,0,1,0],[1,1,0,1,0],[4,1,0,1,0],[2,1,1,1,0],[5,1,1,1,0]],
                  [[5,0,0,0,1],[4,0,0,0,1],[2,1,0,0,1],[2,0,1,0,1],
                   [5,1,1,0,1],[4,1,1,0,1],[1,1,0,0,1],[1,0,1,0,1]],
                  [[2,0,0,0,1],[1,0,0,0,1],[5,1,0,0,1],[5,0,1,0,1],
                   [2,1,1,0,1],[1,1,1,0,1],[4,1,0,0,1],[4,0,1,0,1]],
                  [[1,0,0,1,1],[0,0,1,1,1],[2,0,0,1,1],[4,0,0,1,1],[3,0,1,1,1],[5,0,0,1,1],
                   [0,1,0,1,1],[3,1,0,1,1],[1,0,1,1,1],[4,0,1,1,1],[2,1,0,1,1],[5,1,0,1,1]],
                  [[1,1,1,1,1],[0,1,1,1,1],[2,0,1,1,1],[4,1,1,1,1],[3,1,1,1,1],[5,0,1,1,1],
                   [0,0,0,1,1],[3,0,0,1,1],[4,1,0,1,1],[1,1,0,1,1],[2,1,1,1,1],[5,1,1,1,1]]];

AGClasses[(96,5)]=[[[0,0,0,0]],[[3,0,0,0]],
                 [[0,1,0,0],[0,0,1,0],[3,1,0,0],[3,0,1,0],[0,1,1,0],[3,1,1,0]],
                 [[1,0,0,0],[1,1,0,0],[1,1,1,0],[1,0,1,0]],[[4,0,0,0],[4,1,0,0],[4,1,1,0],[4,0,1,0]],
                 [[2,0,0,0],[5,0,1,0],[5,1,1,0],[5,1,0,0]],[[5,0,0,0],[2,0,1,0],[2,1,1,0],[2,1,0,0]],
                 [[0,0,0,1]],[[3,0,0,1]],
                 [[0,1,0,1],[0,0,1,1],[3,1,0,1],[3,0,1,1],[0,1,1,1],[3,1,1,1]],
                 [[1,0,0,1],[1,1,0,1],[1,1,1,1],[1,0,1,1]],[[4,0,0,1],[4,1,0,1],[4,1,1,1],[4,0,1,1]],
                 [[2,0,0,1],[5,0,1,1],[5,1,1,1],[5,1,0,1]],[[5,0,0,1],[2,0,1,1],[2,1,1,1],[2,1,0,1]],
                 [[0,0,0,2]],[[3,0,0,2]],
                 [[0,1,0,2],[0,0,1,2],[3,1,0,2],[3,0,1,2],[0,1,1,2],[3,1,1,2]],
                 [[1,0,0,2],[1,1,0,2],[1,1,1,2],[1,0,1,2]],[[4,0,0,2],[4,1,0,2],[4,1,1,2],[4,0,1,2]],
                 [[2,0,0,2],[5,0,1,2],[5,1,1,2],[5,1,0,2]],[[5,0,0,2],[2,0,1,2],[2,1,1,2],[2,1,0,2]],
                 [[0,0,0,3]],[[3,0,0,3]],
                 [[0,1,0,3],[0,0,1,3],[3,1,0,3],[3,0,1,3],[0,1,1,3],[3,1,1,3]],
                 [[1,0,0,3],[1,1,0,3],[1,1,1,3],[1,0,1,3]],[[4,0,0,3],[4,1,0,3],[4,1,1,3],[4,0,1,3]],
                 [[2,0,0,3],[5,0,1,3],[5,1,1,3],[5,1,0,3]],[[5,0,0,3],[2,0,1,3],[2,1,1,3],[2,1,0,3]]];

AGClasses[(96,6)]=[[(cs + [0]) for cs in css] for css in AGClasses[(48,4)]]+[[(cs + [1]) for cs in css] for css in AGClasses[(48,4)]]

AGClasses[(96,7)]=[[[0,0,0,0]],[[0,0,3,0]],[[2,0,3,0]],[[2,0,0,0]],
                 [[1,0,0,0],[3,0,0,0],[0,1,3,0],[2,1,3,0],[3,1,0,0],[1,1,0,0]],
                 [[1,0,3,0],[3,0,3,0],[0,1,0,0],[2,1,0,0],[3,1,3,0],[1,1,3,0]],
                 [[0,0,1,0],[3,0,4,0],[2,1,1,0],[1,1,4,0],[0,0,5,0],[0,1,5,0],[3,1,2,0],[1,0,2,0]],
                 [[0,0,4,0],[3,0,1,0],[2,1,4,0],[1,1,1,0],[0,0,2,0],[0,1,2,0],[3,1,5,0],[1,0,5,0]],
                 [[2,0,4,0],[1,0,1,0],[0,1,4,0],[3,1,1,0],[2,0,2,0],[2,1,2,0],[1,1,5,0],[3,0,5,0]],
                 [[2,0,1,0],[1,0,4,0],[0,1,1,0],[3,1,4,0],[2,0,5,0],[2,1,5,0],[1,1,2,0],[3,0,2,0]],
                 [[2,1,0,1],[1,0,3,1],[0,1,4,1],[1,1,1,1],[3,0,5,1],[3,1,5,1]],
                 [[2,1,3,1],[1,0,0,1],[0,1,1,1],[1,1,4,1],[3,0,2,1],[3,1,2,1]],
                 [[0,1,3,1],[3,0,0,1],[2,1,1,1],[3,1,4,1],[1,0,2,1],[1,1,2,1]],
                 [[0,1,0,1],[3,0,3,1],[2,1,4,1],[3,1,1,1],[1,0,5,1],[1,1,5,1]],
                 [[2,0,3,1],[3,1,0,1],[2,1,5,1],[2,0,5,1],[3,0,4,1],[2,0,1,1],
                  [0,0,3,1],[1,1,0,1],[0,1,5,1],[0,0,5,1],[1,0,4,1],[0,0,1,1]],
                 [[2,0,0,1],[3,1,3,1],[2,1,2,1],[2,0,2,1],[3,0,1,1],[2,0,4,1],
                  [0,0,0,1],[1,1,3,1],[0,1,2,1],[0,0,2,1],[1,0,1,1],[0,0,4,1]]];

AGClasses[(96,8)]=[[(cs + [0]) for cs in css] for css in AGClasses[(48,10)]]+[[(cs + [1]) for cs in css] for css in AGClasses[(48,10)]]


AGClasses[(96,9)]=[[[0,0,0,0]],[[0,3,0,0]],[[2,0,0,0]],[[2,3,0,0]],[[0,1,0,0],[0,5,0,0]],
                 [[0,4,0,0],[0,2,0,0]],[[2,1,0,0],[2,5,0,0]],[[2,4,0,0],[2,2,0,0]],
                 [[0,0,1,0],[2,0,1,0]],[[0,3,1,0],[2,3,1,0]],[[0,2,1,0],[2,4,1,0]],
                 [[0,5,1,0],[2,1,1,0]],[[0,1,1,0],[2,5,1,0]],[[0,4,1,0],[2,2,1,0]],
                 [[1,0,0,0],[1,1,0,0],[1,2,0,0],[3,0,0,0],[3,1,0,0],[3,2,0,0],
                  [1,3,0,0],[1,4,0,0],[1,5,0,0],[3,3,0,0],[3,4,0,0],[3,5,0,0]],
                 [[1,0,1,0],[1,1,1,0],[1,2,1,0],[3,0,1,0],[3,1,1,0],[3,2,1,0],
                  [1,3,1,0],[1,4,1,0],[1,5,1,0],[3,3,1,0],[3,4,1,0],[3,5,1,0]],
                 [[0,0,0,1],[0,3,0,1]],[[2,0,0,1],[2,3,0,1]],[[0,1,0,1],[0,2,0,1]],
                 [[0,4,0,1],[0,5,0,1]],[[2,1,0,1],[2,2,0,1]],[[2,4,0,1],[2,5,0,1]],
                 [[0,0,1,1],[2,3,1,1]],[[0,3,1,1],[2,0,1,1]],[[0,2,1,1],[2,1,1,1]],
                 [[0,5,1,1],[2,4,1,1]],[[0,1,1,1],[2,2,1,1]],[[0,4,1,1],[2,5,1,1]],
                 [[1,0,0,1],[1,1,0,1],[1,2,0,1],[3,0,0,1],[3,1,0,1],[3,2,0,1],
                  [1,3,0,1],[1,4,0,1],[1,5,0,1],[3,3,0,1],[3,4,0,1],[3,5,0,1]],
                 [[1,0,1,1],[1,1,1,1],[1,2,1,1],[3,0,1,1],[3,1,1,1],[3,2,1,1],
                  [1,3,1,1],[1,4,1,1],[1,5,1,1],[3,3,1,1],[3,4,1,1],[3,5,1,1]]];

AGClasses[(96,10)]=[[[0,0,0]],[[6,0,0]],[[0,2,0]],[[6,2,0]],[[3,0,1],[9,0,1],[3,2,1],[9,2,1]],[[4,0,0],[8,0,0]],[[4,2,0],[8,2,0]],[[2,0,0],[10,0,0]],[[2,2,0],[10,2,0]],[[1,0,1],[5,2,1],[7,2,1],[11,0,1]],[[1,2,1],[5,0,1],[7,0,1],[11,2,1]],[[1,1,0],[3,1,0],[5,1,0],[7,1,0],[9,1,0],[11,1,0],[1,3,0],[3,3,0],[5,3,0],[7,3,0],[9,3,0],[11,3,0]],[[0,1,1],[4,1,1],[8,1,1],[0,3,1],[4,3,1],[8,3,1]],[[2,1,1],[6,1,1],[10,1,1],[2,3,1],[6,3,1],[10,3,1]],[[0,0,1],[6,2,1]],[[0,2,1],[6,0,1]],[[3,0,0],[9,0,0],[3,2,0],[9,2,0]],[[2,0,1],[4,2,1],[8,2,1],[10,0,1]],[[2,2,1],[4,0,1],[8,0,1],[10,2,1]],[[1,0,0],[5,2,0],[7,2,0],[11,0,0]],[[1,2,0],[5,0,0],[7,0,0],[11,2,0]],[[0,1,0],[2,1,0],[4,1,0],[6,1,0],[8,1,0],[10,1,0]],[[0,3,0],[2,3,0],[4,3,0],[6,3,0],[8,3,0],[10,3,0]],[[1,1,1],[3,1,1],[5,1,1],[7,1,1],[9,1,1],[11,1,1],[1,3,1],[3,3,1],[5,3,1],[7,3,1],[9,3,1],[11,3,1]]];

AGClasses[(192,1)]=[[[0,0,0,0]],[[3,2,0,0]],[[6,0,0,0]],[[9,2,0,0]],
                  [[6,0,1,0],[0,3,0,0],[9,3,1,0],[0,2,1,0],[6,1,0,0],[3,1,1,0]],
                  [[9,2,1,0],[3,1,0,0],[0,1,1,0],[3,0,1,0],[9,3,0,0],[6,3,1,0]],
                  [[11,1,0,0],[8,2,0,0],[5,2,1,0],[2,3,1,0],[7,0,1,0],[4,2,0,0],[4,1,1,0],[1,3,0,0]],
                  [[2,3,0,0],[11,0,0,0],[8,0,1,0],[5,1,1,0],[10,2,1,0],[7,0,0,0],[7,3,1,0],[4,1,0,0]],
                  [[5,1,0,0],[2,2,0,0],[11,2,1,0],[8,3,1,0],[1,0,1,0],[10,2,0,0],[10,1,1,0],[7,3,0,0]],
                  [[8,3,0,0],[5,0,0,0],[2,0,1,0],[11,1,1,0],[4,2,1,0],[1,0,0,0],[1,3,1,0],[10,1,0,0]],
                  [[0,0,0,1],[8,3,1,1],[10,1,1,1],[3,2,1,1],[2,2,0,1],[1,0,1,1]],
                  [[3,2,0,1],[11,1,1,1],[1,3,1,1],[6,0,1,1],[5,0,0,1],[4,2,1,1]],
                  [[6,0,0,1],[2,3,1,1],[4,1,1,1],[9,2,1,1],[8,2,0,1],[7,0,1,1]],
                  [[9,2,0,1],[5,1,1,1],[7,3,1,1],[0,0,1,1],[11,0,0,1],[10,2,1,1]],
                  [[4,2,0,1],[7,1,0,1],[11,1,0,1],[9,3,0,1],[5,2,1,1],[6,3,1,1],
                   [10,0,0,1],[1,3,0,1],[5,3,0,1],[3,1,0,1],[11,0,1,1],[0,1,1,1]],
                  [[7,0,0,1],[10,3,0,1],[2,3,0,1],[0,1,0,1],[8,0,1,1],[9,1,1,1],
                   [1,2,0,1],[4,1,0,1],[8,1,0,1],[6,3,0,1],[2,2,1,1],[3,3,1,1]],
                  [[0,2,0,0]],[[3,0,0,0]],[[6,2,0,0]],[[9,0,0,0]],
                  [[6,2,1,0],[0,1,0,0],[9,1,1,0],[0,0,1,0],[6,3,0,0],[3,3,1,0]],
                  [[9,0,1,0],[3,3,0,0],[0,3,1,0],[3,2,1,0],[9,1,0,0],[6,1,1,0]],
                  [[11,3,0,0],[8,0,0,0],[5,0,1,0],[2,1,1,0],[7,2,1,0],[4,0,0,0],[4,3,1,0],[1,1,0,0]],
                  [[2,1,0,0],[11,2,0,0],[8,2,1,0],[5,3,1,0],[10,0,1,0],[7,2,0,0],[7,1,1,0],[4,3,0,0]],
                  [[5,3,0,0],[2,0,0,0],[11,0,1,0],[8,1,1,0],[1,2,1,0],[10,0,0,0],[10,3,1,0],[7,1,0,0]],
                  [[8,1,0,0],[5,2,0,0],[2,2,1,0],[11,3,1,0],[4,0,1,0],[1,2,0,0],[1,1,1,0],[10,3,0,0]],
                  [[0,2,0,1],[8,1,1,1],[10,3,1,1],[3,0,1,1],[2,0,0,1],[1,2,1,1]],
                  [[3,0,0,1],[11,3,1,1],[1,1,1,1],[6,2,1,1],[5,2,0,1],[4,0,1,1]],
                  [[6,2,0,1],[2,1,1,1],[4,3,1,1],[9,0,1,1],[8,0,0,1],[7,2,1,1]],
                  [[9,0,0,1],[5,3,1,1],[7,1,1,1],[0,2,1,1],[11,2,0,1],[10,0,1,1]],
                  [[4,0,0,1],[7,3,0,1],[11,3,0,1],[9,1,0,1],[5,0,1,1],[6,1,1,1],
                   [10,2,0,1],[1,1,0,1],[5,1,0,1],[3,3,0,1],[11,2,1,1],[0,3,1,1]],
                  [[7,2,0,1],[10,1,0,1],[2,1,0,1],[0,3,0,1],[8,2,1,1],[9,3,1,1],
                   [1,0,0,1],[4,3,0,1],[8,3,0,1],[6,1,0,1],[2,0,1,1],[3,1,1,1]]];

AGClasses[(192,2)]=[[[0,0,0,0]],[[4,0,0,0]],[[0,0,3,0]],[[4,0,3,0]],
                  [[0,1,0,0],[2,1,0,0],[2,0,0,0],[4,1,0,0],[6,1,0,0],[6,0,0,0]],
                  [[0,1,3,0],[2,1,3,0],[2,0,3,0],[4,1,3,0],[6,1,3,0],[6,0,3,0]],
                  [[0,0,0,1],[0,0,3,1]],[[4,0,0,1],[4,0,3,1]],
                  [[0,1,0,1],[2,1,0,1],[2,0,0,1],[4,1,0,1],[6,1,0,1],[6,0,0,1],[0,1,3,1],[2,1,3,1],
                   [2,0,3,1],[4,1,3,1],[6,1,3,1],[6,0,3,1]],
                  [[4,0,4,0],[4,1,4,0],[2,1,4,0],[2,0,4,0],[4,0,2,0],[6,1,2,0],[6,0,2,0],[0,1,2,0]],
                  [[0,0,4,0],[0,1,4,0],[6,1,4,0],[6,0,4,0],[0,0,2,0],[2,1,2,0],[2,0,2,0],[4,1,2,0]],
                  [[4,0,1,0],[4,1,1,0],[2,1,1,0],[2,0,1,0],[4,0,5,0],[6,1,5,0],[6,0,5,0],[0,1,5,0]],
                  [[0,0,1,0],[0,1,1,0],[6,1,1,0],[6,0,1,0],[0,0,5,0],[2,1,5,0],[2,0,5,0],[4,1,5,0]],
                  [[1,0,0,0],[7,0,0,0],[3,0,4,0],[3,1,4,0],[7,1,2,0],[5,0,2,0],
                   [1,0,3,0],[7,0,3,0],[3,0,1,0],[3,1,1,0],[7,1,5,0],[5,0,5,0]],
                  [[5,0,0,0],[3,0,0,0],[7,0,4,0],[7,1,4,0],[3,1,2,0],[1,0,2,0],
                   [5,0,3,0],[3,0,3,0],[7,0,1,0],[7,1,1,0],[3,1,5,0],[1,0,5,0]],
                  [[1,1,0,0],[3,1,0,0],[5,0,4,0],[5,1,4,0],[7,0,2,0],[5,1,2,0],[1,1,3,0],[3,1,3,0],
                   [5,0,1,0],[5,1,1,0],[7,0,5,0],[5,1,5,0],[5,1,0,0],[7,1,0,0],[1,0,4,0],[1,1,4,0],
                   [3,0,2,0],[1,1,2,0],[5,1,3,0],[7,1,3,0],[1,0,1,0],[1,1,1,0],[3,0,5,0],[1,1,5,0]],
                  [[4,0,4,1],[4,1,4,1],[2,1,4,1],[2,0,4,1],[4,0,5,1],[6,1,5,1],[6,0,5,1],[0,1,5,1]],
                  [[0,0,4,1],[0,1,4,1],[6,1,4,1],[6,0,4,1],[0,0,5,1],[2,1,5,1],[2,0,5,1],[4,1,5,1]],
                  [[4,0,1,1],[4,1,1,1],[2,1,1,1],[2,0,1,1],[4,0,2,1],[6,1,2,1],[6,0,2,1],[0,1,2,1]],
                  [[0,0,1,1],[0,1,1,1],[6,1,1,1],[6,0,1,1],[0,0,2,1],[2,1,2,1],[2,0,2,1],[4,1,2,1]],
                  [[1,0,3,1],[7,0,3,1],[3,0,1,1],[3,1,1,1],[5,0,5,1],[7,1,5,1],
                   [1,0,0,1],[7,0,0,1],[3,0,4,1],[3,1,4,1],[5,0,2,1],[7,1,2,1]],
                  [[5,0,3,1],[3,0,3,1],[7,0,1,1],[7,1,1,1],[1,0,5,1],[3,1,5,1],
                   [5,0,0,1],[3,0,0,1],[7,0,4,1],[7,1,4,1],[1,0,2,1],[3,1,2,1]],
                  [[1,1,3,1],[3,1,3,1],[5,0,1,1],[5,1,1,1],[7,0,5,1],[5,1,5,1],[1,1,0,1],[3,1,0,1],
                   [5,0,4,1],[5,1,4,1],[7,0,2,1],[5,1,2,1],[5,1,3,1],[7,1,3,1],[1,0,1,1],[1,1,1,1],
                   [3,0,5,1],[1,1,5,1],[5,1,0,1],[7,1,0,1],[1,0,4,1],[1,1,4,1],[3,0,2,1],[1,1,2,1]]];