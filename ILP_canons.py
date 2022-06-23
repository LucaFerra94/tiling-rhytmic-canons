
import numpy as np
import logging
from gurobipy import Model, quicksum, GRB



logging.basicConfig(
	filename='ILP_alltest.log',
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')



def MakeMatrix(A, n):
    P = np.zeros((n, n))

    for j in range(n):
        for i in A:
            P[j, (i + j) % n] = 1.0

    return np.transpose(P)


def ILP(A, n, D, num_test):

    # Size of complementar
    c = len(A)
    m = int(n // c)

    # Modello Gurobi    
    mod = Model()
    mod.setParam(GRB.Param.OutputFlag, 0)

    mod.setAttr(GRB.Attr.ModelSense, GRB.MINIMIZE)
    mod.setParam(GRB.Param.Method, 1)

    # Partitioning variables
    x = {}
    for i in range(n):
        x[i] = mod.addVar(vtype=GRB.BINARY, obj=1.0)

    U = {}
    y = {}
    z = {}
    for j, d in enumerate(D):
        for i in range(d):
            U[i,j] =  mod.addVar(vtype=GRB.BINARY)
            y[i,j] =  mod.addVar(vtype=GRB.BINARY)
            z[i,j] =  mod.addVar(vtype=GRB.BINARY)
    
    mod.update()

    # Constraints: Complementarity
    T = MakeMatrix(A, n)
    for row in T:
        mod.addConstr(
            quicksum(x[i] for i, a in enumerate(row) if a > 0.0) == 1.0)

    mod.addConstr(x[0] == 1.0)

    # w.l.o.g., in order to reduce translations it finished with c zeros"
    for j in range(n-c, n):
        mod.addConstr(x[j] == 0.0)
        
    # Aperiodicity constraints
    for j, d in enumerate(D):
         for i in range(d):


             mod.addConstr(quicksum(x[i + k*d] for k in range(n//d)) - (n // d) * y[i,j] <= (n // d) - 1) 
             mod.addConstr(quicksum(x[i + k*d] for k in range(n//d)) - (n // d) * y[i,j] >= 0) 
    
    
             mod.addConstr(quicksum(1 - x[i + k*d] for k in range(n//d)) - (n // d) * z[i,j] <= (n // d) - 1)
             mod.addConstr(quicksum(1 - x[i + k*d] for k in range(n//d)) - (n // d) * z[i,j] >= 0)              
            
             
             mod.addConstr(y[i,j] + z[i,j] == U[i,j])
                   
         mod.addConstr(quicksum(U[i,j] for i in range(d)) <= d-1)
            
    
    
    # Limit how many solutions to collect
    mod.setParam(GRB.Param.PoolSolutions, 2000000000)
    # Limit the search space by setting a gap for the worst possible solution
    # that will be accepted
    mod.setParam(GRB.Param.PoolGap, 0.02)
    # do a systematic search for the k-best solutions
    mod.setParam(GRB.Param.PoolSearchMode, 2)

    mod.optimize()

    # Status checking
    status = mod.Status
    if status in (GRB.INF_OR_UNBD, GRB.INFEASIBLE, GRB.UNBOUNDED):
        print('The model cannot be solved because it is infeasible or '
              'unbounded')
        return None

    if status != GRB.OPTIMAL:
        print('Optimization was stopped with status ' + str(status))
        return None

    # Print number of solutions stored
    nSolutions = mod.SolCount
    print('Number of solutions found: ' + str(nSolutions), 'runtime:', mod.runtime)

    # Print objective values of solutions
    # Ls = []
    # for e in range(nSolutions):
    #     mod.setParam(GRB.Param.SolutionNumber, e)
    #     # print('pattern:', e, [j for j in lam if lam[j].Xn > 1e-09])
    #     Ls.append([j for j in x if x[j].Xn > 1e-09])
        # print('x:', e, [x[j].Xn for j in x if x[j].Xn > 1e-09])

    # for l in sorted(Ls):
    #     print(l)
    print(nSolutions)
    # lambar = mod.getAttr('x', lam)

    #

ALLTEST = [
	{'t': 1, 'n': 72, 'D': [24,36], 'A': [0,8,16,18,26,34]},
	{'t': 2, 'n': 108, 'D': [35,54], 'A': [0, 12, 24, 27, 39, 51]},
	{'t': 3, 'n': 120, 'D': [24,40,60], 'A': [0, 8, 16, 30, 38, 46]},
	{'t': 4, 'n': 120, 'D': [24,40,60], 'A': [0, 8, 16, 24, 30, 32, 38, 46, 54, 62]},
	{'t': 5, 'n': 144, 'D': [48,72], 'A': [0, 16, 18, 32, 34, 50]},
	{'t': 6, 'n': 144, 'D': [48,72], 'A': [0, 16, 32, 36, 52, 68]},
	{'t': 7, 'n': 144, 'D': [48,72], 'A': [0, 9, 16, 25, 32, 36, 41, 45, 52, 61, 68, 77]},
	{'t': 8, 'n': 144, 'D': [48,72], 'A': [0, 16, 18, 32, 34, 36, 50, 52, 54, 68, 70, 86]},
	{'t': 9, 'n': 168, 'D': [24,56,84], 'A': [0, 8, 16, 42, 50, 58]},
	{'t': 10, 'n': 168, 'D': [24,56,84], 'A': [0, 8, 16, 24, 32, 40, 42, 48, 50, 58, 66, 74, 82, 90]},
	{'t': 11, 'n': 180, 'D': [36,60,90], 'A': [0, 12, 24, 45, 57, 69]},
	{'t': 12, 'n': 180, 'D': [36,60,90], 'A': [0, 18, 20, 38, 40, 58]},
	{'t': 13, 'n': 180, 'D': [36,60,90], 'A': [0, 12, 24, 36, 45, 48, 57, 69, 81, 93]},
	{'t': 14, 'n': 180, 'D': [36,60,90], 'A': [0, 18, 20, 36, 38, 40, 54, 56, 58, 72, 74, 76, 92, 94, 112]},
	{'t': 15, 'n': 180, 'D': [36,60,90], 'A': [0, 20, 40, 45, 65, 85]},
	{'t': 16, 'n': 420, 'D': [60,84,140,210], 'A': [0, 12, 24, 36, 48, 70, 82, 94, 106, 118]},
	{'t': 17, 'n': 420, 'D': [60,84,140,210], 'A': [0, 12, 24, 36, 48, 60, 70, 72, 82, 94, 106, 118, 130, 142]},
	{'t': 18, 'n': 420, 'D': [60,84,140,210], 'A': [0, 12, 24, 36, 48, 70, 82, 94, 106, 118, 140, 152, 164, 176, 188]},
	{'t': 19, 'n': 420, 'D': [60,84,140,210], 'A': [0, 12, 24, 36, 48, 60, 70, 72, 82, 94, 106, 118, 130, 140, 142, 152, 164, 176, 188, 200, 212]},
	{'t': 20, 'n': 420, 'D': [60,84,140,210], 'A': [0, 20, 40, 42, 62, 82]},
	{'t': 21, 'n': 420, 'D': [60,84,140,210], 'A': [0, 20, 40, 42, 60, 62, 80, 82, 100, 102, 120, 122, 142, 162]},
	{'t': 22, 'n': 420, 'D': [60,84,140,210], 'A': [0, 20, 40, 42, 62, 82, 84, 104, 124, 126, 146, 166, 168, 188, 208]},
	{'t': 23, 'n': 420, 'D': [60,84,140,210], 'A': [0, 20, 40, 42, 60, 62, 80, 82, 84, 100, 102, 104, 120, 122, 124, 126, 142, 144, 146, 162, 164, 166, 168, 184, 186, 188, 204, 206, 208,  226, 228, 246, 248, 268, 288]},
	{'t': 24, 'n': 420, 'D': [60,84,140,210], 'A': [0, 28, 30, 56, 58, 86]},
	{'t': 25, 'n': 420, 'D': [60,84,140,210], 'A': [0, 28, 30, 56, 58, 60, 86, 88, 90, 116, 118, 120, 146, 148, 150, 176, 178, 180, 206, 208, 236]},
	{'t': 26, 'n': 420, 'D': [60,84,140,210], 'A': [0, 28, 30, 56, 58, 84, 86, 112, 114, 142]},
	{'t': 27, 'n': 420, 'D': [60,84,140,210], 'A': [0, 28, 30, 56, 58, 60, 84, 86, 88, 90, 112, 114, 116, 118, 120, 142, 144, 146, 148, 150, 172,174, 176, 178, 180, 202, 204, 206, 208, 232, 234, 236, 262, 264, 292]},
	{'t': 28, 'n': 420, 'D': [60,84,140,210], 'A': [0, 12, 24, 36, 48, 60, 72, 105, 117, 129, 141, 153, 165, 177]},
	{'t': 29, 'n': 900, 'D': [180,300,450], 'A': [0, 18, 100, 118, 200, 218]},
	{'t': 30, 'n': 900, 'D': [180,300,450], 'A': [0, 18, 36, 54, 72, 90, 100, 108, 118, 126, 136, 144, 154, 162, 
                                                   172, 180, 190, 198, 200, 208, 216, 218, 226, 234, 236, 244, 252, 254, 262, 
                                                    270, 272, 280, 288, 290, 298, 306, 308, 316, 324, 326, 334, 342, 344,
                                                    352, 360, 362, 370, 378, 380, 388, 396, 398, 406, 414, 416, 424, 432, 
                                                    434, 442, 452, 460, 470, 478, 488, 496, 506, 514, 524, 532, 542, 560, 
                                                    578, 596, 614, 632]},
	{'t': 31, 'n': 900, 'D': [180,300,450], 'A': [0, 18, 36, 54, 72, 100, 118, 136, 154, 172, 200, 218, 236, 254, 272]},
	{'t': 32, 'n': 900, 'D': [180,300,450], 'A': [0, 18, 36, 54, 72, 90, 100, 108, 118, 126, 136, 144, 154, 162, 172, 
                                               190, 200, 208, 218, 226, 236, 244, 254, 262, 272, 290, 308, 326, 344, 362]},
	{'t': 33, 'n': 900, 'D': [180,300,450], 'A': [0, 36, 50, 72, 86, 108, 122, 144, 158, 194]},
	{'t': 34, 'n': 900, 'D': [180,300,450], 'A': [0, 36, 50, 72, 86, 100, 108, 122, 136, 144, 150, 158, 172, 186, 194, 
                                                   200, 208, 222, 236, 244, 250, 258, 272, 286, 294, 300, 308, 322, 336, 
                                                   344, 350, 358, 372, 386, 394, 400, 408, 422, 436, 444, 458, 472, 494,  508, 544]},
	{'t': 35, 'n': 900, 'D': [180,300,450], 'A': [0, 36, 50, 72, 86, 100, 108, 122, 136, 144, 158, 172, 194, 208, 244]},
	{'t': 36, 'n': 900, 'D': [180,300,450], 'A': [0, 36, 50, 72, 86, 100, 108, 122, 136, 144, 150, 158, 172, 186, 
                                               194, 200, 208, 222, 236, 244, 250, 258, 272, 286, 294, 308, 322, 344, 358, 394]},
]

if __name__ == '__main__':
	for test in ALLTEST:
		ILP(test['A'], test['n'], test['D'], test['t'])

	logging.shutdown()