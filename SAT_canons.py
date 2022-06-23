

from pysat.pb import *
from pysat.formula import CNF
from pysat.solvers import Maplesat

from time import perf_counter

import itertools

import logging

logging.basicConfig(
	filename='pysat_alltest.log',
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def BuildCNF(A, n, D, num_test):
    
    cnf = CNF()
    
    for i in range(n):
           A_1 = [((n + i -j)%n) + 1 for j in A]
           cnf.append([j for j in A_1])
           for comb in itertools.combinations(A_1, 2):
                   cnf.append([-j for j in comb])
    cnf.append([1])
                 
    U = {}
    c = n+1
    for j, d in enumerate(D):
        for i in range(d):
            U[i,j] = c
            c = c+1
            
    for j, d in enumerate(D):
        for i in range(d):
            
            a = [(i + k*d) + 1 for k in range(n//d)]
            b = [-((i + k*d) + 1) for k in range(n//d)]
            
            cnf.extend([a +[U[i,j]],b + [U[i,j]]])
            
            
            cnf.extend([[x,y] + [-U[i,j]] for (x,y) in list(itertools.product(a,b))])
        cnf.append([-(U[i,j]) for i in range(d)])
    
    for j in range(n-len(A), n):
         cnf.append([-(j+1)])
    
    g = Maplesat() 
    g.append_formula(cnf)

    def interrupt():
        g.interrupt()
  
    t0 = perf_counter()
    g.solve()

    num = 0
    for m in g.enum_models():
        num = num+1
    runtime = perf_counter() - t0
    
    logging.info('canons: {}, n: {}, t: {}, time: {}'.format(num, n, num_test, round(runtime, 4)))


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
		BuildCNF(test['A'], test['n'], test['D'], test['t'])

	logging.shutdown()