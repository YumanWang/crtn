import multiprocessing
import operator
import numpy as np
import math
# gillespie算法
import os
from collections import defaultdict
from copy import deepcopy
import matplotlib.pyplot as plt
from scipy.special import comb
from tqdm.autonotebook import tqdm  # 关于进度条
delta2=[i/100 for i in range(15,50,10)]
def simu(delta2):
    class Reaction:
        '''
        Cell division/differentiation type

        Args:
            rate:
                reaction rate function
            num_lefts:
                Cell numbers before reaction
            num_right:
                Cell numbers after reaction
            index:
                Reaction index
        '''

        def __init__(self, rate: callable = None, num_lefts: list = None, num_rights: list = None, index: int = None):
            self.rate = rate
            assert len(num_lefts) == len(num_rights)
            self.num_lefts = np.array(num_lefts)
            self.num_rights = np.array(num_rights)
            self.num_diff = self.num_rights - self.num_lefts
            self.index = index
            if 2 * sum(num_lefts) == sum(num_rights):
                self.type = "proliferate"
            else:
                self.type = "differentiation"

        def combine(self, n, s):
            return np.prod(comb(n, s))

        def propensity(self, n, t):
            return self.rate(t) * self.combine(n, self.num_lefts)


    class Gillespie:
        '''
        Gillespie simulation

        Args:
            num_elements:
                Cell type number
            inits:
                Initial cell number
            max_cell_num:
                Maximum cell number
        '''

        def __init__(
                self,
                num_elements: int,
                inits: list = None,
                max_cell_num: int = 20000,
                max_t: float = 1000
        ):

            assert num_elements > 0
            self.num_elements = num_elements
            self.max_t = max_t
            self.reactions = []
            if inits is None:
                self.n = [np.ones(self.num_elements)]
            else:
                self.n = [np.array(inits)]

            self.generation_time = [0]
            self.max_cell_num = max_cell_num

        def add_reaction(self, rate: callable = None, num_lefts: list = None, num_rights: list = None, index: int = None):
            '''
            Add reactions to simulation

            Args:
                rate:
                    reaction rate function
                num_lefts:
                    Cell numbers before reaction
                num_right:
                    Cell numbers after reaction
                index:
                    Reaction index
            '''
            assert len(num_lefts) == self.num_elements
            assert len(num_rights) == self.num_elements
            self.reactions.append(Reaction(rate, num_lefts, num_rights, index))

        def evolute(self, steps: int):
            '''
            Run simulation

            Args:
                steps:
                    How many steps to evolute before step
            '''
            self.t = [0]
            self.log = []

            #         with tqdm(total=self.max_cell_num) as pbar:

            for _ in range(steps):
                all_cell_num = sum(self.n[-1])
                #                 pbar.update(all_cell_num - pbar.n)
                #                 if all_cell_num > self.max_cell_num:
                #                     print("\n maximum cell number reached")
                #                     break
                if self.t[-1] >= self.max_t:
                    break

                A = np.array(
                    [
                        rec.propensity(self.n[-1], self.t[-1])
                        for rec in self.reactions
                    ]
                )

                A0 = A.sum()
                A /= A0
                t0 = -np.log(np.random.random()) / A0
                self.t.append(self.t[-1] + t0)
                react = np.random.choice(self.reactions, p=A)

                self.log.append(react.index)
                self.n.append(self.n[-1] + react.num_diff)

    import matplotlib.pyplot as plt
    import numpy as np
    n_jobs = 100
    N0=100
    crA = []
    crB = []
    Nt_breed = []
    T = []
    for i in range(n_jobs):
        m_1=np.random.normal(0.5,0.1)
        num_elements = 2
        system = Gillespie(
            num_elements,
            inits=[N0*m_1, N0*(1-m_1)],
            max_cell_num=9000,
            max_t=40,

        )
        r1 = 0.4
        r2 = 0.2
        p1 = 0.25
        #delta2 = 0.25
    
        paa = lambda t: r1*p1
        pab = lambda t: r1*(1-p1)
        pbb = lambda t: r2*(1-delta2)
        pba = lambda t: r2*delta2
    
    
        system.add_reaction(paa, [1, 0], [2, 0])
        system.add_reaction(pab, [1, 0], [0, 2])
        system.add_reaction(pbb, [0, 1], [0, 2])
        system.add_reaction(pba, [0, 1], [2, 0])

        system.evolute(20000000)
        t = system.t
        cell_num_traj = np.array(system.n)

        c0 = cell_num_traj[:, 0]
        c1 = cell_num_traj[:, 1]
        Nt_breed_1 = c0 + c1
        cr = c0 / (c0 + c1)
        crA.append(cr)
        crB.append(1 - cr)
        T.append(t)
        Nt_breed.append(Nt_breed_1)

    def find_i(t, t0):
        i = 0
        while i < len(t):
            if t[i] < t0:
                i = i + 1
            else:
                break
        return i

    huizong_A = []
    huizong_B = []
    huizong = []
    Run_time = []
    huizong1_len = 0
    huizong2_len = 0
    huizong3_len = 0
    Run_time_len = 0
    huizong_A_str = 'huizong_A_'
    huizong_B_str = 'huizong_B_'
    huizong_str = 'huizong_'
    Run_time_str = 'Run_time_'
    t_max = min([max(T[k]) for k in range(len(T))]) - 1
    for j in np.arange(0, t_max, 0.1):
        # print('T=',j)
        xxx1 = []
        xxx2 = []
        z = []
        time_once = []
        for k in range(100):
            m = find_i(T[k], j)
            xxx1.append(crA[k][m])
            xxx2.append(crB[k][m])
            z.append(Nt_breed[k][m])
            time_once.append(T[k])
            print("time_once", "{:^10d}", sep=':', end='\n', flush=True)

        huizong_A.append(xxx1)
        huizong_B.append(xxx2)
        huizong.append(z)
        Run_time.append(time_once)
        del xxx1
        del xxx2
        del z
        del time_once
        huizong_A_len = len(huizong_A)
        print("huizong1_len", "{:^10d}", sep=':', end='\n', flush=True)
        huizong_B_len = len(huizong_B)
        print("huizong2_len", "{:^10d}", sep=':', end='\n', flush=True)
        huizong_len = len(huizong)
        print("huizong3_len", "{:^10d}", sep=':', end='\n', flush=True)
        Run_time_len = len(Run_time)
        print("Run_time_len", "{:^10d}", sep=':', end='\n', flush=True)
    np.savetxt(huizong_A_str + str(delta2), huizong_A)
    np.savetxt(huizong_B_str + str(delta2), huizong_B)
    np.savetxt(huizong_str + str(delta2), huizong)
    # np.savetxt(Run_time_str+str(delta),t_max)
    # del Run_time
    # del huizong1
    # del huizong2
    # del huizong3
    return 0


p = multiprocessing.Pool(5)
b = p.map(simu, delta2)
p.close()
p.join()
