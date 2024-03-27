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
p1=[i/100 for i in range(60,72,5)]
def simu(p1):
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

    init=[]
    for i in range(200):
      m1 = np.random.normal(0.5,0.1)
      ini=[100*m1,100*(1-m1)*0.5,100*(1-m1)*0.5]
      init.append(ini)
    #原始模型
    n_jobs = 200
    N0 = 100
    crOA = []
    crOB = []
    crOC = []
    Nt_Obreed = []
    T = []
    for i in range(n_jobs):
        #m1 = np.random.normal(0.5,0.1)
        num_elements = 3
        system = Gillespie(
            num_elements,
            inits=init[i],
            max_cell_num=9000,
            max_t=40,

        )
        r1=0.6
        r2=0.3
        p2=0.4
        d=0.01

        paa = lambda t: r1*p1
        pab = lambda t: r1*(1-p1)
        pa0 = lambda t: d
        pbb = lambda t: r2*p2
        pbc = lambda t: r2*(1-p2)
        pb0 = lambda t: d
        pc0 = lambda t: d


        system.add_reaction(paa, [1, 0, 0], [2, 0, 0])
        system.add_reaction(pab, [1, 0, 0], [0, 2, 0])
        system.add_reaction(pa0, [1, 0, 0], [0, 0, 0])
        system.add_reaction(pbb, [0, 1, 0], [0, 2, 0])
        system.add_reaction(pbc, [0, 1, 0], [0, 0, 2])
        system.add_reaction(pb0, [0, 1, 0], [0, 0, 0])
        system.add_reaction(pc0, [0, 0, 1], [0, 0, 0])

        system.evolute(20000000)
        t = system.t
        cell_num_traj = np.array(system.n)

        c0 = cell_num_traj[:, 0]
        c1 = cell_num_traj[:, 1]
        c2 = cell_num_traj[:, 2]
        Nt_breed_1 = c0 + c1+c2
        cr_A = c0 / (c0 + c1+c2)
        cr_B = c1 / (c0 + c1+c2)
        cr_C = c2 / (c0 + c1+c2)
        crOA.append(cr_A)
        crOB.append(cr_B)
        crOC.append(cr_C)
        T.append(t)
        Nt_Obreed.append(Nt_breed_1)

    def find_i(t, t0):
        i = 0
        while i < len(t):
            if t[i] < t0:
                i = i + 1
            else:
                break
        return i

    huizong_OA = []
    huizong_OB = []
    huizong_OC = []
    huizongO = []
    Run_Otime = []
    huizong1_len = 0
    huizong2_len = 0
    huizong3_len = 0
    Run_Otime_len = 0
    huizong_OA_str = 'huizong_OA_'
    huizong_OB_str = 'huizong_OB_'
    huizong_OC_str = 'huizong_OC_'
    huizongO_str = 'huizongO_'
    Run_Otime_str = 'Run_Otime_'
    t_max = min([max(T[k]) for k in range(len(T))]) - 1
    for j in np.arange(0, t_max, 0.1):
        # print('T=',j)
        xxx1 = []
        xxx2 = []
        xxx3 = []
        z = []
        time_once = []
        for k in range(200):
            m = find_i(T[k], j)
            xxx1.append(crOA[k][m])
            xxx2.append(crOB[k][m])
            xxx3.append(crOC[k][m])
            z.append(Nt_Obreed[k][m])
            time_once.append(T[k])
            print("time_once", "{:^10d}", sep=':', end='\n', flush=True)

        huizong_OA.append(xxx1)
        huizong_OB.append(xxx2)
        huizong_OC.append(xxx3)
        huizongO.append(z)
        Run_Otime.append(time_once)
        del xxx1
        del xxx2
        del xxx3
        del z
        del time_once
        huizong_OA_len = len(huizong_OA)
        print("huizong1_len", "{:^10d}", sep=':', end='\n', flush=True)
        huizong_OB_len = len(huizong_OB)
        print("huizong2_len", "{:^10d}", sep=':', end='\n', flush=True)
        huizong_OC_len = len(huizong_OC)
        print("huizong3_len", "{:^10d}", sep=':', end='\n', flush=True)
        huizongO_len = len(huizongO)
        print("huizong_len", "{:^10d}", sep=':', end='\n', flush=True)
        Run_Otime_len = len(Run_Otime)
        print("Run_time_len", "{:^10d}", sep=':', end='\n', flush=True)
    np.savetxt(huizong_OA_str + str(p1), huizong_OA)
    np.savetxt(huizong_OB_str + str(p1), huizong_OB)
    np.savetxt(huizong_OC_str + str(p1), huizong_OC)
    np.savetxt(huizongO_str + str(p1), huizongO)
    # np.savetxt(Run_time_str+str(delta),t_max)
    # del Run_time
    # del huizong1
    # del huizong2
    # del huizong3
    return 0
    
    #简化之后的模型
    n_jobs = 200
    N0 = 100
    crSA = []
    crSB = []
    Nt_Sbreed = []
    T = []
    for i in range(n_jobs):
        #m1 = np.random.normal(0.5,0.1)
        num_elements = 2
        system = Gillespie(
            num_elements,
            inits=[init[i][0],init[i][1]+init[i][2]],
            max_cell_num=9000,
            max_t=40,

        )
        r1=0.6
        r2=0.3
        d=0.01
        b=2*r1*p1*(1-p1)-2*r1*(1-p1)*(1-p1)
        c=2*r1*p1*(1-p1)-2*r1*(1-p1)*(1-p1)+4*r2*p1*(1-p1)

        paa = lambda t: r1*p1
        paq = lambda t: r1*(1-p1)
        pa0 = lambda t: d
        pqq = lambda t: r2*(b/c)
        pq0 = lambda t: 2*d
    
    
        system.add_reaction(paa, [1, 0], [2, 0])
        system.add_reaction(paq, [1, 0], [0, 2])
        system.add_reaction(pa0, [1, 0], [0, 0])
        system.add_reaction(pqq, [0, 1], [0, 2])
        system.add_reaction(pq0, [0, 1], [0, 0])
        system.evolute(20000000)
        t = system.t
        cell_num_traj = np.array(system.n)

        c0 = cell_num_traj[:, 0]
        c1 = cell_num_traj[:, 1]
        Nt_breed_1 = c0 + c1
        cr_A = c0 / (c0 + c1)
        cr_B = c1 / (c0 + c1)
        crSA.append(cr_A)
        crSB.append(cr_B)
        T.append(t)
        Nt_Sbreed.append(Nt_breed_1)

    def find_i(t, t0):
        i = 0
        while i < len(t):
            if t[i] < t0:
                i = i + 1
            else:
                break
        return i

    huizong_SA = []
    huizong_SB = []
    huizongS = []
    Run_timeS = []
    huizong1_len = 0
    huizong2_len = 0
    Run_time_len = 0
    huizong_SA_str = 'huizong_SA_'
    huizong_SB_str = 'huizong_SB_'
    huizongS_str = 'huizongS_'
    Run_timeS_str = 'Run_timeS_'
    t_max = min([max(T[k]) for k in range(len(T))]) - 1
    for j in np.arange(0, t_max, 0.1):
        # print('T=',j)
        xxx1 = []
        xxx2 = []

        z = []
        time_once = []
        for k in range(200):
            m = find_i(T[k], j)
            xxx1.append(crA[k][m])
            xxx2.append(crB[k][m])
            z.append(Nt_breed[k][m])
            time_once.append(T[k])
            print("time_once", "{:^10d}", sep=':', end='\n', flush=True)

        huizong_SA.append(xxx1)
        huizong_SB.append(xxx2)
        huizongS.append(z)
        Run_timeS.append(time_once)
        del xxx1
        del xxx2
        del z
        del time_once
        huizong_SA_len = len(huizong_SA)
        print("huizong1_len", "{:^10d}", sep=':', end='\n', flush=True)
        huizong_SB_len = len(huizong_SB)
        print("huizong2_len", "{:^10d}", sep=':', end='\n', flush=True)
        huizong_Slen = len(huizong)
        print("huizong_len", "{:^10d}", sep=':', end='\n', flush=True)
        Run_timeS_len = len(Run_timeS)
        print("Run_time_len", "{:^10d}", sep=':', end='\n', flush=True)
    np.savetxt(huizong_SA_str + str(p1), huizong_SA)
    np.savetxt(huizong_SB_str + str(p1), huizong_SB)
    np.savetxt(huizongS_str + str(p1), huizongS)
    # np.savetxt(Run_time_str+str(delta),t_max)
    # del Run_time
    # del huizong1
    # del huizong2
    # del huizong3
    return 0

p = multiprocessing.Pool(5)
b = p.map(simu, p1)
p.close()
p.join()
