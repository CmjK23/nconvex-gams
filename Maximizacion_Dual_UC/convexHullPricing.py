#%%
import os
import numpy as np
import pyomo.environ as pyo
from scipy.optimize import minimize
from functools import partial

from modelos_esclavos import crear_esclavos_generacion, crear_esclavo_red
from modelo_UC_primal import modeloUnitCommitment

#%% Funciones auxiliares del maestro
def excdt_rule(m, node):
    return sum(
        m.pi[t, node]*m.Demanda[t, node]
        for t in m.t
    )

def congestion_rule(m, i, j):
    return sum(
        -m.Fmax[i, j]*(m.mu[t, i, j] + m.nu[t, i, j])
        for t in m.t
    )
#%% Problema maestro
def crear_maestro(primal, price_cap=150):
    '''
    Crea el problema maestro del dual a partir de su problema primal
    '''
    master = pyo.ConcreteModel()

    # ---------------- SETS ----------------

    master.t = pyo.Set(initialize=list(primal.t), ordered=True)

    master.N = pyo.Set(initialize=list(primal.N), ordered=True)

    master.linea = pyo.Set(
        initialize=list(primal.linea),
        dimen=2
    )

    master.ut = pyo.Set(
        initialize=list(primal.ut),
        ordered=True
    )

    master.map = pyo.Set(
        initialize=list(primal.map),
        dimen=2
    )

    # ---------------- SUBCONJUNTOS AUXILIARES ----------------

    master.NodesIn = pyo.Set(
        master.N,
        initialize=lambda m, node:
            [i for i, j in m.linea if j == node]
    )

    master.NodesOut = pyo.Set(
        master.N,
        initialize=lambda m, node:
            [j for i, j in m.linea if i == node]
    )

    # ---------------- PARÁMETROS ----------------

    master.B = pyo.Param(
        master.linea,
        initialize={(i, j): pyo.value(primal.B[i, j])
                    for (i, j) in primal.linea}
    )

    master.Fmax = pyo.Param(
        master.linea,
        initialize={(i, j): pyo.value(primal.Fmax[i, j])
                    for (i, j) in primal.linea}
    )

    master.Demanda = pyo.Param(
        master.t,
        master.N,
        initialize={(t, n): pyo.value(primal.Demanda[t, n])
                    for t in primal.t
                    for n in primal.N}
    )

    # ---------------- VARIABLES DUALES ----------------

    master.pi = pyo.Var(
        master.t,
        master.N,
        bounds=(0, price_cap)
    )

    master.mu = pyo.Var(
        master.t,
        master.linea,
        domain=pyo.NonNegativeReals
    )

    master.nu = pyo.Var(
        master.t,
        master.linea,
        domain=pyo.NonNegativeReals
    )

    master.lambda_ = pyo.Var(
        master.t,
        master.linea,
        domain=pyo.Reals
    )

    master.theta = pyo.Var(
        master.ut,
        domain=pyo.NonNegativeReals
    )

    # ---------------- EXPRESIONES ----------------

    master.excdtTotalN = pyo.Expression(
        master.N,
        rule=excdt_rule
    )

    master.congestionL = pyo.Expression(
        master.linea,
        rule=congestion_rule
    )

    # ---------------- OBJETIVO ----------------

    master.obj = pyo.Objective(
        expr=
        sum(master.excdtTotalN[n] for n in master.N)
        + sum(master.congestionL[i, j]
              for i, j in master.linea)
        - sum(master.theta[ut] for ut in master.ut),
        sense=pyo.maximize
    )
    
    # --------------- RESTRICCIONES --------------
    master.constr1 = pyo.Constraint(master.t*master.linea,
                                   rule=lambda m,t,i,j:
                                   m.pi[t,i]-m.pi[t,j]+m.mu[t,i,j]\
                                    -m.nu[t,i,j]+m.lambda_[t,i,j]==0
                                   )
    master.constr2 =pyo.Constraint(master.t*master.N,rule=lambda m,t,nodo:
                                   sum(m.lambda_[t,i,nodo]*m.B[i,nodo] for i in m.NodesIn[nodo])\
                                   -sum(m.lambda_[t,nodo,j]*m.B[nodo,j] for j in m.NodesOut[nodo])\
                                   == 0
                                   )

    # ---------------- CORTES ----------------

    master.cuts = pyo.ConstraintList()

    return master


#%%
class CHPframework(object):
    
    def __init__(self,filename,solver='highs',price_cap=150):

        abstract_modeloUC = modeloUnitCommitment()

        data = pyo.DataPortal(model=abstract_modeloUC)
        data.load(filename=filename)
        self.modeloUC = abstract_modeloUC.create_instance(data)

        for t in self.modeloUC.t:
            self.modeloUC.ang[t,self.modeloUC.N.first()].fix(0)

        self.price_cap = price_cap
        self.master = crear_maestro(self.modeloUC,self.price_cap)
        self.esclavos_gen = crear_esclavos_generacion(self.modeloUC)
        self.esclavo_red = crear_esclavo_red(self.modeloUC)
        self.solver = pyo.SolverFactory(solver)

    def optimizar(self,alpha=0.2,pi_0=30,epsilon=1e-6,iter=200):
        
        index = {(t,nodo):pos for pos,(t,nodo)
                 in enumerate((t,nodo)
                 for nodo in self.modeloUC.N
                 for t in self.modeloUC.t)}
        
        Demanda = np.zeros(len(index))
        for (t,nodo), pos in index.items():
            Demanda[pos] = self.modeloUC.Demanda[t,nodo]

        map = list(self.modeloUC.map)
        N = list(self.modeloUC.N)
        T = self.modeloUC.t
        nodesOut = {nodo:list(Nodos) for nodo,Nodos in self.modeloUC.NodesOut.items()}
        nodesIn = {nodo:list(Nodos) for nodo,Nodos in self.modeloUC.NodesIn.items()}

        constraints = []
        bounds = [(0,self.price_cap)]*len(index)

        duales = []

        for k in range(iter):
            if k==0:
                pi_k = np.ones(len(index))*pi_0

            for nodo in N:
                for t in T:
                    self.esclavo_red.pi[t,nodo] = pi_k[index[t,nodo]]
            #Solucion de la red
            self.solver.solve(self.esclavo_red)

            for nodo, ut in map:
                #Fijar precios antes de resolver los esclavos
                for t in T:
                    self.esclavos_gen[ut].pi[t] = pi_k[index[t,nodo]]
                #Solucion de problemas esclavos
                self.solver.solve(self.esclavos_gen[ut])

            flujo_neto = np.zeros(len(index))
            for t in T:
                for nodo in N:
                    flujo_neto[index[t,nodo]] = pyo.value(\
                    sum(self.esclavo_red.flujo[t,i,nodo] for i in nodesIn[nodo]) \
                    - sum(self.esclavo_red.flujo[t,nodo,j] for j in nodesOut[nodo])
                    )

            gen_nodo = np.zeros(len(index))
            for t in T:
                for nodo, ut in map:
                    gen_nodo[index[t,nodo]] += pyo.value(self.esclavos_gen[ut].p[t])
            
            #Calculo del dual
            L_k = np.dot(Demanda,pi_k) + pyo.value(self.esclavo_red.obj) \
                  - sum(pyo.value(self.esclavos_gen[ut].profit) for ut in self.modeloUC.ut)
            
            duales.append(L_k)

            #Calculo de supergradiente
            s = Demanda - gen_nodo - flujo_neto

            #Calculo de Upper Bound (UB)
            for nodo, ut in map:
                self.master.cuts.add(sum(pyo.value(self.esclavos_gen[ut].p[t])*self.master.pi[t,nodo] for t in T)\
                                     - pyo.value(self.esclavos_gen[ut].costo) <= self.master.theta[ut])
            
            self.solver.solve(self.master)

            UB = pyo.value(self.master.obj)
            
            #Lower Bound
            LB = max(duales)

            #Level
            level_k = alpha*UB + (1-alpha)*LB

            #Modelo de proyeccion

            obj = lambda x, pi_ref=pi_k.copy(): \
                np.sum((x - pi_ref)**2)

            cut_k = lambda x, level, pi_k=pi_k.copy(), s=s.copy(), L_k=L_k: \
                np.dot(x-pi_k, s) + L_k - level

            constraints.append(cut_k)

            #Solucion del problema de proyeccion

            result = minimize(
                fun=obj,
                x0=pi_k,
                method='SLSQP',
                bounds=bounds,
                constraints=[{'type':'ineq','fun':partial(f,level=level_k)} for f in constraints]
            )

            pi_k = result.x
            
            gap = (UB - LB) / abs(UB) if abs(UB) > epsilon else 0
            print(f"Iter {k}: LB={LB:.6f}, UB={UB:.6f}, Gap={gap:.6f}")
            
            if gap <= epsilon:
                break
            

         
# %%
test = CHPframework('test_multinodal.dat')
test.optimizar(iter=100)
# %%
