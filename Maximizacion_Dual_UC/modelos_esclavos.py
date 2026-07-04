import pyomo.environ as pyo
import numpy as np
#--------------------MODELO DE MAXIMIZACION DE BENEFICIO PARA TERMICAS------------------------#
#Funciones auxiliares
def descartar_t_inicial(m):
    return list(m.t)[1:]
def costo_exp(model):
    return sum([model.y[t]*model.c_arr
               +model.u[t]*model.coef1
               +model.p[t]*model.coef2 for t in model.t])
def ingresos_exp(model):
    return sum([model.p[t]*model.pi[t] for t in model.t]) 

def cond_lim_op1(m,t): #Condicion limite operativa
    return m.pmin*m.u[t] <= m.p[t]

def cond_lim_op2(m,t): #Condicion limite operativa
    return m.p[t] <= m.pmax*m.u[t]
def cond_op1(m,t):
    if t==m.t.first():
        return m.u[t]-m.u0== m.y[t]- m.w[t]
    else:
        t_prev = m.t.prev(t)
        return m.u[t]-m.u[t_prev]== m.y[t]- m.w[t]

def cond_op2(m,t):
    return m.y[t] <= m.u[t]

def cond_op3(m,t):
    if t==m.t.first():
        return m.y[t] <= 1 - m.u0
    else:
        t_prev = m.t.prev(t)
        return m.y[t] <= 1 - m.u[t_prev]

def cond_op4(m,t):
    if t==m.t.first():
        return m.w[t] <= m.u0
    else:
        t_prev = m.t.prev(t)
        return m.w[t] <= m.u[t_prev]

def cond_op5(m,t):
    return m.w[t] <= 1 - m.u[t]
def startup(m,t): #Restriccion de arranque
    if t==m.t.first():
        return m.p[t] - m.p0 <= m.ramp_n*m.u0 + m.ramp_up*m.y[t]
    else:
        t_prev = m.t.prev(t)
        return m.p[t] - m.p[t_prev] <= m.ramp_n*m.u[t_prev] + m.ramp_up*m.y[t]

def shutdown(m,t): #Restriccion de desconexion
    if t==m.t.first():
        return m.p0 - m.p[t] <= m.ramp_n*m.u[t] + m.ramp_dn*m.w[t]
    else:
        t_prev = m.t.prev(t)
        return m.p[t_prev] - m.p[t] <= m.ramp_n*m.u[t] + m.ramp_dn*m.w[t]
def ineq_arrd_max(m):
    return sum(m.y[t] for t in m.t) <= m.arrd_max

def ineq_pard_max(m):
    return sum(m.w[t] for t in m.t) <= m.pard_max
def ecComb(m): #Ecuacion de reserva de combustible
    return sum([m.p[t]*m.heat_rate + m.u[t]*m.no_load_consumption for t in m.t]) <= m.dComb
def ectminUp(m,t):
    T = list(m.t)
    idx = T.index(t)
    start = max(0,idx-m.tminUp+1)
    return sum(m.y[tau] for tau in T[start:idx+1]) <= m.u[t]

def ectminDn(m,t):
    T = list(m.t)
    idx = T.index(t)
    start = max(0,idx-m.tminDn+1)
    return sum(m.w[tau] for tau in T[start:idx+1]) <= 1 - m.u[t]
def modeloTermica():
    #creacion de modelo abstracto
    model = pyo.AbstractModel()
    #SETS
    model.t = pyo.RangeSet(1,24)
    model.t_rest = pyo.Set(within=model.t,initialize=descartar_t_inicial)

    #PARAMETROS
    #Estado inmediatamente antes de empezar el estudio
    model.u0 = pyo.Param(domain=pyo.Binary)
    model.p0 = pyo.Param(domain=pyo.NonNegativeReals)
    #Cantidades limites de generacion
    model.pmin = pyo.Param(domain=pyo.NonNegativeReals)
    model.pmax = pyo.Param(domain=pyo.NonNegativeReals,
                        validate=lambda m, val_pmax: val_pmax > m.pmin)
    #Precio asignado en una iteracion del problema maestro
    model.pi = pyo.Param(model.t,initialize=0,mutable=True)
    #Costo de combustible
    model.cComb = pyo.Param(domain=pyo.NonNegativeReals)
    #Tasa de calor o Heat rate
    model.heat_rate = pyo.Param(domain=pyo.NonNegativeReals)
    #No load consumption
    model.no_load_consumption = pyo.Param(domain=pyo.NonNegativeReals)
    #Combustible disponible
    model.dComb = pyo.Param(domain=pyo.NonNegativeReals)
    #Modelamiento del costo de operacion de termicas
    model.coef1 = pyo.Param(domain=pyo.NonNegativeReals)
    model.coef2 = pyo.Param(domain=pyo.NonNegativeReals)
    model.c_arr = pyo.Param(domain=pyo.NonNegativeReals)
    #Tiempos minimos de encencido y apagado
    model.tminUp = pyo.Param(domain=pyo.PositiveIntegers)
    model.tminDn = pyo.Param(domain=pyo.PositiveIntegers)
    #Cantidades maximas de arranques y paradas
    model.arrd_max = pyo.Param(domain=pyo.PositiveIntegers)
    model.pard_max = pyo.Param(domain=pyo.PositiveIntegers)
    #Parametros de restricciones de rampas
    model.ramp_n = pyo.Param(domain=pyo.NonNegativeReals)
    model.ramp_up = pyo.Param(domain=pyo.NonNegativeReals,
                              validate=lambda m, val_ramp_up: val_ramp_up > m.pmin)
    model.ramp_dn = pyo.Param(domain=pyo.NonNegativeReals,
                              validate=lambda m, val_ramp_dn: val_ramp_dn > m.pmin)

    #VARIABLES
    #Potencia generada en cada bloque horario        
    model.p = pyo.Var(model.t)
    #Estado de encendido/apagado
    model.u = pyo.Var(model.t,domain=pyo.Binary)
    #Arranque
    model.y = pyo.Var(model.t,domain=pyo.Binary)
    #Parada
    model.w = pyo.Var(model.t,domain=pyo.Binary)

    #Funcion de costo operativo
    model.costo = pyo.Expression(rule=costo_exp)

    #Ingresos
    model.ingresos = pyo.Expression(rule=ingresos_exp)

    #Variable objetivo: beneficio del generador
    model.profit = pyo.Objective(expr=model.ingresos-model.costo,sense=pyo.maximize)

    #RESTRICCIONES
    model.bounds1 = pyo.Constraint(model.t,rule=cond_lim_op1)
    model.bounds2 = pyo.Constraint(model.t,rule=cond_lim_op2)

    #Restricciones de operacion por arranques/paradas
    model.cond_op1 = pyo.Constraint(model.t,rule=cond_op1)
    model.cond_op2 = pyo.Constraint(model.t,rule=cond_op2)
    model.cond_op3 = pyo.Constraint(model.t,rule=cond_op3)
    model.cond_op4 = pyo.Constraint(model.t,rule=cond_op4)
    model.cond_op5 = pyo.Constraint(model.t,rule=cond_op5)

    model.startup = pyo.Constraint(model.t,rule=startup)
    model.shutdown = pyo.Constraint(model.t,rule=shutdown)

    model.ineq_arrd_max = pyo.Constraint(rule=ineq_arrd_max)
    model.ineq_pard_max = pyo.Constraint(rule=ineq_pard_max)

    model.ecComb = pyo.Constraint(rule=ecComb)

    #Tiempos minimos de encencido/apagado
    model.ectminUp = pyo.Constraint(model.t_rest,rule=ectminUp)
    model.ectminDn = pyo.Constraint(model.t_rest,rule=ectminDn)

    return model

def crear_esclavos_generacion(modeloUnitCommitment):

    termica_abs = modeloTermica()

    esclavos = {}

    for ut in modeloUnitCommitment.ut:

        data = {None:{

            'u0': {None: pyo.value(modeloUnitCommitment.u0[ut])},
            'p0': {None: pyo.value(modeloUnitCommitment.p0[ut])},

            'pmin': {None: pyo.value(modeloUnitCommitment.pmin[ut])},
            'pmax': {None: pyo.value(modeloUnitCommitment.pmax[ut])},

            'cComb': {None: pyo.value(modeloUnitCommitment.cComb[ut])},
            'heat_rate': {None: pyo.value(modeloUnitCommitment.heat_rate[ut])},
            'no_load_consumption':
                {None: pyo.value(modeloUnitCommitment.no_load_consumption[ut])},
            'dComb': {None: pyo.value(modeloUnitCommitment.dComb[ut])},

            'coef1': {None: pyo.value(modeloUnitCommitment.coef1[ut])},
            'coef2': {None: pyo.value(modeloUnitCommitment.coef2[ut])},
            'c_arr': {None: pyo.value(modeloUnitCommitment.c_arr[ut])},

            'tminUp': {None: pyo.value(modeloUnitCommitment.tminUp[ut])},
            'tminDn': {None: pyo.value(modeloUnitCommitment.tminDn[ut])},

            'arrd_max': {None: pyo.value(modeloUnitCommitment.arrd_max[ut])},
            'pard_max': {None: pyo.value(modeloUnitCommitment.pard_max[ut])},

            'ramp_n': {None: pyo.value(modeloUnitCommitment.ramp_n[ut])},
            'ramp_up': {None: pyo.value(modeloUnitCommitment.ramp_up[ut])},
            'ramp_dn': {None: pyo.value(modeloUnitCommitment.ramp_dn[ut])},
        }}

        esclavos[ut] = termica_abs.create_instance(data)

    return esclavos

#%%Funcion objetivo del modelo esclavo de red
def red_obj(m):
    s = 0
    for nodo in m.N:
        for t in m.t:
            s+= m.pi[t,nodo]*(sum(m.flujo[t,nodo,destino] 
                for destino in m.NodesOut[nodo])
                - sum(m.flujo[t,origen,nodo]
                for origen in m.NodesIn[nodo])
                )
    return s
#%%
def crear_esclavo_red(primal):
    
    #Definicion de sets
    model = pyo.ConcreteModel()

    model.t = pyo.Set(initialize=list(primal.t),ordered=True)

    model.N = pyo.Set(initialize=list(primal.N),ordered=True)

    model.linea = pyo.Set(initialize=list(primal.linea),ordered=True,dimen=2)
    model.NodesIn = pyo.Set(
        model.N,
        initialize=lambda m, node:
            [i for i, j in m.linea if j == node]
    )

    model.NodesOut = pyo.Set(
        model.N,
        initialize=lambda m, node:
            [j for i, j in m.linea if i == node]
    )

    #Parametros
    model.pi = pyo.Param(model.t*model.N,initialize=20,
                         domain=pyo.NonNegativeReals,mutable=True)
    model.Fmax = pyo.Param(model.linea,initialize={(i, j): pyo.value(primal.Fmax[i, j])
                                                   for (i, j) in primal.linea},
                                                   domain=pyo.NonNegativeReals)
    model.B = pyo.Param(
        model.linea,
        initialize={(i, j): pyo.value(primal.B[i, j])
                    for (i, j) in primal.linea}
    )
    
    #Variables
    model.flujo = pyo.Var(model.t*model.linea,domain=pyo.Reals)
    model.ang = pyo.Var(model.t*model.N,bounds=(-np.pi,np.pi)
                           ,domain=pyo.Reals)
    for t in model.t:
        model.ang[t,model.N.first()].fix(0)

    #Restricciones de linea
    model.flujoMax = pyo.Constraint(model.t*model.linea,rule=lambda m,t,i,j:
                                    m.flujo[t,i,j]<=m.Fmax[i,j])
    model.flujoMin = pyo.Constraint(model.t*model.linea,rule=lambda m,t,i,j:
                                    -m.Fmax[i,j]<=m.flujo[t,i,j])
    model.ecflujoDC = pyo.Constraint(model.t*model.linea,
                                     rule=lambda m,t,i,j:
                                     m.flujo[t,i,j]==m.B[i,j]*(m.ang[t,i]-m.ang[t,j]))
    
    #Funcion objetivo
    model.obj = pyo.Objective(rule=red_obj,sense=pyo.minimize)

    return model