#%%
import numpy as np
import os
import pyomo.environ as pyo

#%%
#------------------------------MODELO DE UNIT COMMITMENT BASADO EN YUPANA-----------------------#
model = pyo.AbstractModel()

#Set de tiempo
model.t = pyo.RangeSet(1,24)
def descartar_t_inicial(m):
    return list(m.t)[1:]
model.t_rest = pyo.Set(within=model.t,initialize=descartar_t_inicial)

#Set de unidades de generacion termicas
model.ut = pyo.Set(ordered=True)

model.u0 = pyo.Param(model.ut,domain=pyo.Binary)
model.p0 = pyo.Param(model.ut,domain=pyo.NonNegativeReals)

#Demanda del nodo
model.Demanda = pyo.Param(model.t,domain=pyo.NonNegativeReals,mutable=True)

#Parametros tecnicos de unidades termicas
model.pmin = pyo.Param(model.ut,domain=pyo.NonNegativeReals)
model.pmax = pyo.Param(model.ut,domain=pyo.NonNegativeReals,
                       validate=lambda m, val_pmax, ut: val_pmax > m.pmin[ut])

model.cComb = pyo.Param(model.ut,domain=pyo.NonNegativeReals)
model.heat_rate = pyo.Param(model.ut,domain=pyo.NonNegativeReals)
model.no_load_consumption = pyo.Param(model.ut,domain=pyo.NonNegativeReals)
model.dComb = pyo.Param(model.ut,domain=pyo.NonNegativeReals)

model.coef1 = pyo.Param(model.ut,domain=pyo.NonNegativeReals)
model.coef2 = pyo.Param(model.ut,domain=pyo.NonNegativeReals)
model.c_arr = pyo.Param(model.ut,domain=pyo.NonNegativeReals)

model.tminUp = pyo.Param(model.ut,domain=pyo.PositiveIntegers)
model.tminDn = pyo.Param(model.ut,domain=pyo.PositiveIntegers)

model.arrd_max = pyo.Param(model.ut,domain=pyo.PositiveIntegers)
model.pard_max = pyo.Param(model.ut,domain=pyo.PositiveIntegers)

model.ramp_n = pyo.Param(model.ut,domain=pyo.NonNegativeReals)
model.ramp_up = pyo.Param(model.ut,domain=pyo.NonNegativeReals, validate=lambda m, val_ramp_up,ut: val_ramp_up > m.pmin[ut])
model.ramp_dn = pyo.Param(model.ut,domain=pyo.NonNegativeReals, validate=lambda m, val_ramp_dn,ut: val_ramp_dn > m.pmin[ut])


#VARIABLES
#Potencia generada en cada bloque horario        
model.p = pyo.Var(model.ut,model.t,domain=pyo.NonNegativeReals)
#Estado de encendido/apagado
model.u = pyo.Var(model.ut,model.t,domain=pyo.Binary)
#Arranque
model.y = pyo.Var(model.ut,model.t,domain=pyo.Binary)
#Parada
model.w = pyo.Var(model.ut,model.t,domain=pyo.Binary)

#Funcion de costo
def costo_op(model):
    return sum(sum(model.y[ut,t]*model.c_arr[ut] +
               model.u[ut,t]*model.coef1[ut] +
               model.p[ut,t]*model.coef2[ut] for t in model.t) for ut in model.ut)
    
model.costo = pyo.Objective(rule=costo_op,sense=pyo.minimize)

#RESTRICCIONES

def cond_lim_op1(m,ut,t): #Condicion limite operativa
    return m.pmin[ut]*m.u[ut,t] <= m.p[ut,t]

def cond_lim_op2(m,ut,t): #Condicion limite operativa
    return m.p[ut,t] <= m.pmax[ut]*m.u[ut,t]

model.bounds1 = pyo.Constraint(model.ut,model.t,rule=cond_lim_op1)
model.bounds2 = pyo.Constraint(model.ut,model.t,rule=cond_lim_op2)

#Restricciones de operacion por arranques/paradas
def cond_op1(m,ut,t):
    if t==m.t.first():
        return m.u[ut,t]-m.u0[ut]== m.y[ut,t]- m.w[ut,t]
    else:
        t_prev = m.t.prev(t)
        return m.u[ut,t]-m.u[ut,t_prev]== m.y[ut,t]- m.w[ut,t]

def cond_op2(m,ut,t):
    return m.y[ut,t] <= m.u[ut,t]

def cond_op3(m,ut,t):
    if t==m.t.first():
        return m.y[ut,t] <= 1 - m.u0[ut]
    else:
        t_prev = m.t.prev(t)
        return m.y[ut,t] <= 1 - m.u[ut,t_prev]

def cond_op4(m,ut,t):
    if t==m.t.first():
        return m.w[ut,t] <= m.u0[ut]
    else:
        t_prev = m.t.prev(t)
        return m.w[ut,t] <= m.u[ut,t_prev]

def cond_op5(m,ut,t):
    return m.w[ut,t] <= 1 - m.u[ut,t]

model.cond_op1 = pyo.Constraint(model.ut,model.t,rule=cond_op1)
model.cond_op2 = pyo.Constraint(model.ut,model.t,rule=cond_op2)
model.cond_op3 = pyo.Constraint(model.ut,model.t,rule=cond_op3)
model.cond_op4 = pyo.Constraint(model.ut,model.t,rule=cond_op4)
model.cond_op5 = pyo.Constraint(model.ut,model.t,rule=cond_op5)

def startup(m,ut,t): #Restriccion de arranque
    if t==m.t.first():
        return m.p[ut,t] - m.p0[ut] <= m.ramp_n[ut]*m.u0[ut] + m.ramp_up[ut]*m.y[ut,t]
    else:
        t_prev = m.t.prev(t)
        return m.p[ut,t] - m.p[ut,t_prev] <= m.ramp_n[ut]*m.u[ut,t_prev] + m.ramp_up[ut]*m.y[ut,t]

def shutdown(m,ut,t): #Restriccion de desconexion
    if t==m.t.first():
        return m.p0[ut] - m.p[ut,t] <= m.ramp_n[ut]*m.u[ut,t] + m.ramp_dn[ut]*m.w[ut,t]
    else:
        t_prev = m.t.prev(t)
        return m.p[ut,t_prev] - m.p[ut,t] <= m.ramp_n[ut]*m.u[ut,t] + m.ramp_dn[ut]*m.w[ut,t]

model.startup = pyo.Constraint(model.ut,model.t,rule=startup)
model.shutdown = pyo.Constraint(model.ut,model.t,rule=shutdown)

def ineq_arrd_max(m,ut):
    return sum(m.y[ut,t] for t in m.t) <= m.arrd_max[ut]

def ineq_pard_max(m,ut):
    return sum(m.w[ut,t] for t in m.t) <= m.pard_max[ut]

model.ineq_arrd_max = pyo.Constraint(model.ut,rule=ineq_arrd_max)
model.ineq_pard_max = pyo.Constraint(model.ut,rule=ineq_pard_max)

def ecComb(m,ut): #Ecuacion de reserva de combustible
    return sum([m.p[ut,t]*m.heat_rate[ut] + m.u[ut,t]*m.no_load_consumption[ut] for t in m.t]) <= m.dComb[ut]

model.ecComb = pyo.Constraint(model.ut,rule=ecComb)

#Tiempos minimos de encencido/apagado
def ectminUp(m,ut,t):
    T = list(m.t)
    idx = T.index(t)
    start = max(0,idx-m.tminUp[ut]+1)
    return sum(m.y[ut,tau] for tau in T[start:idx+1]) <= m.u[ut,t]

def ectminDn(m,ut,t):
    T = list(m.t)
    idx = T.index(t)
    start = max(0,idx-m.tminDn[ut]+1)
    return sum(m.w[ut,tau] for tau in T[start:idx+1]) <= 1 - m.u[ut,t]

model.ectminUp = pyo.Constraint(model.ut,model.t_rest,rule=ectminUp)
model.ectminDn = pyo.Constraint(model.ut,model.t_rest,rule=ectminDn)

#Ecuacion de balance nodal
def ecBalNodal(m,t):
    return sum(m.p[ut,t] for ut in m.ut) == m.Demanda[t]

model.ecBalNodal = pyo.Constraint(model.t,rule=ecBalNodal)
# %%
#Se crea una instancia del modelo
#El archivo de prueba se llama test.dat
data = pyo.DataPortal()
data.load(filename='test.dat')
instance = model.create_instance(data)
# %%
#Dependiendo del solver MILP que tengan, se instancia un optimizador
opt = pyo.SolverFactory('glpk',tee=True)
#opt = pyo.SolverFactory('cplex',tee=True)
# %%
#Se resuelve el modelo de prueba
opt.solve(instance)
# %%
#Para obtener los costos marginales horarios como hace el modelo Quipu
#Se fijan las variables binarias previamente calculadas en el despacho
#De este modo el problema se convierte en un problema de LP y se pueden calcular los duales
for ut in instance.ut:
    for t in instance.t:
        instance.u[ut,t].fix(round(pyo.value(instance.u[ut,t])))
        instance.y[ut,t].fix(round(pyo.value(instance.y[ut,t])))
        instance.w[ut,t].fix(round(pyo.value(instance.w[ut,t])))
#%%
instance.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)
#%%
#Se vuelve a resolver el despacho con las variables binarias como constantes
opt.solve(instance)
#%%
#Impresion de costos marginales horarios
for t in instance.t:
    print(instance.dual[instance.ecBalNodal[t]])
# %%
