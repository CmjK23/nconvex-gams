#%%
import numpy as np
import os
import pyomo.environ as pyo

#%%
#--------------------MODELO DE MAXIMIZACION DE BENEFICIO PARA TERMICAS------------------------#
#creacion de modelo abstracto
model = pyo.AbstractModel()

#SETS
model.t = pyo.Set(ordered=True)

def descartar_t_inicial(m):
    return list(m.t)[1:]
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
model.precio_n = pyo.Param(model.t,mutable=True)
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
model.ramp_up = pyo.Param(domain=pyo.NonNegativeReals, validate=lambda m, val_ramp_up: val_ramp_up > m.pmin)
model.ramp_dn = pyo.Param(domain=pyo.NonNegativeReals, validate=lambda m, val_ramp_dn: val_ramp_dn > m.pmin)

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
def costo_exp(model):
    return sum([model.y[t]*model.c_arr
               +model.u[t]*model.coef1
               +model.p[t]*model.coef2 for t in model.t])
model.costo = pyo.Expression(rule=costo_exp)

#Ingresos
def ingresos_exp(model):
    return sum([model.p[t]*model.precio_n[t] for t in model.t]) 
model.ingresos = pyo.Expression(rule=ingresos_exp)

#Variable objetivo: beneficio del generador
model.profit = pyo.Objective(expr=model.ingresos-model.costo,sense=pyo.maximize)

#RESTRICCIONES
def cond_lim_op1(m,t): #Condicion limite operativa
    return m.pmin*m.u[t] <= m.p[t]

def cond_lim_op2(m,t): #Condicion limite operativa
    return m.p[t] <= m.pmax*m.u[t]

model.bounds1 = pyo.Constraint(model.t,rule=cond_lim_op1)
model.bounds2 = pyo.Constraint(model.t,rule=cond_lim_op2)

#Restricciones de operacion por arranques/paradas
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

model.cond_op1 = pyo.Constraint(model.t,rule=cond_op1)
model.cond_op2 = pyo.Constraint(model.t,rule=cond_op2)
model.cond_op3 = pyo.Constraint(model.t,rule=cond_op3)
model.cond_op4 = pyo.Constraint(model.t,rule=cond_op4)
model.cond_op5 = pyo.Constraint(model.t,rule=cond_op5)

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

model.startup = pyo.Constraint(model.t,rule=startup)
model.shutdown = pyo.Constraint(model.t,rule=shutdown)

def ineq_arrd_max(m):
    return sum(m.y[t] for t in m.t) <= m.arrd_max

def ineq_pard_max(m):
    return sum(m.w[t] for t in m.t) <= m.pard_max

model.ineq_arrd_max = pyo.Constraint(rule=ineq_arrd_max)
model.ineq_pard_max = pyo.Constraint(rule=ineq_pard_max)

def ecComb(m): #Ecuacion de reserva de combustible
    return sum([m.p[t]*m.heat_rate + m.u[t]*m.no_load_consumption for t in m.t]) <= m.dComb

model.ecComb = pyo.Constraint(rule=ecComb)

#Tiempos minimos de encencido/apagado
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

model.ectminUp = pyo.Constraint(model.t_rest,rule=ectminUp)
model.ectminDn = pyo.Constraint(model.t_rest,rule=ectminDn)

# %%
#------------------------MODELO DE MAXIMIZACION DUAL-------------------------#
n_iteraciones = 80
alpha = 0.2
epsilon = 1e-8
termicas = ['G1','G2','G3','G4','G5']
workspace = 'Data'
Demanda = np.array([600,520,510,530,585,610,
                    640,620,670,720,820,900,
                    850,910,1000,980,1100,1050,
                    950,900,850,770,710,650]) 
opt = pyo.SolverFactory('cplex')
#opt = pyo.SolverFactory('glpk')
#%%
#Creando instancias para cada generador
modeloTermica = {}

os.chdir(workspace)

for ut in termicas:
    data = pyo.DataPortal()
    data.load(filename='common.dat')
    data.load(filename=f'{ut}.dat')
    
    modeloTermica[ut] = model.create_instance(data)

del data  
os.chdir('..')
# %%
#En cada iteracion del modelo maestro se tendran que restablecer los precios por hora en cada modelo
def restablecer_precios(precios,modelo):
    T = list(modelo.t)
    for i in range(len(T)):
        modelo.precio_n[T[i]] = precios[i]
# %% pruebas
Duales = []
precios = np.ones(24)*20

master = pyo.ConcreteModel()

master.t = pyo.Set(initialize=range(24),ordered=True)
master.precios = pyo.Var(master.t,bounds=(0,300),domain=pyo.NonNegativeReals)
master.ut = pyo.Set(initialize=termicas,ordered=True)
master.theta = pyo.Var(master.ut,domain=pyo.NonNegativeReals)
master.obj = pyo.Objective(expr=sum(master.precios[t]*Demanda[t] for t in master.t)
                           -sum(master.theta[ut] for ut in master.ut),
                           sense=pyo.maximize)
master.cuts = pyo.ConstraintList()


#----------------------------SUBPROBLEMA DE PROYECCION--------------------------#
projection = pyo.ConcreteModel()
projection.t = pyo.Set(initialize=range(24), ordered=True)
projection.precios = pyo.Var(projection.t, bounds=(0, 300))
# Parametro para el centro de proyeccion (precios de la iteracion anterior)
projection.precios_actuales = pyo.Param(projection.t, mutable=True, initialize=0)
projection.level_val = pyo.Param(mutable=True, initialize=0)

# Objetivo: Minimizar distancia cuadratica (QP)
projection.obj = pyo.Objective(expr=sum((projection.precios[t] - projection.precios_actuales[t])**2 for t in projection.t))

# En el Level Method, la proyeccion debe cumplir TODOS los cortes acumulados
# para estar dentro del "modelo" de la funcion dual por debajo del nivel
projection.cuts = pyo.ConstraintList()
#%%

#%% Bucle de optimización corregido
for k in range(n_iteraciones):
    # 1. Resolver Esclavos
    generacion_total_k = np.zeros(24)
    beneficio_acumulado = 0
    
    for ut in termicas:
        restablecer_precios(precios, modeloTermica[ut])
        opt.solve(modeloTermica[ut])
        
        pot_u = np.array([pyo.value(modeloTermica[ut].p[t]) for t in modeloTermica[ut].t])
        costo_u = pyo.value(modeloTermica[ut].costo)
        generacion_total_k += pot_u
        
        # Guardar corte en el Maestro (Multi-cut)
        # theta_g >= sum(p*pi) - costo
        master.cuts.add(sum(master.precios[t]*pot_u[t] for t in master.t) - costo_u <= master.theta[ut])
        
        # Guardar el mismo corte en la Proyección para definir el modelo f_k(pi)
        projection.cuts.add(sum(projection.precios[t]*pot_u[t] for t in projection.t) - costo_u <= master.theta[ut])

    # 2. Calcular Lk (Lower Bound actual)
    Lk = sum(Demanda * precios) - sum(pyo.value(m.profit) for m in modeloTermica.values())
    Duales.append(Lk)
    LB = max(Duales)

    # 3. Resolver Maestro para obtener Upper Bound (UB)
    opt.solve(master)
    UB = pyo.value(master.obj)

    # 4. Resolver Proyección
    # Actualizar nivel y punto de referencia
    current_level = alpha * UB + (1 - alpha) * LB
    for t in range(24):
        projection.precios_actuales[t] = precios[t]
    
    # El conjunto de nivel se define como: f_dual(pi) >= current_level
    # f_dual es sum(pi*D) - sum(theta_g)
    if k == 0:
        projection.level_const = pyo.Constraint(expr= sum(projection.precios[t]*Demanda[t] for t in projection.t) 
                                                - sum(master.theta[ut] for ut in termicas) >= current_level)
    else:
        # Actualizar el valor del nivel en la restricción existente
        projection.level_const.set_value(sum(projection.precios[t]*Demanda[t] for t in projection.t) 
                                         - sum(master.theta[ut] for ut in termicas) >= current_level)

    opt.solve(projection)
    
    # Actualizar precios para la siguiente iteración
    precios = np.array([pyo.value(projection.precios[t]) for t in projection.t])

    # 5. Criterio de parada
    gap = (UB - LB) / abs(UB) if abs(UB) > 1e-6 else 0
    print(f"Iter {k}: LB={LB:.2f}, UB={UB:.2f}, Gap={gap:.6f}")
    
    if gap <= epsilon:
        break
