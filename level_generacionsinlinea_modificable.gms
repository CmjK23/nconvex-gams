$title Despacho termico IEEE BUS 9 CON EL METODO DE NIVEL (MULTICORTES)

option limcol=10;
option limrow=10;
option decimals=4;

SET
    i             "generador"   /1*5/
    NB            "nodo"        /B001*B009/
    ENL           "Linea"       /E001*E009/
    k             "Iteraciones" /1*50/
    Cut(k)        "Cortes generados"
;

alias (NB,j);

VARIABLE
ang(NB)
p(i)
;

positive variable
prac(NB);
scalar
crac /100000000/;

*Parque Generador
*********************************************
* Termica
*********************************************
SET  GT  /
        T001 G2Steam 18       TG-1
        T002 G3Steam 13.8     TG-2
/;

SET  PtBus(GT,NB) /
        T001.B002
        T002.B003
/;

TABLE PtData(GT,*)
                Pgen      Pmin      Pmax    Costo1    CI1       CI2     Pmax1     Pmax2   Forzada   Sist  Calif
        T001  163.000   50.000   192.000   8000.00   50.00     50.00    192.000   192.000     no      1     101
        T002   85.000   30.000   128.000   9000.00   60.00     60.00    128.000   128.000     no      1     102
;
*********************************************
* Hidraulica
*********************************************
SET  GH  /
        H001 G1Hydro 16.5     G-1
/;

SET  PhBus(GH,NB) /
        H001.B001
/;

TABLE PhData(GH,*)
                Pgen      Pmin      Pmax      CI        Forzada   Sist    Cmgh
        H001   71.600   20.000   247.500   2.00         no        1      yes
;

TABLE  FData(ENL,*)
                  R0       X0         G0      Pmax     Cong  Sist
        E001    0.000000  0.0576   0.00000  0250.000    no     1
        E002    0.000000  0.0625   0.00000  0200.000    no     1
        E003    0.000000  0.0586   0.00000  0150.000    no     1
        E004    0.010000  0.0850   0.00000  0300.000    no     1
        E005    0.017000  0.0920   0.00000  0300.000    no     1
        E006    0.032000  0.1610   0.00000  0300.000    no     1
        E007    0.039000  0.1700   0.00000  0300.000    no     1
        E008    0.008500  0.0720   0.00000  0300.000    no     1
        E009    0.011900  0.1008   0.00000  0300.000    no     1
;

*Conexion de la linea
SET Fbus(ENL,NB,NB)
/    E001. B004.B001
     E002. B004.B005
     E003. B006.B005
     E004. B003.B006
     E005. B006.B007
     E006. B007.B008
     E007. B008.B002
     E008. B008.B009
     E009. B009.B004  /;

PARAMETER b(NB,NB);

b(NB,j) = sum(ENL$Fbus(ENL,NB,j),-Fdata(ENL,'X0')/(FData(ENL,'R0')*FData(ENL,'R0')+FData(ENL,'X0')*FData(ENL,'X0')));

DISPLAY b;

TABLE demanda(NB,*)
         Pc
B002     0.8000
B005     0.7068
B007     0.7854
B009     0.9817
;
DISPLAY demanda;
PARAMETER Loss(ENL);
Loss(ENL)=0;

ang.lo(NB)=-PI;
ang.up(NB)=PI;
*--------------------PARAMETROS DE CONTROL DEL ALGORITMO-----------------------

SCALAR alpha     "Parametro de nivel (step size parameter)" /0.5/;
SCALAR Nivel     "Nivel objetivo para la proyeccion";
SCALAR LB        "Lower Bound (Mejor dual encontrado)" /-1e9/;
SCALAR UB        "Upper Bound (Estimado por el maestro)" /1e9/;
SCALAR Terminado "Bandera de parada" /0/;
SCALAR Gap       "Brecha de convergencia" /1/;

*==============================================================================*
*                           MEMORIA DEL ESCLAVO
*==============================================================================*
* -----------1. MEMORIA PARAMETROS DE CONTROL-----------------------------------

PARAMETER hist_LB(k)     "Historial de Lower Bound";
PARAMETER hist_UB(k)     "Historial de Upper Bound";
PARAMETER hist_Gap(k)    "Historial de Gap";
PARAMETER hist_Nivel(k)  "Historial de Nivel";

* -----------1. MEMORIA PARA GENERADORES TÉRMICOS (GT) -------------------------
* Corresponde a la Ecuación (12b) del paper para el set G_termico
* Almacenamos el Beneficio (H) y la Potencia (p) de cada iteración k

PARAMETER hist_profit_GT(GT, k)    "Termino constante del corte Termico (Beneficio)";
PARAMETER hist_p_GT(GT, k)         "Gradiente del corte Termico (Potencia)";


* ----------- 2. MEMORIA PARA GENERADORES HIDRÁULICOS (GH) ---------------------
* Corresponde a la Ecuación (12b) del paper para el set G_hidro

PARAMETER hist_profit_GH(GH, k)    "Termino constante del corte Hidro (Beneficio)";
PARAMETER hist_p_GH(GH, k)         "Gradiente del corte Hidro (Potencia)";


* --- 3. MEMORIA PARA LA RED (ENLACES) ---
* Corresponde a la Ecuación (12d) del paper.
* Almacenamos la Renta de Congestion (H_net) y la Inyeccion Neta (S)

PARAMETER hist_profit_red(k)       "Termino constante corte Red (Renta Congestion)";
PARAMETER hist_inj_net(NB, k)      "Gradiente corte Red (Inyeccion Neta)";


* --- 4. HISTORIAL DE PRECIOS ---
* Necesario para calcular (pi - pi_k) en las ecuaciones de corte

PARAMETER hist_pi(NB, k)           "Precio nodal usado en la iteracion k";

* --- 5. VARIABLES DE ESTADO DEL ALGORITMO ---
* Precios actuales y centro de estabilidad (Eq. 9 del paper)

PARAMETER pi_iter(NB)              "Precio candidato actual (pi_k)";
PARAMETER pi_stab(NB)              "Centro de estabilidad (pi_stab)";

* Inicialización
pi_iter(NB) = 30;
pi_stab(NB) = 30;
Cut(k) = no;
*==============================================================================*
* Variables de los problema esclavos
*==============================================================================*
VARIABLES
    z_slave_GT           "Beneficio de los Generadores Termico "
    z_slave_GH           "Beneficio del Generadores Hidraulico "
    z_slave_net          "Renta de Congestion Total de la Red"
    z_slave_GT_ind(GT)   "Beneficio individual de cada GT"
    z_slave_GH_ind(GH)   "Beneficio individual de cada GH"

    p_GT(GT)        "Potencia Termica Total"
    p1_GT(GT)       "Potencia Termica Tramo 1"
    p2_GT(GT)       "Potencia Termica Tramo 2"
    u_GT(GT)        "Binaria de compromiso Termico"

    p_GH(GH)        "Potencia Hidraulica Total"

    f_slave(ENL)    "Flujo DC en linea"
    theta_slave(NB) "Angulo de voltaje"
    inj_slave(NB)   "Inyeccion neta (Gen - Carga) vista por la red"
;

BINARY VARIABLE  u_GT;
POSITIVE VARIABLE p1_GT, p2_GT, p_GH;

*==============================================================================*
* Ecuaciones de los esclavos del generador (MIP/LP)
*==============================================================================*
EQUATIONS

    eq_obj_GT_total    "Objetivo total GT como suma de individuales"
    eq_obj_GH_total    "Objetivo total GH como suma de individuales"
    eq_obj_GT(GT)      "Max Beneficio Termico (Eq. 4b Paper)"
    eq_suma_GT(GT)     "Suma de tramos potencia"
    eq_lim_p1(GT)      "Limite tramo 1"
    eq_lim_p2(GT)      "Limite tramo 2"
    eq_obj_GH(GH)      "Max Beneficio Hidro (Eq. 4b Paper)"
    eq_lim_GH(GH)      "Limites potencia hidro"
;

eq_obj_GT_total.. z_slave_GT =e= sum(GT, z_slave_GT_ind(GT));
eq_obj_GH_total.. z_slave_GH =e= sum(GH, z_slave_GH_ind(GH));


* --- 1. ESCLAVO TÉRMICO (Eq. 4b: Max pi*p - Costo(p)) ---
* El generador ve el precio pi_iter y decide cuanto producir para ganar dinero.
eq_obj_GT(GT)..
    z_slave_GT_ind(GT) =e= sum(NB$PtBus(GT,NB), pi_iter(NB)) * p_GT(GT)
                     - ( PtData(GT,'Costo1')*u_GT(GT)
                       + PtData(GT,'CI1')*p1_GT(GT)
                       + PtData(GT,'CI2')*p2_GT(GT) );

* Definición de tramos (Piecewise Linear)
eq_suma_GT(GT).. p_GT(GT) =e= PtData(GT,'Pmin')*u_GT(GT) + p1_GT(GT) + p2_GT(GT);
eq_lim_p1(GT)..  p1_GT(GT) =l= (PtData(GT,'Pmax1') - PtData(GT,'Pmin')) * u_GT(GT);
eq_lim_p2(GT)..  p2_GT(GT) =l= (PtData(GT,'Pmax2') - PtData(GT,'Pmax1')) * u_GT(GT);


* --- 2. ESCLAVO HIDRÁULICO (Eq. 4b Simplificada) ---
* Asumimos costo de arranque 0 y un solo tramo con costo CI (Valor del agua)
eq_obj_GH(GH)..
    z_slave_GH_ind(GH) =e= sum(NB$PhBus(GH,NB), pi_iter(NB)) * p_GH(GH)
                     - ( PhData(GH,'CI') * p_GH(GH) );

eq_lim_GH(GH)..
    p_GH(GH) =l= PhData(GH,'Pmax');
* Nota: Si tienes mínimo técnico hidro, agrega: p_GH =g= PhData(GH,'Pmin');


*==============================================================================*
* Ecuaciones de los esclavos de la red (LP)
*==============================================================================*
EQUATIONS
    eq_obj_net         "Max Renta Congestion (Eq. 4c Paper)"
    eq_flujo_dc(ENL)   "Ley de Kirchhoff DC (Eq. 10d Paper)"
    eq_bal_net(NB)     "Balance Nodal Red"
;

* --- 3. ESCLAVO RED (Eq. 4c y Eq. 10) ---
* La red maximiza el cobro por transportar energía: Sum(Precio * InyeccionNeta)
eq_obj_net..
    z_slave_net =e= sum(NB, pi_iter(NB) * inj_slave(NB));

* Flujo DC: F = B * (Angulo_i - Angulo_j)
eq_flujo_dc(ENL)..
    f_slave(ENL) =e= sum((NB,j)$Fbus(ENL,NB,j), (theta_slave(NB) - theta_slave(j))/FData(ENL,'X0'));

* Balance Nodal: Inyeccion_neta = Flujos que salen - Flujos que entran
* Esta variable 'inj_slave' es el gradiente que usaremos en el corte maestro
eq_bal_net(NB)..
    inj_slave(NB) =e= sum((ENL,j)$Fbus(ENL,NB,j), f_slave(ENL))
                    - sum((ENL,j)$Fbus(ENL,j,NB), f_slave(ENL));

* Límites de transmisión (Eq. 10b/10c Paper)
f_slave.up(ENL) =  FData(ENL,'Pmax');
f_slave.lo(ENL) = -FData(ENL,'Pmax');

* Nodo Slack (Referencia angular)
theta_slave.fx('B001') = 0;


* --- DEFINICIÓN DE MODELOS ---
MODEL  Slave_GT  /eq_obj_GT_total, eq_obj_GT, eq_suma_GT, eq_lim_p1, eq_lim_p2/;
MODEL Slave_GH  /eq_obj_GH_total, eq_obj_GH, eq_lim_GH/;
MODEL Slave_Net /eq_obj_net, eq_flujo_dc, eq_bal_net/;

* Opciones para resolver
Slave_GT.optcr = 0;
* ==============================================================================
* PASO 4: PROBLEMA MAESTRO Y PROYECCIÓN
* ==============================================================================

VARIABLES

    Z_dual             "Valor aproximado de la funcion dual (L_hat)"
    theta_GT(GT)       "Variable auxiliar: Beneficio aproximado Termico"
    theta_GH(GH)       "Variable auxiliar: Beneficio aproximado Hidro"
    theta_net          "Variable auxiliar: Renta de Congestion aproximada"

    pi_new(NB)         "Nuevo precio candidato (Variable de decision)"
    mu(ENL)        >= 0 "dual limite superior flujo"
    nu(ENL)        >= 0 "dual limite inferior flujo"
    lambda(NB)          "dual balance flujo"

    dist               "Distancia cuadratica al centro de estabilidad"
;

EQUATIONS

    eq_obj_master      "Maximizar funcion dual aproximada (Eq. 6 / 12a)"
    eq_cut_GT(GT, k)   "Multicortes para Termicos (Eq. 12b)"
    eq_cut_GH(GH, k)   "Multicortes para Hidros (Eq. 12b)"
    eq_cut_net(k)      "Cortes para la Red (Eq. 12d)"
    eq_network_dual1(ENL) "Restricción dual de flujo DC"
    eq_network_dual2(NB)  "Restricción dual de balance nodal"

    eq_min_dist        "Minimizar distancia (Eq. 9 Objetivo)"
    eq_nivel_rst       "Restriccion de Nivel (Eq. 9 Restriccion)"
;

* ------------------------------------------------------------------------------
* 1. FUNCIÓN OBJETIVO DEL MAESTRO (Eq. 6 y 12a del Paper)
* ------------------------------------------------------------------------------
* L(pi) = Ingreso_Demanda - Beneficio_Gen - Renta_Red
* El objetivo es encontrar el 'pi' que MAXIMIZA esta función cóncava.
* Usamos variables auxiliares (theta) para representar los beneficios convexos.

eq_obj_master..
    Z_dual =e= sum(NB, pi_new(NB) * demanda(NB,'Pc'))
             + sum(ENL, nu_slave(ENL)*FData(ENL,'Pmax')
                    - mu_slave(ENL)*FData(ENL,'Pmax'))
             - sum(GT, theta_GT(GT))
             - sum(GH, theta_GH(GH))
             - theta_net;

* ------------------------------------------------------------------------------
* 2. MULTICORTES PARA GENERADORES TÉRMICOS (Eq. 12b del Paper)
* ------------------------------------------------------------------------------
* theta_g >= Beneficio_Hist + Gradiente * (pi_new - pi_hist)
* El Gradiente es la Potencia Generada (hist_p_GT).
* Se genera un corte para CADA generador en CADA iteración 'Cut'.

eq_cut_GT(GT, Cut)..  theta_GT(GT) =g= hist_profit_GT(GT, Cut) + sum(NB$PtBus(GT,NB),
                         hist_p_GT(GT, Cut) * (pi_new(NB) - hist_pi(NB, Cut))
                     );

* ------------------------------------------------------------------------------
* 3. MULTICORTES PARA GENERADORES HIDRÁULICOS (Eq. 12b del Paper)
* ------------------------------------------------------------------------------
* Misma lógica Eq. 12b, aplicada al conjunto GH.

eq_cut_GH(GH, Cut)..
    theta_GH(GH) =g= hist_profit_GH(GH, Cut)
                   + sum(NB$PhBus(GH,NB),
                         hist_p_GH(GH, Cut) * (pi_new(NB) - hist_pi(NB, Cut))
                     );


eq_network_dual1(ENL)..
    pi_new(FData(ENL,'from_node')) - pi_new(FData(ENL,'to_node'))
    + mu_slave(ENL) - nu_slave(ENL) + lambda_slave(FData(ENL,'from_node'))
    - lambda_slave(FData(ENL,'to_node')) =e= 0;

eq_network_dual2(NB)..
    sum(ENL$(FData(ENL,'to_node') = NB), lambda_slave(NB)/FData(ENL,'X0'))
    - sum(ENL$(FData(ENL,'from_node') = NB), lambda_slave(NB)/FData(ENL,'X0')) =e= 0;


* ------------------------------------------------------------------------------
* 5. PROBLEMA DE PROYECCIÓN (Eq. 9 del Paper)
* ------------------------------------------------------------------------------
* Objetivo: Minimizar la distancia euclidiana al cuadrado respecto al centro pi_stab.
* || pi - pi_stab ||^2

eq_min_dist..    dist =e= sum(NB, sqr(pi_new(NB) - pi_stab(NB)));

* Restricción de Nivel:
* El valor de la función dual (Z_dual) debe ser al menos el Nivel calculado.
* Z_dual >= (alpha*UB + (1-alpha)*LB)

eq_nivel_rst..    Z_dual =g= Nivel;


* ==============================================================================
* DEFINICIÓN DE MODELOS
* ==============================================================================

* Modelo Maestro (LP): Solo busca maximizar el techo construido por los cortes.
MODEL  Maestro /eq_obj_master, eq_cut_GT, eq_cut_GH, eq_cut_net,
                eq_network_dual1, eq_network_dual2/;

* Modelo Proyección (QCP): Busca estabilidad sujeto a mejorar el objetivo.
MODEL Proyeccion /Maestro, eq_min_dist, eq_nivel_rst/;

* Configuración de Solvers
Maestro.optcr = 0;
Proyeccion.optcr = 0;

* Límites para evitar precios infinitos en las primeras iteraciones
pi_new.lo(NB) = 0;
pi_new.up(NB) = 2000;
* ==============================================================================
* PASO 5: BUCLE ITERATIVO (ALGORITMO DE NIVEL)
* ==============================================================================

SCALAR L_actual "Valor real de la funcion dual en la iteracion k";
SCALAR k_count "Contador de iteraciones ejecutadas" /0/;

Loop(k$(ord(k) <= 50 and Terminado = 0),

    k_count=k_count+1;
    Solve Slave_GT using mip maximizing z_slave_GT;
    Solve Slave_GH using lp maximizing z_slave_GH;


    Cut(k) = yes;
    hist_pi(NB, k) = pi_iter(NB);

    hist_p_GT(GT, k) = p_GT.l(GT);
    hist_p_GH(GH, k) = p_GH.l(GH);

    hist_profit_GT(GT, k) = z_slave_GT_ind.l(GT);
    hist_profit_GH(GH, k) = z_slave_GH_ind.l(GH);

    hist_profit_red(k) = z_slave_net.l;
    hist_inj_net(NB, k) = inj_slave.l(NB);

    L_actual = sum(NB, pi_iter(NB) * Demanda(NB,'Pc'))
             - sum(GT, z_slave_GT_ind.l(GT))
             - sum(GH, z_slave_GH_ind.l(GH))
             - z_slave_net.l;
    if(L_actual > LB,
        LB = L_actual;
        pi_stab(NB) = pi_iter(NB);
    );

    Display "Iteracion:", k_count, "LB:", LB, "UB:", UB, "Gap:", gap;

    Solve Maestro using lp maximizing Z_dual;

    UB = Z_dual.l;
    gap = (UB - LB) / (abs(LB) + 1e-5);

    if(gap < 0.001,
        Terminado = 1;
        Display "CONVERGENCIA ALCANZADA EXITOSAMENTE";
    );

    if(Terminado = 0,
        Nivel = alpha * UB + (1 - alpha) * LB;
        Solve Proyeccion using qcp minimizing dist;
        pi_iter(NB) = pi_new.l(NB);
    );
 );
* ==============================================================================
*                                 REPORTE FINAL
* ==============================================================================
Display "Precios Finales (CHP):", pi_iter;
Display "iteraciones:",k_count;
Display "Historial de potencia térmica generada:", hist_p_GT, "en la iteracion",k_count;
Display hist_profit_GH,hist_profit_GT,hist_pi;

