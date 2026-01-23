$title Despacho termico IEEE BUS 9 CON EL METODO DE NIVEL (MULTICORTES)

option limcol=10;
option limrow=10;
option decimals=8;

SET
    i             "generador"   /1*5/
    NB            "nodo"        /B001*B009/
    ENL           "Linea"       /E001*E009/
    k             "Iteraciones" /1*50/
    Cut(k)        "Cortes generados"
;

alias (NB,j);

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

b(NB,j) = sum(ENL$Fbus(ENL,NB,j),
         -Fdata(ENL,'X0')/(FData(ENL,'R0')*FData(ENL,'R0')+
         FData(ENL,'X0')*FData(ENL,'X0')));

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

*==============================================================================*
*--------------------PARAMETROS DE CONTROL DEL ALGORITMO-----------------------
*==============================================================================*
SCALAR
       alpha     "Parametro de nivel (step size parameter)(sensibilidad)" /0.2/
       Nivel     "Nivel objetivo para la proyeccion"
       LB        "Lower Bound (Mejor dual encontrado)" /-1e9/
       UB        "Upper Bound (Estimado por el maestro)" /1e9/
       Terminado "Bandera de parada" /0/
       Gap       "Brecha de convergencia" /0.1/
       epsilon   "Tolerancia de convergencia" /0.001/
       L_actual  "Valor real de la funcion dual";

PARAMETER
          pi_iter(NB)              "Precio candidato actual (pi_k)"
          pi_stab(NB)              "Centro de estabilidad (pi_stab)"
          pi_lo(NB)                "Limite inferior precios"
          pi_up(NB)                "Limite superior precios"
          hist_pi(NB, k)           "Precio nodal usado en la iteracion k"
          hist_profit_GT(GT, k)    "Beneficio termico en iteracion k"
          hist_p_GT(GT, k)         "Potencia termica en iteracion k"
          hist_profit_GH(GH, k)    "Beneficio hidro en iteracion k"
          hist_p_GH(GH, k)         "Potencia hidro en iteracion k";

* Almacenar evolucion del LB= Lower Bound, UB= Upper Bound y el gap *
PARAMETER hist_LB(k)  "Historial de Lower Bound"
          hist_UB(k)  "Historial de Upper Bound"
          hist_Gap(k) "Historial de Gap";

pi_iter(NB) = 30;
pi_stab(NB) = 30;
pi_lo(NB) = 0;
pi_up(NB) = 9000;

*inicializar cortes en vacio*
Cut(k) = no;

*==============================================================================*
*---------------------------------Generacion-----------------------------------
*==============================================================================*
VARIABLES
          z_slave_GT           "Beneficio de los Generadores Termico (Ecuacion 2a) "
          z_slave_GH           "Beneficio del Generadores Hidraulico (Ecuacion 2a) "
          z_slave_GT_ind(GT)   "Beneficio individual de cada GT"
          z_slave_GH_ind(GH)   "Beneficio individual de cada GH"

          p_GT(GT)             "Potencia Termica"
          u_GT(GT)             "Binaria de compromiso Termico"
          p_GH(GH)             "Potencia Hidraulica Total"

;

BINARY VARIABLE u_GT;
POSITIVE VARIABLE p_GT, p_GH;

EQUATIONS
         eq_obj_GT_total    "Objetivo total GT como suma de individuales"
         eq_obj_GH_total    "Objetivo total GH como suma de individuales"
         eq_obj_GT(GT)      "Max Beneficio Termico (Ecuacion 4b )"
         eq_min_GT(GT)      "Limite potencia termo (Pmin)"
         eq_max_GT(GT)      "Limite potencia termo (Pmax)"
         eq_obj_GH(GH)      "Max Beneficio Hidro (Ecuacion 4b)"
         eq_max_GH(GH)      "Limites potencia hidro (Pmax)"
         eq_min_GH(GH)      "Limites potencia hidro (Pmin)"
;

eq_obj_GT_total.. z_slave_GT =e= sum(GT, z_slave_GT_ind(GT));
eq_obj_GH_total.. z_slave_GH =e= sum(GH, z_slave_GH_ind(GH));

eq_obj_GT(GT)..
         z_slave_GT_ind(GT) =e=
          ( sum(NB$PtBus(GT,NB), pi_iter(NB)) * p_GT(GT) ) - (PtData(GT,'Costo1') * u_GT(GT)
          + PtData(GT,'CI1')*(p_GT(GT) - PtData(GT,'Pmin')*u_GT(GT))
          );

eq_min_GT(GT)..  p_GT(GT) =g= PtData(GT,'Pmin');

eq_max_GT(GT)..  p_GT(GT) =l= PtData(GT,'Pmax');

eq_obj_GH(GH)..
    z_slave_GH_ind(GH) =e= sum(NB$PhBus(GH,NB), pi_iter(NB)) * p_GH(GH)
                     - ( PhData(GH,'CI') * p_GH(GH) );
eq_max_GH(GH)..
    p_GH(GH) =l= PhData(GH,'Pmax');
eq_min_GH(GH)..
    P_GH(GH) =g= PhData(GH,'Pmin');
*==============================================================================*
*---------------------------------RED -----------------------------------
*==============================================================================*
VARIABLE
         z_slave_net        "Renta de Congestion Total de la Red"
         f_dc_linea(ENL)    "Flujo DC en linea"
         ang(NB)            "Angulo de voltaje"
         inj_fromto(NB)     "Inyeccion neta (Gen - Carga) vista por la red";

EQUATIONS
         eq_obj_net         " Congestion (Ecuacion 4c)"
         eq_flujo_dc(ENL)   "Ley de Kirchhoff DC (Ecuacion 10d)"
         eq_bal_net(NB)     "Balance Nodal Red"
;

eq_obj_net..    z_slave_net =e= sum(NB, pi_iter(NB) * inj_fromto(NB));

eq_flujo_dc(ENL)..  f_dc_linea(ENL) =e= sum((NB,j)$Fbus(ENL,NB,j), (ang(NB) - ang(j))/FData(ENL,'X0'));

eq_bal_net(NB)..   inj_fromto(NB) =e= sum((ENL,j)$Fbus(ENL,NB,j), f_dc_linea(ENL)) - sum((ENL,j)$Fbus(ENL,j,NB), f_dc_linea(ENL));

f_dc_linea.up(ENL) =  FData(ENL,'Pmax');
f_dc_linea.lo(ENL) = -FData(ENL,'Pmax');

*Nodo de referencia*
ang.fx('B001') = 0;
ang.lo(NB)=-PI;
ang.up(NB)=PI;

Model Slave_GT  /eq_obj_GT_total, eq_obj_GT,eq_min_GT,eq_max_GT/;
Model Slave_GH  /eq_obj_GH_total, eq_obj_GH, eq_min_GH/;
Model Slave_Net /eq_obj_net, eq_flujo_dc, eq_bal_net/;

Slave_GT.optcr = 0;

* ==============================================================================
*                        PROBLEMA MAESTRO Y PROYECCIÓN
* ==============================================================================


VARIABLES
            Z_dual             "Valor aproximado de la funcion dual (L_hat)"
            theta_GT(GT)       "Variable auxiliar: Beneficio aproximado Termico"
            theta_GH(GH)       "Variable auxiliar: Beneficio aproximado Hidro"
            theta_net          "Variable auxiliar: Renta de Congestion aproximada"
            pi_new(NB)         "Nuevo precio candidato (Variable de decision)"
            dist               "Distancia cuadratica al centro de estabilidad"
;

EQUATIONS
            eq_obj_master      "Maximizar funcion dual aproximada (Ecuacion 6 / 12a)"
            eq_cut_GT(GT, k)   "Multicortes para Termicos (Ecuacion 12b)"
            eq_cut_GH(GH, k)   "Multicortes para Hidros (Ecuacion 12b)"
            eq_cut_net(k)      "Cortes para la Red (Ecuacion 12d)"
            eq_min_dist        "Minimizar distancia (Ecuacion 9 Objetivo)"
            eq_nivel_rst       "Restriccion de Nivel (Ecuacion 9 Restriccion)"
;

eq_obj_master..    Z_dual =e= sum(NB, pi_new(NB) * demanda(NB,'Pc'))
             - sum(GT, theta_GT(GT))
             - sum(GH, theta_GH(GH))
             - theta_net;

eq_cut_GT(GT, Cut)..    theta_GT(GT) =g= hist_profit_GT(GT, Cut)
                   + sum(NB$PtBus(GT,NB),
                         hist_p_GT(GT, Cut) * (pi_new(NB) - hist_pi(NB, Cut))
                     );

eq_cut_GH(GH, Cut)..    theta_GH(GH) =g= hist_profit_GH(GH, Cut)
                   + sum(NB$PhBus(GH,NB),
                         hist_p_GH(GH, Cut) * (pi_new(NB) - hist_pi(NB, Cut))
                     );

eq_cut_net(Cut)..
    theta_net =g= 0;

eq_min_dist..
    dist =e= sum(NB, sqr(pi_new(NB) - pi_stab(NB)));

eq_nivel_rst..
    Z_dual =g= Nivel;

Model Maestro    /eq_obj_master, eq_cut_GT, eq_cut_GH, eq_cut_net/;
Model Proyeccion /eq_obj_master, eq_cut_GT, eq_cut_GH, eq_cut_net, eq_min_dist, eq_nivel_rst/;

Maestro.optcr = 0;
Proyeccion.optcr = 0;

pi_new.lo(NB) = 0;
pi_new.up(NB) = 2000;
Scalar k_count "Numero total de iteraciones" /0/;

Loop(k$(ord(k) <= 50 and Terminado = 0),

    k_count=k_count+1;

    Solve Slave_GT using mip maximizing z_slave_GT;

    Solve Slave_GH using lp maximizing z_slave_GH;

    Solve Slave_Net using lp minimizing z_slave_net;

    Cut(k) = yes;

    hist_pi(NB, k) = pi_iter(NB);

    hist_profit_GT(GT, k) = z_slave_GT_ind.l(GT);
    hist_p_GT(GT, k) = p_GT.l(GT);

    hist_profit_GH(GH, k) = z_slave_GH_ind.l(GH);
    hist_p_GH(GH, k) = p_GH.l(GH);

    L_actual = sum(NB, pi_iter(NB) * demanda(NB,'Pc'))
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

    Gap = (UB - LB) / (abs(LB));
    hist_LB(k) = LB;
    hist_UB(k) = UB;
    hist_Gap(k) = Gap;

    if(Gap < epsilon,
        Terminado = 1;
    );

    if(Terminado = 0,
        Nivel = alpha * UB + (1 - alpha) * LB;
        Solve Proyeccion using qcp minimizing dist;
        pi_iter(NB) = pi_new.l(NB);
    );

);

Display "Numero de iteraciones", k_count;
Display "Historial del gap ", hist_Gap;
Display "Historial del Lowwer Bound ", hist_LB;
Display "Historial del Upper  Bound ", hist_UB;
Display hist_profit_GT,hist_p_GT,hist_profit_GH,hist_p_GH;

