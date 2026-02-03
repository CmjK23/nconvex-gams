$title Despacho termico IEEE BUS 39 CON EL METODO DE NIVEL (MULTICORTES)
$if not set DATA $abort "Falta definir DATA"
$if not set PRICE $abort "Falta definir PRICE"
$if not set GEN $abort "Falta definir GEN"
$if not set HIST $abort "Falta definir HIST"
$include %DATA%

option limcol=10;
option limrow=10;
option decimals=8;

SET
     K         "Iteraciones" /1*50/
     cut(k)    "cortes generados";
alias (NB,j);
positive variable
prac(NB);
scalar
crac /100000000/;

PARAMETER b(NB,NB);

b(NB,j) = sum(ENL$Fbus(ENL,NB,j),
         -Fdata(ENL,'X0')/(FData(ENL,'R0')*FData(ENL,'R0')+
         FData(ENL,'X0')*FData(ENL,'X0')));

DISPLAY b;
DISPLAY demanda;
PARAMETER Loss(ENL);
Loss(ENL)=0;

*==============================================================================*
*--------------------PARAMETROS DE CONTROL DEL ALGORITMO------------------------
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
*-----------------------------------RED-----------------------------------------
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

eq_bal_net(NB)..       inj_fromto(NB) =e= (sum(GT$PtBus(GT,NB), p_GT(GT))
                       + sum(GH$PhBus(GH,NB), p_GH(GH))
                       - demanda(NB,'Pc'))
                     + (sum(ENL$Fbus(ENL,NB,NB), f_dc_linea(ENL))
                       - sum(ENL$Fbus(ENL,NB,NB), f_dc_linea(ENL)));

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
* ------------------PROBLEMA MAESTRO Y PROYECCIÓN-------------------------------
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
    theta_net =g=0;

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

    Solve Slave_Net using lp maximizing z_slave_net;

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

* Configuración de archivo para Generacion
FILE f_gen /'%GEN%'/;
* .pc = 5 indica formato separado por comas (CSV estándar)
f_gen.pc = 5;
* .pw = 10000 indica el ancho de pagina (para evitar que corte lineas)
f_gen.pw = 10000;

PUT f_gen;
* Cabecera del CSV
PUT 'Tipo,Generador,Nodo,Potencia_MW,Estado_Commitment,Precio_Nodal_USD' /;

* 1. Exportar Generadores Térmicos
LOOP((GT, NB)$PtBus(GT,NB),
    PUT 'Termica', GT.tl, NB.tl, p_GT.l(GT):0:4, u_GT.l(GT):0:0, pi_iter(NB):0:2 /;
);

* 2. Exportar Generadores Hidráulicos
LOOP((GH, NB)$PhBus(GH,NB),
    PUT 'Hidraulica', GH.tl, NB.tl, p_GH.l(GH):0:4, '1', pi_iter(NB):0:2 /;
);


* ------------------------------------------------------------------------------

* Configuración de archivo para Precios Nodales
FILE f_price /'%PRICE%'/;
f_price.pc = 5;
f_price.pw = 10000;

PUT f_price;
PUT 'Nodo,Precio_Final' /;

LOOP(NB,
    PUT NB.tl, pi_iter(NB):0:4 /;
);


* ------------------------------------------------------------------------------

*Configuración de archivo para Historial de Convergencia
FILE f_conv /'%HIST%'/;
f_conv.pc = 5;
f_conv.pw = 10000;

PUT f_conv;
PUT 'Iteracion,Lower_Bound,Upper_Bound,Gap' /;

*Solo iteramos hasta k_count (las iteraciones que realmente ocurrieron)
LOOP(k$(ord(k) <= k_count),
    PUT k.tl, hist_LB(k):0:4, hist_UB(k):0:4, hist_Gap(k):0:6 /;
);
