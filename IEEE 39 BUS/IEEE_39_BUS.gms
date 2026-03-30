$title Despacho IEEE BUS 39 - Level Method Multicortes (VERSION FINAL)
*==============================================================================*
* Paper base: Stevens & Papavasiliou (IEEE Trans. Power Syst. 2022)            *
*==============================================================================*
* CORRECCIONES APLICADAS EN ESTA VERSION:                                      *
*                                                                               *
* [F1]  alias(NB,j) ELIMINADO — j ya existe en el DAT como alias interno       *
*        Se usa NBB como alias de NB en eq_flujo_dc y eq_bal_net               *
* [F2]  PI renombrado a PIVAL — PI es palabra reservada en algunos contextos   *
*        de GAMS y puede generar el error 141 "symbol not initialized"         *
* [F3]  Acoplamiento binario en GT: Pmin*u y Pmax*u (UC estandar)             *
* [F4]  p_GT y p_GH fijados (.fx) antes de Slave_Net, restaurados despues     *
* [F5]  Supergradiente de red = inj_fromto.l (teorema del envelope)           *
* [F6]  Signo correcto L_actual: +z_slave_net                                 *
* [F7]  Signo correcto eq_cut_net: gradiente negativo                         *
* [F8]  activo(k) como filtro $(activo(kk)) en ecuaciones de corte            *
* [F9]  Gap = (UB-LB)/(abs(UB)+1e-6) segun Eq.(8) del paper                  *
* [F10] Todos los abort removidos — reemplazados por Display de alerta        *
*        Los abort cortaban ejecucion antes de llegar a los CSV               *
* [F11] PUTCLOSE en todos los archivos CSV                                    *
* [F12] pi_iter inicializado en 20 USD/MWh segun paper Seccion IV             *
* [F13] Display fuera del loop para no ralentizar ejecucion                   *
* [F14] Encogimiento de caja Q en iter 10, 20, 30 (paper Sec. III-C)         *
*==============================================================================*

$INCLUDE /39_FINAL_bus_QUIPU_input.DAT

option limcol = 10;
option limrow = 10;
option decimals = 8;
option solprint = on;

*==============================================================================*
*--------------------------------- CONJUNTOS ----------------------------------*
*==============================================================================*
SET K "Iteraciones del Level Method" /1*50/;

alias(k, kk);

* [F1] Se usa NBB en lugar de j para evitar conflicto con alias del DAT
alias(NB, NBB);

*==============================================================================*
*--------------------- PARAMETROS DE CONTROL DEL ALGORITMO -------------------*
*==============================================================================*

* [F2] PIVAL en lugar de PI para evitar conflicto con nombre reservado
SCALAR
       PIVAL     "Numero pi"                              /3.14159265/
       alpha     "Parametro de nivel"                    /0.2/
       Nivel     "Nivel objetivo para la proyeccion"     /0/
       LB        "Lower Bound: mejor valor dual"         /-1e9/
       UB        "Upper Bound: estimado por el maestro"  /1e9/
       Terminado "Bandera de parada del loop"            /0/
       Gap       "Brecha relativa de convergencia"       /1/
       epsilon   "Tolerancia de convergencia"            /0.001/
       L_actual  "Valor de L(pi) en iteracion actual"    /0/
       total_dem "Demanda total del sistema"             /0/
       diferencia "Balance generacion menos demanda"     /0/
       flag_err  "Flag de error en algun solver"         /0/;

PARAMETER
       activo(k)              "1 si el corte k fue generado"
       pi_iter(NB)            "Precio candidato actual pi_k"
       pi_stab(NB)            "Centro de estabilidad"
       hist_pi(NB,k)          "Precio nodal en iteracion k"
       hist_profit_GT(GT,k)   "Beneficio termico en iteracion k"
       hist_p_GT(GT,k)        "Potencia termica en iteracion k"
       hist_profit_GH(GH,k)   "Beneficio hidro en iteracion k"
       hist_p_GH(GH,k)        "Potencia hidro en iteracion k"
       hist_dual_net(NB,k)    "Supergradiente de red en iteracion k"
       hist_z_net(k)          "Valor N(pi) en iteracion k"
       hist_LB(k)             "Historial Lower Bound"
       hist_UB(k)             "Historial Upper Bound"
       hist_Gap(k)            "Historial Gap relativo"
       hist_gen_GT(GT,k)      "Potencia termica por iteracion"
       hist_gen_GH(GH,k)      "Potencia hidro por iteracion"
       hist_u_GT(GT,k)        "Commitment termico por iteracion";

* Inicializaciones
activo(k)   = 0;
* [F12] Inicializar en 20 USD/MWh segun paper Seccion IV
pi_iter(NB) = 20;
pi_stab(NB) = 20;

*==============================================================================*
*------------------------------ SLAVES: GENERACION ---------------------------*
*==============================================================================*
VARIABLES
       z_slave_GT           "Beneficio total termico"
       z_slave_GH           "Beneficio total hidro"
       z_slave_GT_ind(GT)   "Beneficio individual GT"
       z_slave_GH_ind(GH)   "Beneficio individual GH";

BINARY   VARIABLE u_GT(GT)  "Compromiso binario GT";
POSITIVE VARIABLE p_GT(GT)  "Potencia termica [MW]";
POSITIVE VARIABLE p_GH(GH)  "Potencia hidro [MW]";

EQUATIONS
       eq_obj_GT_total
       eq_obj_GH_total
       eq_obj_GT(GT)
       eq_min_GT(GT)
       eq_max_GT(GT)
       eq_obj_GH(GH)
       eq_min_GH(GH)
       eq_max_GH(GH);

eq_obj_GT_total..
    z_slave_GT =e= sum(GT, z_slave_GT_ind(GT));

eq_obj_GH_total..
    z_slave_GH =e= sum(GH, z_slave_GH_ind(GH));

eq_obj_GT(GT)..
    z_slave_GT_ind(GT) =e=
          sum(NB$PtBus(GT,NB), pi_iter(NB)) * p_GT(GT)
        - PtData(GT,'Costo1') * u_GT(GT)
        - PtData(GT,'CI1') * (p_GT(GT) - PtData(GT,'Pmin') * u_GT(GT));

* [F3] Acoplamiento binario correcto: p=0 si u=0, p en [Pmin,Pmax] si u=1
eq_min_GT(GT)..
    p_GT(GT) =g= PtData(GT,'Pmin') * u_GT(GT);

eq_max_GT(GT)..
    p_GT(GT) =l= PtData(GT,'Pmax') * u_GT(GT);

eq_obj_GH(GH)..
    z_slave_GH_ind(GH) =e=
          sum(NB$PhBus(GH,NB), pi_iter(NB)) * p_GH(GH)
        - PhData(GH,'CI') * p_GH(GH);

eq_min_GH(GH)..
    p_GH(GH) =g= PhData(GH,'Pmin');

eq_max_GH(GH)..
    p_GH(GH) =l= PhData(GH,'Pmax');

Model Slave_GT /eq_obj_GT_total, eq_obj_GT, eq_min_GT, eq_max_GT/;
Model Slave_GH /eq_obj_GH_total, eq_obj_GH, eq_min_GH, eq_max_GH/;

Slave_GT.optcr = 0;
Slave_GH.optcr = 0;

*==============================================================================*
*------------------------------ SLAVE: RED DC --------------------------------*
*==============================================================================*
VARIABLE
       z_slave_net       "N(pi): congestion de la red (Eq. 4c)"
       f_dc_linea(ENL)   "Flujo DC en linea [MW]"
       ang(NB)           "Angulo de voltaje nodal [rad]"
       inj_fromto(NB)    "Inyeccion neta nodal [MW]";

EQUATIONS
       eq_obj_net
       eq_flujo_dc(ENL)
       eq_bal_net(NB);

eq_obj_net..
    z_slave_net =e= sum(NB, pi_iter(NB) * inj_fromto(NB));

* [F1] Uso de NBB como alias de NB para evitar conflicto con j del DAT
eq_flujo_dc(ENL)..
    f_dc_linea(ENL) =e= sum((NB,NBB)$Fbus(ENL,NB,NBB),
                            (ang(NB) - ang(NBB)) / FData(ENL,'X0'));

* Balance nodal: flujos salientes - entrantes + gen - demanda
* p_GT y p_GH estaran fijados con .fx antes de este solve (ver [F4])
eq_bal_net(NB)..
    inj_fromto(NB) =e=
          sum((ENL,NBB)$Fbus(ENL,NB,NBB), f_dc_linea(ENL))
        - sum((ENL,NBB)$Fbus(ENL,NBB,NB), f_dc_linea(ENL))
        + sum(GT$PtBus(GT,NB), p_GT(GT))
        + sum(GH$PhBus(GH,NB), p_GH(GH))
        - demanda(NB,'Pc');

* Limites de flujo por linea
f_dc_linea.up(ENL) =  FData(ENL,'Pmax');
f_dc_linea.lo(ENL) = -FData(ENL,'Pmax');

* Nodo de referencia angular (bus slack = B001)
ang.fx('B001') = 0;
ang.lo(NB)     = -PIVAL;
ang.up(NB)     =  PIVAL;

Model Slave_Net /eq_obj_net, eq_flujo_dc, eq_bal_net/;
Slave_Net.optcr = 0;

*==============================================================================*
*---------------------- PROGRAMA MAESTRO Y PROYECCION ------------------------*
*==============================================================================*
VARIABLES
       Z_dual        "L_hat(pi,k): aproximacion superior de L(pi)"
       theta_GT(GT)  "Cota inferior beneficio termico GT"
       theta_GH(GH)  "Cota inferior beneficio hidro GH"
       theta_net     "Cota inferior congestion de red"
       pi_new(NB)    "Nuevo precio candidato"
       dist          "Distancia cuadratica al centro de estabilidad";

EQUATIONS
       eq_obj_master
       eq_cut_GT(GT,kk)
       eq_cut_GH(GH,kk)
       eq_cut_net(kk)
       eq_min_dist
       eq_nivel_rst;

* Eq. 12a: L_hat = D*pi - sum(theta_GT) - sum(theta_GH) - theta_net
eq_obj_master..
    Z_dual =e= sum(NB, pi_new(NB) * demanda(NB,'Pc'))
             - sum(GT, theta_GT(GT))
             - sum(GH, theta_GH(GH))
             - theta_net;

* [F8] Filtro $(activo(kk)): solo genera cortes para iteraciones ya realizadas
eq_cut_GT(GT,kk)$(activo(kk))..
    theta_GT(GT) =g= hist_profit_GT(GT,kk)
                   + sum(NB$PtBus(GT,NB),
                         hist_p_GT(GT,kk) * (pi_new(NB) - hist_pi(NB,kk)));

eq_cut_GH(GH,kk)$(activo(kk))..
    theta_GH(GH) =g= hist_profit_GH(GH,kk)
                   + sum(NB$PhBus(GH,NB),
                         hist_p_GH(GH,kk) * (pi_new(NB) - hist_pi(NB,kk)));

* [F7] Signo negativo en gradiente: N concava => -N(pi') >= -N(pi_k) - s^T(pi'-pi_k)
eq_cut_net(kk)$(activo(kk))..
    theta_net =g= -hist_z_net(kk)
                - sum(NB, hist_dual_net(NB,kk) * (pi_new(NB) - hist_pi(NB,kk)));

eq_min_dist..
    dist =e= sum(NB, sqr(pi_new(NB) - pi_stab(NB)));

eq_nivel_rst..
    Z_dual =g= Nivel;

Model Maestro    /eq_obj_master, eq_cut_GT, eq_cut_GH, eq_cut_net/;
Model Proyeccion /eq_obj_master, eq_cut_GT, eq_cut_GH, eq_cut_net,
                  eq_min_dist, eq_nivel_rst/;

Maestro.optcr    = 0;
Proyeccion.optcr = 0;

* Caja Q inicial: paper usa +/- 300 USD/MWh
pi_new.lo(NB) = 0;
pi_new.up(NB) = 300;

* Demanda total del sistema
total_dem = sum(NB, demanda(NB,'Pc'));
Display 'Demanda total del sistema:', total_dem;

Scalar k_count "Contador de iteraciones" /0/;

*==============================================================================*
*========================== LOOP PRINCIPAL ====================================*
*==============================================================================*
Loop(k$(ord(k) <= 50 and Terminado = 0 and flag_err = 0),

    k_count = k_count + 1;

*--- PASO 1: Slave Termico (MIP) ----------------------------------------------
    Solve Slave_GT using mip maximizing z_slave_GT;

*   [F10] Display de alerta en lugar de abort — no corta ejecucion ni CSVs
    if(Slave_GT.modelstat > 2,
        Display 'ALERTA: Slave_GT no optimo en iteracion:', k_count;
        Display 'Slave_GT modelstat:', Slave_GT.modelstat;
        flag_err = 1;
    );

*--- PASO 2: Slave Hidro (LP) -------------------------------------------------
    if(flag_err = 0,
        Solve Slave_GH using lp maximizing z_slave_GH;

        if(Slave_GH.modelstat > 2,
            Display 'ALERTA: Slave_GH no optimo en iteracion:', k_count;
            Display 'Slave_GH modelstat:', Slave_GH.modelstat;
            flag_err = 1;
        );
    );

*--- PASO 3: Slave Red (LP) ---------------------------------------------------
*   [F4] Fijar p_GT y p_GH antes del solve para que Slave_Net no los
*        re-optimice. Son datos de entrada para la red, no variables de red.
    if(flag_err = 0,
        p_GT.fx(GT) = p_GT.l(GT);
        p_GH.fx(GH) = p_GH.l(GH);

        Solve Slave_Net using lp minimizing z_slave_net;

*       Restaurar limites para la siguiente iteracion
        p_GT.lo(GT) = 0;
        p_GT.up(GT) = PtData(GT,'Pmax');
        p_GH.lo(GH) = 0;
        p_GH.up(GH) = PhData(GH,'Pmax');

        if(Slave_Net.modelstat > 2,
            Display 'ALERTA: Slave_Net no optimo en iteracion:', k_count;
            Display 'Slave_Net modelstat:', Slave_Net.modelstat;
            flag_err = 1;
        );
    );

*--- PASO 4: Registrar corte k -----------------------------------------------
    if(flag_err = 0,
        activo(k) = 1;

        hist_pi(NB,k)        = pi_iter(NB);
        hist_profit_GT(GT,k) = z_slave_GT_ind.l(GT);
        hist_p_GT(GT,k)      = p_GT.l(GT);
        hist_profit_GH(GH,k) = z_slave_GH_ind.l(GH);
        hist_p_GH(GH,k)      = p_GH.l(GH);
        hist_z_net(k)        = z_slave_net.l;

*       [F5] Supergradiente correcto = inyeccion neta optima (teorema del envelope)
*            NO usar eq_bal_net.m(NB) que por KKT del LP equivale a pi_iter
        hist_dual_net(NB,k)  = inj_fromto.l(NB);

        hist_gen_GT(GT,k) = p_GT.l(GT);
        hist_gen_GH(GH,k) = p_GH.l(GH);
        hist_u_GT(GT,k)   = u_GT.l(GT);

*--- PASO 5: Evaluar L(pi_k) — candidato a Lower Bound -----------------------
*       [F6] N(pi) se SUMA: L(pi) = D*pi - Pi_GT - Pi_GH + N(pi)
        L_actual = sum(NB, pi_iter(NB) * demanda(NB,'Pc'))
                 - sum(GT, z_slave_GT_ind.l(GT))
                 - sum(GH, z_slave_GH_ind.l(GH))
                 + z_slave_net.l;

*--- PASO 6: Actualizar LB y centro de estabilidad ---------------------------
        if(L_actual > LB,
            LB = L_actual;
            pi_stab(NB) = pi_iter(NB);
        );

*--- PASO 7: Resolver Maestro (LP) — nuevo UB --------------------------------
        Solve Maestro using lp maximizing Z_dual;

        if(Maestro.modelstat > 2,
            Display 'ALERTA: Maestro no optimo en iteracion:', k_count;
            Display 'Maestro modelstat:', Maestro.modelstat;
            flag_err = 1;
        else
            UB = Z_dual.l;

*--- PASO 8: Gap y criterio de parada ----------------------------------------
*           [F9] Denominador abs(UB) segun Eq.(8) del paper
            Gap = (UB - LB) / (abs(UB) + 1e-6);

            hist_LB(k)  = LB;
            hist_UB(k)  = UB;
            hist_Gap(k) = Gap;

            diferencia = sum(GT, p_GT.l(GT)) + sum(GH, p_GH.l(GH)) - total_dem;

            if(Gap < epsilon,
                Terminado = 1;
            );

*--- PASO 9: Proyeccion — nuevo precio candidato pi_{k+1} -------------------
            if(Terminado = 0,
                Nivel = alpha * UB + (1 - alpha) * LB;

*               [F14] Encogimiento de la caja Q (paper Seccion III-C)
                if(k_count = 10 or k_count = 20 or k_count = 30,
                    pi_new.lo(NB) = max(0,    pi_iter(NB) - 25);
                    pi_new.up(NB) = min(2000, pi_iter(NB) + 25);
                );

                Solve Proyeccion using qcp minimizing dist;

                if(Proyeccion.modelstat > 2,
                    Display 'ALERTA: Proyeccion no optima en iteracion:', k_count;
                    Display 'Proyeccion modelstat:', Proyeccion.modelstat;
                    flag_err = 1;
                else
                    pi_iter(NB) = pi_new.l(NB);
                );
            );
        );
    );

);
* ============================= FIN DEL LOOP ==================================

*==============================================================================*
*------------------------------ REPORTES FINALES ------------------------------*
*==============================================================================*
* [F13] Todos los Display fuera del loop
Display 'Numero de iteraciones realizadas:', k_count;
Display 'Lower Bound final:', LB;
Display 'Upper Bound final:', UB;
Display 'Gap final:', Gap;
Display 'Convergencia (1=si 0=no):', Terminado;
Display 'Flag de error (1=hubo error):', flag_err;
Display hist_Gap;
Display hist_LB;
Display hist_UB;
Display hist_profit_GT;
Display hist_p_GT;
Display hist_profit_GH;
Display hist_p_GH;
Display hist_u_GT;

*==============================================================================*
*------------------------ EXPORTACION DE RESULTADOS --------------------------*
*==============================================================================*

* [F11] Todos los archivos tienen PUTCLOSE al final

*--- Archivo 1: Generacion final ---------------------------------------------*
FILE f_gen '\Resultados_Generacion.csv\';
f_gen.pc = 5;
f_gen.pw = 10000;
PUT f_gen;
PUT 'Tipo,Generador,Nodo,Potencia_MW,Estado_Commitment,Precio_Nodal_USD_MWh' /;
LOOP((GT,NB)$PtBus(GT,NB),
    PUT 'Termica', GT.tl, NB.tl,
        p_GT.l(GT):0:4,
        u_GT.l(GT):0:0,
        pi_iter(NB):0:4 /;
);
LOOP((GH,NB)$PhBus(GH,NB),
    PUT 'Hidraulica', GH.tl, NB.tl,
        p_GH.l(GH):0:4,
        '1',
        pi_iter(NB):0:4 /;
);
PUTCLOSE f_gen;

*--- Archivo 2: Precios nodales finales --------------------------------------*
FILE f_price '\Resultados_Precios.csv\';
f_price.pc = 5;
f_price.pw = 10000;
PUT f_price;
PUT 'Nodo,Precio_Final_USD_MWh' /;
LOOP(NB,
    PUT NB.tl, pi_iter(NB):0:6 /;
);
PUTCLOSE f_price;

*--- Archivo 3: Historial de convergencia ------------------------------------*
FILE f_conv '\Historial_Convergencia.csv\';
f_conv.pc = 5;
f_conv.pw = 10000;
PUT f_conv;
PUT 'Iteracion,Lower_Bound,Upper_Bound,Gap_Relativo' /;
LOOP(k$(activo(k)),
    PUT k.tl,
        hist_LB(k):0:6,
        hist_UB(k):0:6,
        hist_Gap(k):0:8 /;
);
PUTCLOSE f_conv;

*--- Archivo 4: Historial de generacion por iteracion ------------------------*
FILE f_geniter '\Historial_Generacion.csv\';
f_geniter.pc = 5;
f_geniter.pw = 10000;
PUT f_geniter;
PUT 'Iteracion,Generador,Tipo,Potencia_MW,Commitment' /;
LOOP(k$(activo(k)),
    LOOP(GT,
        PUT k.tl, GT.tl, 'Termica',
            hist_gen_GT(GT,k):0:4,
            hist_u_GT(GT,k):0:0 /;
    );
    LOOP(GH,
        PUT k.tl, GH.tl, 'Hidraulica',
            hist_gen_GH(GH,k):0:4,
            '1' /;
    );
);
PUTCLOSE f_geniter;

*--- Archivo 5: Precios nodales por iteracion --------------------------------*
FILE f_pi '\Precios_Nodales_por_Iteracion.csv\';
f_pi.pc = 5;
f_pi.pw = 10000;
PUT f_pi;
PUT 'Iteracion,Nodo,Precio_USD_MWh' /;
LOOP(k$(activo(k)),
    LOOP(NB,
        PUT k.tl, NB.tl, hist_pi(NB,k):0:6 /;
    );
);
PUTCLOSE f_pi;
