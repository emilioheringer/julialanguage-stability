## DADOS DO SISTEMA
##DADOS BASE
s_base=100;
ngen=3;

##TIPOS DE BARRAS

##DADOS DA BARRA
##TIPO 0 -> PQ    TIPO 1 -> PV   TIPO 2 -> VTETA
##               NUM        TYPE            V           TETA         PG            QG         PD           QD           GSH          BSHJ
bus_data = [     1            2            1.04         0.000       71.95         24.07        0            0            0            0;
                 2            1            1.025        9.669      163.00         76.61        0            0            0            0;   
                 3            1            1.025        4.771       85.00         14.00        0            0            0            0;
                 4            0            0.987       -2.407       00.00         0         0            0            0            0;
                 5            0            0.958       -4.350          0          0          125           50            0            0;                            
                 6            0            0.975       -4.017          0          0           90           30            0            0;                            
                 7            0            0.996        3.799       00.00         00.00        0            0            0            0;
                 8            0            0.986        0.622       00.00         00.00      100           35            0            0;
                 9            0            1.003        1.926       00.00         00.00        0            0            0            0;
           ];


##DADOS DAS LINHAS
##             FROM     TO          R             X           BSH(TOTAL)  TAP
line_data =  [  1        4        0.0000        0.0576       0.0000        0;
                4        5        0.0100        0.0850       0.1760        0;
                4        6        0.0170        0.0920       0.1580        0;
                5        7        0.0320        0.1610       0.3060        0;
                6        9        0.0390        0.1700       0.3580        0;
                2        7        0.0000        0.0625       0.0000        0;
                7        8        0.0085        0.0720       0.1490        0;
                8        9        0.0119        0.1008       0.2090        0;
                3        9        0.0000        0.0586       0.0000        0
             ];

##Reatancia transitoria de eixo direto dos geradores 1 ate n;            
##                     Bus   MODEL  BASE      H     Xdl      D       Xd       Xq      Xql      Tdol      R
global generator_data =[1      2     100   23.64   0.0608   0      0.1460   0.0969   0.0969   8.96      0;
                        2      1     100   6.4     0.1198    0      0.8958   0.8645   0.1969   6.00      0;
                        3      1     100   3.01    0.1813    0      1.3125   1.2578   0.2500   5.89      0];
            #Bus  Ka   Ta
global exc = [1   20   0.2 ;
              2   20   0.2 ;
              3   20   0.2 ];



global saidas        = [7     1    0    0;
                        7     2    0    0;
                        7     3    0    0;
];

global ns=size(saidas,1);