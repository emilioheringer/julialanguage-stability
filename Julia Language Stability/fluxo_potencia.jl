
##Fluxo de POTENCIA
global V
##CONSTANTES
toler_P = 1e-20;           # Tolerância para P
toler_Q = 1e-20;           # Tolerância para q
maxit = 20;               # Número máximo de iterações

##ATRIBUIÇÃO DE CO@NSTANTES PARA OS DADOS DO SISTEMA

num = bus_data[:,1];     ##NÚMERO DA BARRA
type = bus_data[:,2];     ##TÍPO DE BARRA
vesp = bus_data[:,3];    ##TENSÕES INICIAIS ESPERADAS PARA GERAÇÃO E CARGA

teta = (pi / 180) * bus_data[:,4];  ##ANGULOS EM RADIANOS
PG = bus_data[:,5] / s_base;      ##POTÊNCIA ATIVA GERADA
QG = bus_data[:,6] / s_base;      ##POTÊNCIA REATIVA GERADA
PD = bus_data[:,7] / s_base;      ##POTÊNCIA ATIVA DEMANDADA
QD = bus_data[:,8] / s_base;      ##POTÊNCIA REATIVA DEMANDADA
GSH = bus_data[:,9] / s_base;     ##SHUNT DE BARRA -> RESISTIVO
BSH = bus_data[:,10] / s_base;    ##SHUNT DE BARRA -> REATIVO


from_bus = line_data[:,1];      ##BARRA DE
to_bus = line_data[:,2];        ##BARRA PARA 
R = line_data[:,3];             ##RESISTÊNCIA
X = line_data[:,4];             ##REATÂNCIA
BSHKM = line_data[:,5] / 2;       ##SUSCENPTÂNCIA SHUNT EM CADA EXTREMIDADE
TAP = line_data[:,6];
indx= findall(x -> x == 0, TAP);
TAP[indx] .= 1;
TAP = 1 ./TAP;




##INICIALIZAÇÃO DE VARIÁVEIS
global nbus = size(num, 1);     ##NÚMERO DE BARRAS
nlin = size(X, 1);       ##NÚMERO DE CIRCUITOS
YKM = zeros(nlin);
YKM = complex.(YKM);
for i = 1:nlin
    YKM[i] = 1 / (R[i] + im * X[i]);                ##SUSCENPTÂNCIA DE CADA CIRCUITO []
end

Y_bus = zeros(nbus, nbus);       ##INICIALIZAÇÃO DA MATRIZ DE ADMITÂNCIA

p_esp = PG - PD;                  ##INJEÇÃO DE POTÊNCIA ATIVA LÍQUIDA
q_esp = QG - QD;                  ##INJEÇÃO DE POTÊNCIA REATIVA LÍQUIDA

p_cal = zeros(nbus, 1);              ##INICIALIZAÇÃO ARRAY INJEÇÃO ATIVA CALCULADA
q_cal = zeros(nbus, 1);              ##INICIALIZAÇÃO ARRAY INJEÇÃO REATIVA CALCULADA
dp = zeros(nbus, 1);                 ##INICIALIZAÇÃO ARRAY DE RESÍDUOS DE ATIVO
dq = zeros(nbus, 1);                 ##INICIALIZAÇÃO ARRAY DE RESÍDUOS DE REATIVO

V = vesp;

pos_PV = findall(x -> x == 1, type);   ##POSIÇÃO DAS BARRAS PV
pos_VT = findall(x -> x == 2, type);   ##POSIÇÃO DAS BARRAS VT
pos_PQ = findall(x -> x == 0, type);   ##POSIÇÃO DAS BARRAS PQ

#Q_err = zeros(maxit, 2);
#P_err = zeros(maxit, 2);
##TRANSFORMAÇÃO DE VETORES EM MATRIZES COLUNAS
GSH = reshape(GSH, length(GSH), 1);
BSH = reshape(BSH, length(BSH), 1);
BSHKM = reshape(BSHKM, length(BSHKM), 1);
YKM = reshape(YKM, length(YKM), 1);
V = reshape(V, length(V), 1);
teta = reshape(teta, length(teta), 1);
##COVENTENDO TIPOS DE DADOS

from_bus = Int.(from_bus);
to_bus = Int.(to_bus);
Y_bus = Complex.(Y_bus);


##FORMAÇÃO DA MATRIX YBARRA

for i = 1:nbus
    Y_bus[i,i] = Y_bus[i,i] + (GSH[i,1] + im * BSH[i,1]);
end


for a = 1:nlin
    i = from_bus[a];
    j = to_bus[a];
    akm = TAP[a];
    Y_bus[i,i] = Y_bus[i,i] + im*BSHKM[a,1] + YKM[a,1]*akm^2;
    Y_bus[i,j] = Y_bus[i,j] - akm*YKM[a,1];
    Y_bus[j,i] = Y_bus[j,i] - akm*YKM[a,1];
    Y_bus[j,j] = Y_bus[j,j] + YKM[a,1] + im*BSHKM[a,1];
end

G = real(Y_bus);
B = imag(Y_bus);

##SUBSISTEMA 1
it = 1;           ##COLOCANDO ITERAÇÕES NO VALOR 1;

P_error = 0;    ##ERRO ATIVO
Q_error = 0;    ##ERRO REATIVO
max_error = max(P_error, Q_error);     ##MAIOR VALOR DE ERRO

##INJEÇÃO DE POTÊNCIA NAS BARRAS
Pcalc = zeros(nbus);
Qcalc = zeros(nbus);
dP = zeros(nbus);
dQ = zeros(nbus);
global Bbus, Gbus;
Bbus=B;
Gbus=G;
for i = 1:nbus
    Vk = V[i,1];
    tk = teta[i,1];

    for j = 1:nbus
        Vm = V[j,1];
        tm = teta[j,1];

        Gkm = G[i,j];
        Bkm = B[i,j];

        Pcalc[i,1] = Pcalc[i,1] + Vk * Vm * (Gkm * cos(tk - tm) + Bkm * sin(tk - tm));
        Qcalc[i,1] = Qcalc[i,1] + Vk * Vm * (Gkm * sin(tk - tm) - Bkm * cos(tk - tm));
    
    end

    dP[i] = p_esp[i] - Pcalc[i,1];
    dQ[i] = q_esp[i] - Qcalc[i,1];
end

##ANULA RESÍDUOS ATIVOS PARA VTETA E RESIDOS REATIVOS PARA VTETA E PV
dP[pos_VT] .= 0;
dQ[pos_VT] .= 0;
dQ[pos_PV] .= 0;

P_error = maximum(abs.(dP));
Q_error = maximum(abs.(dQ));
max_error = max(P_error, Q_error);


##PROCESSO ITERATIVO

while it <= maxit && ((P_error > toler_P) || (Q_error > toler_Q))
    global it += 1;
    H = zeros(nbus, nbus);
    N = zeros(nbus, nbus);
    M = zeros(nbus, nbus);
    L = zeros(nbus, nbus);
   ##MATRIX JACOBIANA
    
   
    for i = 1:nbus
        Vk = V[i,1];
        tk = teta[i,1];
        
        for j = 1:nbus
            Vm = V[j,1];
            tm = teta[j,1];
            
            Gkm = G[i,j];
            Bkm = B[i,j];

            if i != j       ##ELEMENTOS FORA DA DIAGONAL

                H[i,j] = Vk * Vm * (Gkm * sin(tk - tm) - Bkm * cos(tk - tm));
                N[i,j] = Vk * (Gkm * cos(tk - tm) + Bkm * sin(tk - tm));
                M[i,j] = -Vk * Vm * (Gkm * cos(tk - tm) + Bkm * sin(tk - tm))
                L[i,j] = Vk * (Gkm * sin(tk - tm) - Bkm * cos(tk - tm));
                
            else          ##ELEMENTOS DA DIAGONAL
                H[i,i] = -Vk^2 * Bkm - Qcalc[i];
                N[i,i] = (Pcalc[i] + Vk^2 * Gkm) / Vk;
                M[i,i] = -Vk^2 * Gkm + Pcalc[i];
                L[i,i] = (Qcalc[i] - Vk^2 * Bkm) / Vk; 
            end

        end
           

        ##ELIMINAÇÃO DE EQUAÇÕES
        if type[i] == 2   ##VTETA
            H[i,i] = 1e10;
            L[i,i] = 1e10;

        elseif type[i] == 1
            L[i,i] = 1e10;
        end
         

    end

    J = ([H N ; M L]);   
    J = sparse(J);   
    TV = [teta ; V];
    dPQ = [dP;dQ];  
    dPQ = reshape(dPQ, length(dPQ), 1);  
    dTV = J \ dPQ;
    TV  = TV + dTV;
   
   ##ATUALIZANDO TETA
    for k = 1:nbus
        teta[k,1] = TV[k];
    end

   ##ATUALIZANDO v
    m = 1;
    for k = nbus + 1:length(TV)
        V[m,1] = TV[k];
        m=m+1;
    end
    V[pos_VT,1] = bus_data[pos_VT,3];
    teta[pos_VT,1] = bus_data[pos_VT,4];
   ##INJEÇÃO DE POTÊNCIA NAS BARRAS 
   ##ANULA RESÍDUOS ATIVOS PARA VTETA E RESIDOS REATIVOS PARA VTETA E PV
   
   ##INJEÇÃO DE POTÊNCIA NAS BARRAS
    global Pcalc = zeros(nbus);
    global Qcalc = zeros(nbus);
    global dP = zeros(nbus);
    global dQ = zeros(nbus);

    for i = 1:nbus
        Vk = V[i,1];
        tk = teta[i,1];
    
        for j = 1:nbus
            Vm = V[j,1];
            tm = teta[j,1];
    
            Gkm = G[i,j];
            Bkm = B[i,j];
    
            Pcalc[i,1] = Pcalc[i,1] + Vk * Vm * (Gkm * cos(tk - tm) + Bkm * sin(tk - tm));
            Qcalc[i,1] = Qcalc[i,1] + Vk * Vm * (Gkm * sin(tk - tm) - Bkm * cos(tk - tm));
        end
    
        dP[i] = p_esp[i] - Pcalc[i,1];
        dQ[i] = q_esp[i] - Qcalc[i,1];
    end

    ##ANULA RESÍDUOS ATIVOS PARA VTETA E RESIDOS REATIVOS PARA VTETA E PV
    dP[pos_VT] .= 0;
    dQ[pos_VT] .= 0;
    dQ[pos_PV] .= 0;

    global P_error = maximum(abs.(dP));
    global Q_error = maximum(abs.(dQ));
    global max_error = max(P_error, Q_error);
    
end 

##SUBSISTEMA 2

##CALCULA Pk e Qk PARA VTETA
##CALCULA  Qk PARA VTETA

p_esp[pos_VT] .= 0;
q_esp[pos_VT] .= 0;
q_esp[pos_PV] .= 0;   

for i = 1:nbus
    Vk = V[i,1];
    tk = teta[i,1];

    for j = 1:nbus
        Vm = V[j,1];
        tm = teta[j,1];

        Gkm = G[i,j];
        Bkm = B[i,j];

        if type[i] == 1 || type[i] == 2  ##PV OU VTETA
            q_esp[i] = q_esp[i] + Vk * Vm * (Gkm * sin(tk - tm) - Bkm * cos(tk - tm)); 

            if type[i] == 2   ##VTETA
                p_esp[i] = p_esp[i] + Vk * Vm * (Gkm * cos(tk - tm) + Bkm * sin(tk - tm));
                
            end
        end

       
    end

end
q_esp = reshape(q_esp, length(q_esp), 1);
QD = reshape(QD, length(QD), 1);
p_esp = reshape(p_esp, length(p_esp), 1);
DP = reshape(PD, length(PD), 1);
PG = p_esp + PD;
QG = q_esp + QD;

##SHUNT DA BARRA
bus_sh = BSH .* (V .* V);


println("#####  Relatório de Iterações\n\n")
println("Máximo erro de potência ativa:   ",s_base * P_error)
println("Máximo erro de potência reativa:   ",s_base * Q_error)
println("Tolerância:   ",toler_P * s_base);
println("Total de iterações:   ",it);
println("\n\n####Dados do Sistema")
println("Total de carga ativa:       ",s_base * sum(PD))
println("Total de carga reativa:     ",s_base * sum(QD))
println("Total de geração ativa:     ",s_base * sum(PG))
println("Total de geração reativa:   ",s_base * sum(QG))
println("\n\n                        ##### Fluxo de potência #####")
power_Flow = [type V teta*180/pi PG QG PD QD bus_sh];

header=["Tipo" "Tensão(V)" "Ângulo(°)" "PG" "QG" "PD" "QD" "Bus Shunt"];
pretty_table(power_Flow,header, border_crayon = crayon"bold yellow", header_crayon = crayon"bold white")

