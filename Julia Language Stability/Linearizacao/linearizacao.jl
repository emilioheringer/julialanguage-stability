
function linearizacao(Y_bus)
    ##Criacao dos vetores que será perturbado
    x=[delta' dWpu' Eql' Efd'];
    z=[teta' V'];
    u=[Prp' Vref'];
    xzu0=[x z u];
    tol=1*10^-5;
    ns=10;
    ##inicializacao de matrizes
    mp=zeros(4*ngen+2*nbus+ns,size(xzu0)[2]);
    xzu=zeros(size(xzu0));
    F1_0=zeros(ngen,1);
    F2_0=zeros(ngen,1);
    F3_0=zeros(ngen,1);
    F4_0=zeros(ngen,1);
    F1_P=zeros(ngen,1);
    F2_P=zeros(ngen,1);
    F3_P=zeros(ngen,1);
    F4_P=zeros(ngen,1);
    G3_0=zeros(nbus,1);
    G4_0=zeros(nbus,1);
    G3_P=zeros(nbus,1);
    G4_P=zeros(nbus,1);
    H0=zeros(ns,1);
    HP=zeros(ns,1);

    ##Calculo da potencia injetada
    Th = xzu0[1+4*ngen:4*ngen+nbus];
    Vt = xzu0[1+4*ngen+nbus:4*ngen+2*nbus];

    #Adicinando cargas na ybus

    Vrec=Vt.*cos.(Th)+Vt.*sin.(Th)im;
    IBUS=Y_bus*Vrec;
    Sinj=Vrec.*conj(IBUS);
    Pinj=real(Sinj);
    Qinj=imag(Sinj);

    global e=0;
    ##Equacoes de estado
    for i=1:nbus
        if bus_data[i,2] == 2 || bus_data[i,2] == 1
            global e=e+1;
            dt,vel,ed,ef=estado(xzu0,e);
            F1_0[e,1]=dt;
            F2_0[e,1]=vel;
            F3_0[e,1]=ed;
            F4_0[e,1]=ef;
        end

    end

    #Equacoes Algebricas
    global e=0;
    for i=1:nbus
        if bus_data[i,2] == 2 || bus_data[i,2] == 1   ##Possui geradores
            global e=e+1;
            pg,qg=algebricas_gerador(xzu0,e);
            G3_0[i,1] = pg-0*PD[i]-Pinj[i];
            G4_0[i,1] = qg-0*QD[i]-Qinj[i];
        else                                          ##Não Possui Gerador
            G3_0[i,1] = -0*PD[i]-Pinj[i];
            G4_0[i,1] = -0*QD[i]-Qinj[i];
        end
    end

    for i=1:ns
       dY= saida(xzu0,i);
       H0[i]=dY;
    end
    #Vetor com os calculos inicias
    F0=[F1_0;F2_0;F3_0;F4_0;];
    G0=[G3_0;G4_0];
    FG0=[F0;G0;H0];
 
    ##Criando as perturbacoes
   
    #Percorre todas as variaveis
    for j=1:size(xzu0)[2]

        #Restaura as variaveis para o valor inicial
        for k=1:size(xzu0)[2]
            xzu[k]=xzu0[k];
        end
        
        #Pertubacao
        xzu[j]=xzu0[j]+tol;

        #Calculo da potencia injetada
        Th = xzu[1+4*ngen:4*ngen+nbus];
        Vt = xzu[1+4*ngen+nbus:4*ngen+2*nbus];

        Vrec=Vt.*cos.(Th)+Vt.*sin.(Th)im;
        IBUS=Y_bus*Vrec;
        Sinj=Vrec.*conj(IBUS);
        Pinj=real(Sinj);
        Qinj=imag(Sinj);
        global e=0;
        #Estado
        for i=1:nbus
            if bus_data[i,2] == 2 || bus_data[i,2] == 1
                global e=e+1;
                dt,vel,ed,ef=estado(xzu,e);
                F1_P[e,1]=dt;
                F2_P[e,1]=vel;
                F3_P[e,1]=ed;
                F4_P[e,1]=ef;   
            end
         
        end
        global e=0;
        #Algebricas
        for i=1:nbus           
            if bus_data[i,2] == 2 || bus_data[i,2] == 1   ##Possui geradores
                global e=e+1;
                pg,qg=algebricas_gerador(xzu,e);
                G3_P[i,1] = pg-0*PD[i]-Pinj[i];
                G4_P[i,1] = qg-0*QD[i]-Qinj[i];
            else                                          ##Não Possui Gerador
                G3_P[i,1] = -0*PD[i]-Pinj[i];
                G4_P[i,1] = -0*QD[i]-Qinj[i];
            end       
         
        end
        
        for i=1:ns
            dY= saida(xzu,i);
            HP[i]=dY;
        end
        #Vetores apos perturbacao
        FP=[F1_P;F2_P;F3_P;F4_P];
        GP=[G3_P;G4_P];
        FGP=[FP;GP;HP];

        #Formacao da matriz de perturbacao
        dfg = (FGP-FG0)/tol;
        mp[:,j]=dfg;
    end
    return mp;
end
    


    
