delta=zeros(ngen);
Id=zeros(ngen);
Iq=zeros(ngen);
Vd=zeros(ngen);
Vq=zeros(ngen);
Eql=zeros(ngen);
Edl=zeros(ngen);
Efd=zeros(ngen);
Vr=zeros(ngen);
Rf=zeros(ngen);
Vref=zeros(ngen);
Prp=zeros(ngen);
Qrp=zeros(ngen);
tetag=zeros(ngen);
Vg=zeros(ngen);
dWpu=zeros(ngen);
Vn=zeros(nbus-ngen);
tetan=zeros(nbus-ngen);
Ig = Complex.(zeros(ngen));
global e=0;
for i=1:nbus
    
        ##Angulo interno
        if bus_data[i,2] == 2 || bus_data[i,2] == 1
                global e=e+1;
                ##Corrente de geracao
                Ig[e] = conj(((PG[i,1]) + im *(QG[i,1])) / ((V[i,1] * cos(teta[i,1])) + im*(V[i,1] * sin(teta[i,1]))));
                
                delta[e]=angle(V[i,1]*(cos(teta[i,1]) + im*sin(teta[i,1])) + im*(generator_data[e,8])*(Ig[e]));
                ##Corrente na referencia do gerador
                Id[e]= abs(Ig[e])*sin(delta[e]-angle(Ig[e]));
                Iq[e]= abs(Ig[e])*cos(delta[e]-angle(Ig[e]));
                #Tensoes na referencia do gerador
                Vd[e]= V[i,1]*sin(delta[e]-teta[i,1]);
                Vq[e]= V[i,1]*cos(delta[e]-teta[i,1]);
        
                #Tensao transitoria
                Eql[e] = Vq[e] + Id[e]*generator_data[e,5];
                Edl[e] = (generator_data[e,8]-generator_data[e,9])*Iq[e];
                #Tensao de campo
                Efd[e]= Eql[e] + (generator_data[e,7]-generator_data[e,5])*Id[e];
             
                #Tensao de referencia
                Vref[e] = V[i,1] + Efd[e]/exc[e,2];

                ##Potencia em regime permanente
                Prp[e] = bus_data[i,5]/s_base;# Vd[e]*Iq[e] + Vq[e]*Iq[e];
                Qrp[e] = bus_data[i,6]/s_base;#Vq[e]*Id[e] - Vd[e]*Iq[e];

        end
       # Sl[i] = PD[i] + im*QD[i];
        #Yl[i]=conj(Sl[i])/V[i]^2;
    
end

#Adicionando as cargas na ybus
#Yl=Diagonal(Yl);
#Y_bus=Y_bus+Yl;
        
        