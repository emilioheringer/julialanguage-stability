function saida(xzu,out)
    
    type = saidas[out,1];
    dy1 = saidas[out,2];
    dy2 = saidas[out,3];
    dy3 = saidas[out,4];
    Th = xzu[1+4*ngen:4*ngen+nbus];
    Vt = xzu[1+4*ngen+nbus:4*ngen+2*nbus];
    if type==1  #tensao
        dY= V[dy1];
    elseif type==7     #velocidade gerador 1
        geni= dy1;
        dY = dWpu[geni];
    elseif type==6     #Dv = Vk-Vm
        dY = Vt[dy1]-Vt[dy2];
    elseif type==2     #Dt = tk - tm
        geni= dy1;
        dY = dWpu[geni];
    elseif type==3     #Dwk - dWm
        gen1 = dy1;
        gen2 = dy2;
        dY = dWpu[gen1] - dWpu[gen2];
    elseif type==8     #dY = Pkm
        s=dy1;
        fr=line_data[s,1];
        to=line_data[s,2];
        Vk=Vt[fr];
        tk=Th[fr];
        R = line_data[s,3];
        Vm=Vt[to];
        tm=Th[to];
        X = line_data[s,3];
        y=inv(R+im*X);
        gkm = real(y);
        bsh = line_data[s,5];
        bkm = imag(y);
        akm = line_data[s,6];
        Pkm = (akm^2)*(Vk^2)*gkm - akm*Vk*Vm*gkm*cos(tk-tm) - akm*Vk*Vm*bkm*sin(tk-tm);
        Pmk = Vm^2*gkm - akm*Vk*Vm*gkm*cos(tk-tm) + akm*Vk*Vm*bkm*sin(tk-tm);
        dY = Pkm;
    end
    return dY;
end
