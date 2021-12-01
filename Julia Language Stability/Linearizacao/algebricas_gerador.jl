function algebricas_gerador(xzu,i)

    #Dados Geradores
    Xdl=generator_data[i,5];
    Xq=generator_data[i,8];
    Xdl=generator_data[i,5];


    ##separacao de variaveis
    #Estado
    Delta = xzu[1:ngen];
    eql   = xzu[1+2*ngen:3*ngen];


    #Algebrica
    Th = xzu[1+4*ngen:4*ngen+nbus];
    Vt = xzu[1+4*ngen+nbus:4*ngen+2*nbus];


    Eqli=eql[i];
    Vi=Vt[i];
    Thi=Th[i];
    Deltai=Delta[i];
    phi=Deltai-Thi;
    

    Pg = Eqli*Vi*sin(phi)/Xdl + 0.5*(1/Xq - 1/Xdl)*(Vi^2)*sin(2*phi);
    Qg = Eqli*Vi*cos(phi)/Xdl - (Vi^2)/Xdl - 0.5*(1/Xq - 1/Xdl)*(Vi^2)*(1-cos(2*phi));
    return Pg,Qg;
end