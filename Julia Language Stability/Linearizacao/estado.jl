function estado(xzu,i)

  #Dados Geradores
  H=generator_data[i,4];
  Xdl=generator_data[i,5];
  D=generator_data[i,6];
  Xd=generator_data[i,7];
  Xq=generator_data[i,8];
  Xdl=generator_data[i,5];
  Td0l=generator_data[i,10];
  KA=exc[i,2];
  TA=exc[i,3];


  ##separacao de variaveis
  #Estado
  Delta = xzu[1:ngen];
  dwpu  = xzu[1+1*ngen:2*ngen];
  eql   = xzu[1+2*ngen:3*ngen];
  efd   = xzu[1+3*ngen:4*ngen];

  #Algebrica
  Th = xzu[1+4*ngen:4*ngen+nbus];
  Vt = xzu[1+4*ngen+nbus:4*ngen+2*nbus];
   
  #Entrada
  pm = xzu[1+4*ngen+2*nbus:5*ngen+2*nbus];
  vref = xzu[1+5*ngen+2*nbus:6*ngen+2*nbus];


  ws=377;
  Eqli=eql[i];
  Efdi = efd[i]
  Pmi = pm[i];
  Vi=Vt[i];
  Thi=Th[i];
  Deltai=Delta[i];
  phi=Deltai-Thi;
  Vrefi=vref[i];
  
    
  Pgi = Eqli*Vi*sin(phi)/Xdl + 0.5*(1/Xq - 1/Xdl)*Vi^2*sin(2*phi);
  Idi  = (Eqli - Vi*cos(phi))/Xdl;

  #Equacoes
    
  delta=ws*dwpu[i,1];
  velocidade = (Pmi- Pgi - D*dwpu[i])/(2*H)
  ed= -Eqli/Td0l - (Xd-Xdl)*Idi/Td0l + Efdi/Td0l;
  ef= -Efdi/TA + KA/TA*(Vrefi - Vi);
    
  return delta,velocidade,ed,ef;
    
end